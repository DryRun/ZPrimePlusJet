import ROOT
from ROOT import *
from multiprocessing import Process
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array
import os

from DAZSLE.ZPrimePlusJet.xbb_config import analysis_parameters as params 
import DAZSLE.PhiBBPlusJet.analysis_config as config

from rhalphabet_builder import GetSF

def writeDataCard(boxes,workspace_path,sigs,bkgs,histograms,options,datacard_path):
	obsRate = {}
	for box in boxes:
		obsRate[box] = histograms['data_obs_%s'%box].Integral()
	nBkgd = len(bkgs)
	nSig = len(sigs)

	rates = {}
	lumiErrs = {}
	hbb125ptErrs = {}
	mcStatErrs = {}
	veffErrs = {}
	bbeffErrs = {}
	znormEWErrs = {}
	znormQErrs = {}
	wznormEWErrs = {}
	mutriggerErrs = {}
	muidErrs = {}
	muisoErrs = {}
	jesErrs = {}
	jerErrs = {}
	puErrs = {}
	for proc in sigs+bkgs:
		for box in boxes:
			process_name = "{}_{}".format(proc, box)
			error = array.array('d',[0.0])
			rate = histograms[process_name].IntegralAndError(1,histograms[process_name].GetNbinsX(),error)
			rates[process_name]  = rate
			lumiErrs[process_name] = 1.025
			if proces == "hbb":
				hbb125ptErrs[process_name] = 1.3
			else:
				hbb125ptErrs[process_name] = 1.0
			if proc=='wqq' or proc=='zqq' or 'hbb' in proc or 'Sbb' in proc:
				veffErrs[process_name] = 1.0+params[options.jet_type]["V_SF_ERR"]/params[options.jet_type]["V_SF"]
				if box=='pass':
					bbeffErrs[process_name] = 1.0+params[options.jet_type]["BB_SF_ERR"]/params[options.jet_type]["BB_SF"]
				else:
					ratePass = histograms['%s_%s'%(proc,'pass')].Integral()
					rateFail = histograms['%s_%s'%(proc,'fail')].Integral()
					if rateFail>0:
						bbeffErrs[process_name] = 1.0-params[options.jet_type]["BB_SF_ERR"]*(ratePass/rateFail)
					else:
						bbeffErrs[process_name] = 1.0
					
			else:
				veffErrs[process_name] = 1.
				bbeffErrs[process_name] = 1.
			mutriggerErrs[process_name] = 1
			muidErrs[process_name] = 1
			muisoErrs[process_name] = 1
			#jesErrs[process_name] = 1
			#jerErrs[process_name] = 1
			if proc=='wqq':
				wznormEWErrs[process_name] = 1.05
			else:
				wznormEWErrs[process_name] = 1.
			if proc=='zqq' or proc=='wqq':
				znormQErrs[process_name] = 1.1
				znormEWErrs[process_name] = 1.15
			else:
				znormQErrs[process_name] = 1.
				znormEWErrs[process_name] = 1.
				
			if rate>0:
				mcStatErrs[process_name] = 1.0+(error[0]/rate)
			else:
				mcStatErrs[process_name] = 1.0
				
			if rate>0:
				rateJESUp = histograms['{}_JESUp'.format(process_name)].Integral()
				rateJESDown = histograms['{}_JESDown'.format(process_name)].Integral()
				rateJERUp = histograms['{}_JERUp'.format(process_name)].Integral()
				rateJERDown = histograms['{}_JERDown'.format(process_name)].Integral()
				ratePuUp = histograms["{}_PuUp".format(process_name)].Integral()
				ratePuDown = histograms["{}_PuDown".format(process_name)].Integral()
				jesErrs[process_name] =  1.0+(abs(rateJESUp-rate)+abs(rateJESDown-rate))/(2.*rate)   
				jerErrs[process_name] =  1.0+(abs(rateJERUp-rate)+abs(rateJERDown-rate))/(2.*rate)
				puErrs[process_name] =  1.0+(abs(ratePuUp-rate)+abs(ratePuDown-rate))/(2.*rate)
			else:
				jesErrs[process_name] =  1.0
				jerErrs[process_name] =  1.0
				puErrs[process_name] =  1.0

	divider = '------------------------------------------------------------\n'
	datacard_file = open(datacard_path,'w')

	datacard_file.write('imax 2 number of channels\n')
	datacard_file.write('jmax * number of processes minus 1\n')
	datacard_file.write('kmax * number of nuisance parameters\n')
	datacard_file.write(divider)
	datacard_file.write('bin fail_muonCR pass_muonCR\n')
	datacard_file.write('observation %.3f %.3f\n'%(obsRate['fail'],obsRate['pass']))
	datacard_file.write(divider)
	datacard_file.write('shapes * pass_muonCR %s w_muonCR:$PROCESS_pass w_muonCR:$PROCESS_pass_$SYSTEMATIC\n'%os.path.basename(workspace_path))
	datacard_file.write('shapes * fail_muonCR %s w_muonCR:$PROCESS_fail w_muonCR:$PROCESS_fail_$SYSTEMATIC\n'%os.path.basename(workspace_path))
	datacard_file.write(divider)
	binString = 'bin'
	processString = 'process'
	processNumberString = 'process'
	rateString = 'rate'
	lumiString = 'lumi\tlnN'
	hqq125ptString = 'hqq125pt\tlnN'
	veffString = 'veff\tlnN'
	bbeffString = 'bbeff\tlnN'
	znormEWString = 'znormEWmuonCR\tlnN'
	znormQString = 'znormQ\tlnN'    
	wznormEWString = 'wznormEWmuonCR\tlnN'
	muidString = 'muid\tshape'   
	muisoString = 'muiso\tshape'   
	mutriggerString = 'mutrigger\tshape'  
	#jesString = 'JES\tshape'    
	#jerString = 'JER\tshape'
	jesString = 'JES\tlnN'
	jerString = 'JER\tlnN'
	puString = 'Pu\tlnN'
	mcStatErrString = {}
	for proc in sigs+bkgs:
		for box in boxes:
			mcStatErrString[process_name] = '%s%smuonCRmcstat\tlnN'%(proc,box)

	for box in boxes:
		i = -1
		for proc in sigs+bkgs:
			i+=1
			if rates[process_name] <= 0.0: continue
			binString +='\t%s_muonCR\n'%box
			processString += '\t%s\n'%(proc)
			processNumberString += '\t%i\n'%(i-nSig+1)
			rateString += '\t%.3f\n'%rates[process_name]
			lumiString += '\t%.3f\n'%lumiErrs[process_name]
			hqq125ptString += '\t%.3f\n'%hqq125ptErrs[process_name]
			veffString += '\t%.3f\n'%veffErrs[process_name]
			bbeffString += '\t%.3f\n'%bbeffErrs[process_name]
			znormEWString += '\t%.3f\n'%znormEWErrs[process_name]
			znormQString += '\t%.3f\n'%znormQErrs[process_name]
			wznormEWString += '\t%.3f\n'%wznormEWErrs[process_name]
			mutriggerString += '\t%.3f\n'%mutriggerErrs[process_name]
			muidString += '\t%.3f\n'%muidErrs[process_name]
			muisoString += '\t%.3f\n'%muisoErrs[process_name]
			jesString += '\t%.3f\n'%jesErrs[process_name]
			jerString += '\t%.3f\n'%jerErrs[process_name]
			puString += '\t%.3f\n'%puErrs[process_name]
			for proc1 in sigs+bkgs:
				for box1 in boxes:
					if proc1==proc and box1==box:
						mcStatErrString['%s_%s'%(proc1,box1)] += '\t%.3f\n'% mcStatErrs[process_name]
					else:                        
						mcStatErrString['%s_%s'%(proc1,box1)] += '\t-\n'
			
	datacard_file.write(binString)
	datacard_file.write(processString)
	datacard_file.write(processNumberString)
	datacard_file.write(rateString)
	datacard_file.write(divider)

	# now nuisances
	datacard_file.write(lumiString)
	datacard_file.write(hqq125ptString)
	datacard_file.write(veffString)
	datacard_file.write(bbeffString)
	datacard_file.write(znormEWString)
	datacard_file.write(znormQString)
	datacard_file.write(wznormEWString)
	datacard_file.write(mutriggerString)
	datacard_file.write(muidString)
	datacard_file.write(muisoString)
	datacard_file.write(jesString)
	datacard_file.write(jerString)
	datacard_file.write(puString)

	for proc in (sigs+bkgs):
		for box in boxes:
			process_name = "{}_{}".format(proc, box)
			if rates[process_name] <= 0.0: 
				continue
			datacard_file.write(mcStatErrString[process_name])

	# now top rate params
	tqqeff = histograms['tqq_pass'].Integral()/(histograms['tqq_pass'].Integral()+histograms['tqq_fail'].Integral())

	
	datacard_file.write('tqqpassmuonCRnorm rateParam pass_muonCR tqq (@0*@1) tqqnormSF,tqqeffSF\n')
	datacard_file.write('tqqfailmuonCRnorm rateParam fail_muonCR tqq (@0*(1.0-@1*%.4f)/(1.0-%.4f)) tqqnormSF,tqqeffSF\n'%(tqqeff,tqqeff))
	datacard_file.write('tqqnormSF extArg 1.0 [0.0,10.0]\n')
	datacard_file.write('tqqeffSF extArg 1.0 [0.0,10.0]\n')
	datacard_file.close()

	
def main(options, args):
	boxes = ['pass', 'fail']
	bkgs = ['zqq','wqq','qcd','tqq','vvqq','stqq','wlnu','zll', "hbb"]
	for signal_name in config.signal_names:
		sigs = [signal_name]
		systs = ['JER','JES','mutrigger','muid','muiso','Pu']
		
		input_file = TFile.Open(options.idir+'/histograms_muCR_{}.root'.format(options.jet_type),'read')
		if not input_file.IsOpen():
			print "ERROR : Couldn't open file at " + options.idir+'/histograms_muCR_{}.root'.format(options.jet_type)
			sys.exit(1)
		input_file.ls()
		histograms = {}
		data_histograms = {}
		
		for proc in (bkgs+sigs+['data_obs']):
			for box in boxes:
				process_name = "{}_{}".format(proc, box)
				print 'getting histogram for process: {}'.format(process_name)

				if not input_file.Get(process_name):
					print "ERROR: Can't load histogram {} from file {}".format(process_name, input_file.GetPath())
					sys.exit(1)
				histograms[process_name] = input_file.Get(process_name).Clone()
				histograms[process_name].Scale(GetSF(proc,box,input_file))
				for syst in systs:
					if proc!='data_obs':
						for direction in ["Up", "Down"]:
							process_name_syst = "{}_{}{}".format(process_name, syst, direction)
							print 'getting histogram for process: {}'.format(process_name_syst)
							if not input_file.Get(process_name_syst):
								print "ERROR : Can't load histogram {} from file {}".format(process_name_syst, input_file.GetPath())
								sys.exit(1)
							histograms[process_name_syst] = input_file.Get(process_name_syst).Clone()
							histograms[process_name_syst].Scale(GetSF(proc,box,input_file))
						
		outFile = 'datacard_muonCR.root'
		
		workspace_path = "{}/{}/{}".format(options.odir, sig, "workspace_muonCR.root")
		datacard_path = "{}/{}/{}".format(options.odir, sig, "datacard_muonCR.txt")
		outputFile = rt.TFile.Open(workspace_path,'recreate')
		outputFile.cd()
		w = rt.RooWorkspace('w_muonCR')
		#w.factory('y[40,40,201]')
		#w.var('y').setBins(1)
		w.factory('x[%i,%i,%i]'%(params[options.jet_type]["MSD"][0],params[options.jet_type]["MSD"][0],params[options.jet_type]["MSD"][1]))
		w.var('x').setBins(params[options.jet_type]["MASS_BIN"])
		for key, histo in histograms.iteritems():
			#histo.Rebin(23)
			#ds = rt.RooDataHist(key,key,rt.RooArgList(w.var('y')),histo)
			ds = rt.RooDataHist(key,key,rt.RooArgList(w.var('x')),histo)
			getattr(w,'import')(ds, rt.RooCmdArg())
		w.Write()
		outputFile.Close()

		writeDataCard(boxes,workspace_path,[sig],bkgs,histograms,options,datacard_path)
		print '\ndatacard:\n'
		os.system('cat {}'.format(datacard_path))



if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
	parser.add_option('--lumi', dest='lumi', type=float, default = 20,help='lumi in 1/fb ', metavar='lumi')
	parser.add_option('-i','--idir', dest='idir', default = './',help='directory with data', metavar='idir')
	parser.add_option('-o','--odir', dest='odir', default = './',help='directory to write cards', metavar='odir')
	parser.add_option('-j', '--jet_type', dest='jet_type', default="AK8", help="AK8 or CA15")
	(options, args) = parser.parse_args()

	main(options, args)
