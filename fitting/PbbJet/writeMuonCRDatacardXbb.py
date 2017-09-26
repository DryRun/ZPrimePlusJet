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

import DAZSLE.ZPrimePlusJet.xbb_config as config
from buildRhalphabetXbb import GetSF

# fake_signal: instead of loading the signal template, insert a fake template, gaussian with tiny normalization. Using this until we can come up with a good template.
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
	print sigs+bkgs
	for proc in sigs+bkgs:
		for box in boxes:
			process_name = "{}_{}".format(proc, box)
			error = array.array('d',[0.0])
			rate = histograms[process_name].IntegralAndError(1,histograms[process_name].GetNbinsX(),error)
			rates[process_name]  = rate
			lumiErrs[process_name] = 1.025
			if proc == "hqq125":
				hbb125ptErrs[process_name] = 1.3
			else:
				hbb125ptErrs[process_name] = 1.0
			if proc in ["wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"] or ("Sbb" in proc) or ("ZPrime" in proc):
				veffErrs[process_name] = 1.0+config.analysis_parameters[options.jet_type]["V_SF_ERR"]/config.analysis_parameters[options.jet_type]["V_SF"]
				if box=='pass':
					bbeffErrs[process_name] = 1.0+config.analysis_parameters[options.jet_type]["BB_SF_ERR"]/config.analysis_parameters[options.jet_type]["BB_SF"]
				else:
					ratePass = histograms['%s_%s'%(proc,'pass')].Integral()
					rateFail = histograms['%s_%s'%(proc,'fail')].Integral()
					if rateFail>0:
						bbeffErrs[process_name] = 1.0-config.analysis_parameters[options.jet_type]["BB_SF_ERR"]*(ratePass/rateFail)
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
				ratePuUp = histograms["{}_PUUp".format(process_name)].Integral()
				ratePuDown = histograms["{}_PUDown".format(process_name)].Integral()
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
	hbb125ptString = 'hqq125pt\tlnN'
	veffString = 'veff\tlnN'
	bbeffString = 'bbeff\tlnN'
	znormEWString = 'znormEWmuonCR\tlnN'
	znormQString = 'znormQ\tlnN'    
	wznormEWString = 'wznormEWmuonCR\tlnN'
	muidString = 'MuID\tshape'   
	muisoString = 'MuIso\tshape'   
	mutriggerString = 'MuTrigger\tshape'  
	#jesString = 'JES\tshape'    
	#jerString = 'JER\tshape'
	jesString = 'JES\tlnN'
	jerString = 'JER\tlnN'
	puString = 'PU\tlnN'
	mcStatErrString = {}
	for proc in sigs+bkgs:
		for box in boxes:
			process_name = "{}_{}".format(proc, box)
			mcStatErrString[process_name] = '%s%smuonCRmcstat\tlnN'%(proc,box)

	for box in boxes:
		i = -1
		for proc in sigs+bkgs:
			process_name = "{}_{}".format(proc, box)
			i+=1
			if rates[process_name] <= 0.0: 
				print "Skipping process {} in datacard because rate = {}".format(process_name, rates[process_name])
				continue
			binString +='\t%s_muonCR'%box
			processString += '\t%s'%(proc)
			processNumberString += '\t%i'%(i-nSig+1)
			rateString += '\t%.3f'%rates[process_name]
			lumiString += '\t%.3f'%lumiErrs[process_name]
			hbb125ptString += '\t%.3f'%hbb125ptErrs[process_name]
			veffString += '\t%.3f'%veffErrs[process_name]
			bbeffString += '\t%.3f'%bbeffErrs[process_name]
			znormEWString += '\t%.3f'%znormEWErrs[process_name]
			znormQString += '\t%.3f'%znormQErrs[process_name]
			wznormEWString += '\t%.3f'%wznormEWErrs[process_name]
			mutriggerString += '\t%.3f'%mutriggerErrs[process_name]
			muidString += '\t%.3f'%muidErrs[process_name]
			muisoString += '\t%.3f'%muisoErrs[process_name]
			jesString += '\t%.3f'%jesErrs[process_name]
			jerString += '\t%.3f'%jerErrs[process_name]
			puString += '\t%.3f'%puErrs[process_name]
			for proc1 in sigs+bkgs:
				for box1 in boxes:
					if proc1==proc and box1==box:
						print mcStatErrs
						mcStatErrString['%s_%s'%(proc1,box1)] += '\t%.3f'% mcStatErrs[process_name]
					else:                        
						mcStatErrString['%s_%s'%(proc1,box1)] += '\t-'
	binString += "\n"
	processString += "\n"
	processNumberString += "\n"
	rateString += "\n"
	lumiString += "\n"
	hbb125ptString += "\n"
	veffString += "\n"
	bbeffString += "\n"
	znormEWString += "\n"
	znormQString += "\n"
	wznormEWString += "\n"
	mutriggerString += "\n"
	muidString += "\n"
	muisoString += "\n"
	jesString += "\n"
	jerString += "\n"
	puString += "\n"
	for proc in (sigs+bkgs):
		for box in boxes:
			mcStatErrString['%s_%s'%(proc,box)] += '\n'
	datacard_file.write(binString)
	datacard_file.write(processString)
	datacard_file.write(processNumberString)
	datacard_file.write(rateString)
	datacard_file.write(divider)

	# now nuisances
	datacard_file.write(lumiString)
	datacard_file.write(hbb125ptString)
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

	# Copy to subdirs
	print "Wrote muon CR datacards to " + datacard_path

	
def main(options, args):
	boxes = ['pass', 'fail']
	bkgs = ['zqq','wqq','qcd','tqq','vvqq','stqq','wlnu','zll', "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"]
	for signal_name in config.limit_signal_names[options.jet_type]:
		output_directory = config.get_datacard_directory(signal_name, options.jet_type, qcd=True)
		sigs = [signal_name]
		systs = ['JER','JES','MuTrigger','MuID','MuIso','PU']
		print config.interpolated_signal_names
		if signal_name in config.interpolated_signal_names:
			signal_file = TFile.Open(options.idir+'/interpolations_muCR_{}.root'.format(options.jet_type),'read')
		else:
			signal_file = TFile.Open(options.idir+'/histograms_muCR_{}.root'.format(options.jet_type),'read')
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

				if proc in sigs:
					this_file = signal_file
				else:
					this_file = input_file
				if not this_file.Get(process_name):
					print "[writeMuonCRDatacardXbb::main] ERROR : Can't load histogram {} from file {}".format(process_name, this_file.GetPath())
					sys.exit(1)
				histograms[process_name] = this_file.Get(process_name)
				histograms[process_name].Scale(GetSF(proc,box,options.jet_type,this_file))
				for syst in systs:
					if proc!='data_obs':
						for direction in ["Up", "Down"]:
							process_name_syst = "{}_{}{}".format(process_name, syst, direction)
							print 'getting histogram for process: {}'.format(process_name_syst)
							if not this_file.Get(process_name_syst):
								print "[writeMuonCRDatacardXbb::main] ERROR : Can't load histogram {} from file {}".format(process_name_syst, this_file.GetPath())
								sys.exit(1)
							histograms[process_name_syst] = this_file.Get(process_name_syst)
							histograms[process_name_syst].Scale(GetSF(proc,box,options.jet_type,this_file))

		if options.fake_signal:
			for proc in sigs:
				for box in boxes:
					process_name = "{}_{}".format(proc, box)
					sig_mass = config.signal_masses[proc]
					sig_width = sig_mass / 10.
					histograms[process_name].Reset()
					gaussian = TF1("mygaus","TMath::Gaus(x,{},{})".format(sig_mass, sig_width), 0., 1000.);
					histograms[process_name].FillRandom("mygaus", 1000)
					histograms[process_name].Scale(0.01 / histograms[process_name].Integral())
					
					for syst in systs:
						for direction in ["Up", "Down"]:
							process_name_syst = "{}_{}{}".format(process_name, syst, direction)
							histograms[process_name_syst].Reset()
							histograms[process_name_syst].FillRandom("mygaus", 1000)
							histograms[process_name_syst].Scale(0.01 / histograms[process_name_syst].Integral())


		outFile = 'datacard_muonCR.root'
		
		os.system("mkdir -pv {}".format(output_directory))
		workspace_path = "{}/{}".format(output_directory, "workspace_muonCR.root")
		datacard_path = "{}/{}".format(output_directory, "datacard_muonCR.txt")
		outputFile = TFile.Open(workspace_path,'recreate')
		outputFile.cd()
		w = RooWorkspace('w_muonCR')
		#w.factory('y[40,40,201]')
		#w.var('y').setBins(1)
		w.factory('x[%i,%i,%i]'%(config.analysis_parameters[options.jet_type]["MSD"][0], config.analysis_parameters[options.jet_type]["MSD"][0],config.analysis_parameters[options.jet_type]["MSD"][1]))
		w.var('x').setBins(config.analysis_parameters[options.jet_type]["MASS_BINS"])
		for key, histo in histograms.iteritems():
			#histo.Rebin(23)
			#ds = RooDataHist(key,key,RooArgList(w.var('y')),histo)
			ds = RooDataHist(key,key,RooArgList(w.var('x')),histo)
			getattr(w,'import')(ds, RooCmdArg())
		w.Write()
		outputFile.Close()

		writeDataCard(boxes,workspace_path,[signal_name],bkgs,histograms,options,datacard_path)
		print '\ndatacard:\n'
		os.system('cat {}'.format(datacard_path))



if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-i','--idir', dest='idir', default='/uscms/home/dryu/DAZSLE/data/LimitSetting/Xbb_inputs/',help='directory with data', metavar='idir')
	parser.add_option('-j', '--jet_type', dest='jet_type', default="AK8", help="AK8 or CA15")
	parser.add_option('--fake_signal', action='store_true', help="Use fake signal template (gaussian with tiny normalization)")
	(options, args) = parser.parse_args()

	main(options, args)
