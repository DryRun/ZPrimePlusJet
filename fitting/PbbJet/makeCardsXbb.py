#!/usr/bin/env python
import sys
import math
import array
import os
from ROOT import TFile, TTree, TChain, gPad, gDirectory, gSystem
from multiprocessing import Process
from optparse import OptionParser
from operator import add
import time
import array

# including other directories
from DAZSLE.ZPrimePlusJet.tools import *
from DAZSLE.ZPrimePlusJet.xbb_config import analysis_parameters as params

gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libDAZSLEPhiBBPlusJet.so"))
import DAZSLE.PhiBBPlusJet.analysis_configuration as config

do_syst_mcstat=False
##-------------------------------------------------------------------------------------
def main(options,args):
	for signal_name in config.signal_names:
		input_file = TFile.Open(options.ifile)
		boxes = ['pass', 'fail']
		sigs = [signal_name]
		bkgs = ['zqq','wqq','qcd','tqq','hbb']
		systs = ["JER", "JES", "Pu"] #['JER','JES']

		nBkgd = len(bkgs)
		nSig = len(sigs)
		len(params[options.jet_type]["PT_BINS"])-1 = len(cuts[options.jet_type]["PT_BINS"])-1

		histograms = {}

		for proc in (sigs+bkgs):
			for box in boxes:
				process_name = "{}_{}".format(proc,box)
				print 'getting histogram for process: {}'.format(process_name)
				histograms[process_name] = input_file.Get(process_name)
				if not histograms[process_name]:
					print "[makeCardsXbb] ERROR : Couldn't find histogram {} in file {}".format(process_name, options.ifile)
					sys.exit(1)
				for syst in systs:
					for direction in ["Up", "Down"]:
						process_name_syst = "{}_{}{}".format(process_name, syst, direction)
						print 'getting histogram for process: {}'%(process_name_syst)
						histograms[process_name_syst] = input_file.Get('{}'.format(process_name_syst))
						if not histograms[process_name_syst]:
							print "[makeCardsXbb] ERROR : Couldn't find histogram {} in file {}".format(process_name_syst, options.ifile)
							sys.exit(1)

		dctpl = open("datacardPhibb.tpl")
		#dctpl = open("datacardZbb.tpl")

		linel = [];
		for line in dctpl: 
			#print line.strip().split()
			line = line.replace("SIGNALNAMESIGNALMASS", signal_name)
			linel.append(line.strip())

		for i in range(1,len(params[options.jet_type]["PT_BINS"])-1+1):

			errs = {}
			#jesErrs = {}
			#jerErrs = {}
			#bbErrs = {}
			#vErrs = {}
			#mcstatErrs = {}
			#scaleptErrs = {}
			for box in boxes:
				for proc in (sigs+bkgs):
					#print "Taking integral of {}".format('%s_%s'%(proc,box))
					process_name = "{}_{}".format(proc, box)
					errs[process_name] = {}
					rate_central = histograms[process_name].Integral(1, params[options.jet_type]["MASSBINS"], i, i)
					if rate>0:
						for syst in systs:
							rate_up = histograms["{}_{}Up".format(process_name, syst)].Integral(1, params[options.jet_type]["MASSBINS"], i, i)
							rate_down = histograms["{}_{}Down".format(process_name, syst)].Integral(1, params[options.jet_type]["MASSBINS"], i, i)
							errs[process_name][syst] = 1.0 + (abs(rate_up - rate_central) + abs(rate_down - rate_central)) / (2. * rate_central)
					else:
						for syst in systs:
							errs[process_name][syst] = 1.0
					if i == 2:
						errs[process_name]["scalept"] =  0.05
					elif i == 3:
						errs[process_name]["scalept"] =  0.1
					elif i == 4:
						errs[process_name]["scalept"] =  0.2
					elif i == 5:
						errs[process_name]["scalept"] =  0.3
					elif i == 6:
						errs[process_name]["scalept"] =  0.4
					errs[process_name]["veff"] = 1.0+params[options.jet_type]["V_SF_ERR"] / params[options.jet_type]["V_SF"]
					if box=='pass':
						errs[process_name]["bbeff"] = 1.0+params[options.jet_type]["BB_SF_ERR"] / params[options.jet_type]["BB_SF"]
					else:
						ratePass = histograms['%s_%s'%(proc,'pass')].Integral()
						rateFail = histograms['%s_%s'%(proc,'fail')].Integral()
						if rateFail>0:
							errs[process_name]["bbeff"] = 1.0-params[options.jet_type]["BB_SF_ERR"]*(ratePass/rateFail)
						else:
							errs[process_name]["bbeff"] = 1.0

					# MC stat							
					errs[process_name]["mcstat"] = {}
					for j in range(1,params[options.jet_type]["MASSBINS"]):
						if options.noMcStatShape:
							error = array.array('d',[0.0])
							rate = histograms[process_name].IntegralAndError(1,histograms[process_name].GetNbinsX(),i,i,error)                 
							#mcstatErrs['%s_%s'%(proc,box),i,j] = 1.0+histograms[process_name].GetBinError(j,i)/histograms[process_name].Integral()
							errs[process_name]["mcstat"][j] = 1.0+(error[0]/rate)
						else:
							errs[process_name]["mcstat"][j] = 1.0
							

			jesString = 'JES lnN'
			jerString = 'JER lnN'
			bbString = 'bbeff lnN'
			vString = 'veff lnN'
			scaleptString = "scalept shape"
			mcStatStrings = {}
			mcStatGroupString = 'mcstat group ='
			qcdGroupString = 'qcd group = qcdeff'
			for box in boxes:
				for proc in sigs+bkgs:
					process_name = "{}_{}".format(proc, box)
					for j in range(1,params[options.jet_type]["MASSBINS"]):
						if options.noMcStatShape:
							mcStatStrings[process_name,i,j] = '{}cat{}mcstat{} lnN'%(process_name,i,j)
						else:
							mcStatStrings[process_name,i,j] = '{}cat{}mcstat{} shape'%(process_name,i,j)

			for box in boxes:
				for proc in sigs+bkgs:
					process_name = "{}_{}".format(proc, box)
					if proc=='qcd':
						jesString += ' -'
						jerString += ' -'
						jerString += ' -'
					else:
						jesString += " {:.3f}".format(errs[process_name]["JES"])
						jerString += " {:.3f}".format(errs[process_name]["JER"])
						puString += " {:.3f}".format(errs[process_name]["Pu"])

					if proc in ["qcd", "tqq"]:
						if i > 1:
							scaleptString += " -"
					else:
						if i > 1:
							scaleptString += ' {:.3f}'.format(errs[process_name]["scalept"])

					if proc in ['qcd', 'tqq', 'wqq']:
						bbString += ' -'
					else:
						bbString += ' {:.3f}'.format(errs[process_name]["bbeff"])
					if proc in ['qcd','tqq']:
						vString += ' -'
					else:
						vString += ' {:.3f}'.format(errs[process_name]["veff"])
					for j in range(1,params[options.jet_type]["MASSBINS"]):
						for box1 in boxes:
							for proc1 in sigs+bkgs:
								if proc1==proc and box1==box:
									mcStatStrings['%s_%s'%(proc1,box1),i,j] += '\t{:.3f}'.format(errs[process_name]["mcstat"][j])
								else:									mcStatStrings['%s_%s'%(proc1,box1),i,j] += '\t-'

			tag = "cat"+str(i)
			os.system("mkdir -pv " + options.odir + "/{}".format(signal_name))
			dctmp = open(options.odir+"{}/card_rhalphabet_{}.txt".format(signal_name, tag), 'w')
			for l in linel:
				if 'JES' in l:
					newline = jesString
				elif 'JER' in l:
					newline = jerString
				elif "Pu" in l:
					newline = puString
				elif 'bbeff' in l:
					newline = bbString
				elif 'veff' in l:
					newline = vString
				elif "scalept" in l and i>1:
					newline = scaleptString
				elif 'TQQEFF' in l:
					if (histograms['tqq_pass'].Integral() + histograms['tqq_fail'].Integral()) > 0:
						tqqeff = histograms['tqq_pass'].Integral() / (histograms['tqq_pass'].Integral() + histograms['tqq_fail'].Integral())
					else:
						print "[makeCardsXbb] WARNING : tqq histograms have zero integral! Is this right?"
						tqqeff = 0.
					newline = l.replace('TQQEFF','%.4f'%tqqeff)
				elif 'wznormEW' in l:
					if i==4:
						newline = l.replace('1.05','1.15')
					elif i==5:
						newline = l.replace('1.05','1.15')
					elif i==6:
						newline = l.replace('1.05','1.15')
					else:
						newline = l
				elif 'znormEW' in l:
					if i==3:
						newline = l.replace('1.15','1.25')
					elif i==4:
						newline = l.replace('1.15','1.35')
					elif i==5:
						newline = l.replace('1.15','1.35')
					elif i==6:
						newline = l.replace('1.15','1.35')
					else:
						newline = l              
				else:
					newline = l
				if "CATX" in l:
					newline = newline.replace('CATX',tag)
				dctmp.write(newline + "\n")
			for box in boxes:
				for proc in sigs+bkgs:
					process_name = "{}_{}".format(proc, box)
					if options.noMcStatShape and proc != "qcd":
						#print 'include {}cat{}mcstat'%(process_name, i)
						dctmp.write(mcStatStrings[process_name,i,1].replace('mcstat1','mcstat') + "\n")
						mcStatGroupString += ' {}cat{}mcstat'.format(process_name, i)
					else:
						for j in range(1,params[options.jet_type]["MASSBINS"]):                    
							# if stat. unc. is greater than 50% 
							if abs(histograms[process_name].GetBinContent(j,i)) > 0 and histograms[process_name].GetBinError(j,i) > 0.5*histograms[process_name].GetBinContent(j,i) and proc!='qcd':
								massVal = histograms[process_name].GetXaxis().GetBinCenter(j)
								ptVal = histograms[process_name].GetYaxis().GetBinLowEdge(i) + 0.3*(histograms[process_name].GetYaxis().GetBinWidth(i))
								rhoVal = TMath.Log(massVal*massVal/ptVal/ptVal)
								if not(options.blind and massVal > BLIND_LO and massVal < BLIND_HI) and not (rhoVal < params[options.jet_type]["RHO"][0] or rhoVal > params[options.jet_type]["RHO"][1]):
									dctmp.write(mcStatStrings[process_name,i,j] + "\n")
									#print 'include %s%scat%imcstat%i'%(proc,box,i,j)
									mcStatGroupString += ' %s%scat%imcstat%i'%(proc,box,i,j)
								else:
									#print 'do not include %s%scat%imcstat%i'%(proc,box,i,j)
									pass
							else:
								#print 'do not include %s%scat%imcstat%i'%(proc,box,i,j)
								pass
							
			for im in range(params[options.jet_type]["MASSBINS"]):
				dctmp.write("qcd_fail_%s_Bin%i flatParam \n" % (tag,im+1))
				qcdGroupString += ' qcd_fail_%s_Bin%i'%(tag,im+1)
			if do_syst_mcstat:
				dctmp.write(mcStatGroupString + "\n")
			dctmp.write(qcdGroupString + "\n")


###############################################################


	
##-------------------------------------------------------------------------------------
if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
	parser.add_option('-i','--ifile', dest='ifile', default = 'hist_1DZbb.root',help='file with histogram inputs', metavar='ifile')
	parser.add_option('-o','--odir', dest='odir', default = 'cards/',help='directory to write cards', metavar='odir')
	parser.add_option('--pseudo', action='store_true', dest='pseudo', default =False,help='signal comparison', metavar='isData')
	parser.add_option('--blind', action='store_true', dest='blind', default =False,help='blind signal region', metavar='blind')
	parser.add_option('--remove-unmatched', action='store_true', dest='removeUnmatched', default =False,help='remove unmatched', metavar='removeUnmatched')
	parser.add_option('--no-mcstat-shape', action='store_true', dest='noMcStatShape', default =False,help='change mcstat uncertainties to lnN', metavar='noMcStatShape')
	parser.add_option('--jet_type', type=str, help='AK8 or CA15')

	(options, args) = parser.parse_args()

	import tdrstyle
	tdrstyle.setTDRStyle()
	r.gStyle.SetPadTopMargin(0.10)
	r.gStyle.SetPadLeftMargin(0.16)
	r.gStyle.SetPadRightMargin(0.10)
	r.gStyle.SetPalette(1)
	r.gStyle.SetPaintTextFormat("1.1f")
	r.gStyle.SetOptFit(0000)
	r.gROOT.SetBatch()

	main(options,args)
##-------------------------------------------------------------------------------------
