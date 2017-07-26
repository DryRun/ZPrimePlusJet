#!/usr/bin/env python

import ROOT as r,sys,math,array,os
from ROOT import TFile, TTree, TChain, gPad, gDirectory, gSystem
from multiprocessing import Process
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array

# including other directories
sys.path.insert(0, '../.')
from tools import *

#from buildRhalphabetXbb import cuts[options.jet_type]["MASS_BINS"],cuts[options.jet_type]["MSD"][0],cuts[options.jet_type]["MSD"][1],cuts[options.jet_type]["RHO"][0],cuts[options.jet_type]["RHO"][1]
#from rhalphabet_builder import sfs[options.jet_type]["BB"],sfs[options.jet_type]["BB_ERR"],sfs[options.jet_type]["V"],sfs[options.jet_type]["V_ERR"],GetSF
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libDAZSLEPhiBBPlusJet.so"))
from DAZSLE.PhiBBPlusJet.cuts import cuts
from DAZSLE.PhiBBPlusJet.sfs import sfs

import DAZSLE.PhiBBPlusJet.analysis_configuration as config
from DAZSLE.PhiBBPlusJet.cuts import *

do_syst_mcstat=False
##-------------------------------------------------------------------------------------
def main(options,args):
	for signal_name in config.signal_names:
		tfile = r.TFile.Open(options.ifile)
		boxes = ['pass', 'fail']
		sigs = [signal_name]
		bkgs = ['zqq','wqq','qcd','tqq','hbb']
		systs = [] #['JER','JES']

		nBkgd = len(bkgs)
		nSig = len(sigs)
		numberOfMassBins = cuts[options.jet_type]["MASS_BINS"]   
		numberOfPtBins = len(cuts[options.jet_type]["PT_BINS"])-1

		histoDict = {}

		for proc in (sigs+bkgs):
			for box in boxes:
				print 'getting histogram for process: %s_%s'%(proc,box)
				histoDict['%s_%s'%(proc,box)] = tfile.Get('%s_%s'%(proc,box))
				if not histoDict['%s_%s'%(proc,box)]:
					print "[makeCardsXbb] ERROR : Couldn't find histogram {} in file {}".format('%s_%s'%(proc,box), options.ifile)
					sys.exit(1)
				for syst in systs:
					print 'getting histogram for process: %s_%s_%sUp'%(proc,box,syst)
					histoDict['%s_%s_%sUp'%(proc,box,syst)] = tfile.Get('%s_%s_%sUp'%(proc,box,syst))
					if not histoDict['%s_%s_%sUp'%(proc,box,syst)]:
						print "[makeCardsXbb] ERROR : Couldn't find histogram {} in file {}".format('%s_%s_%sUp'%(proc,box,syst), options.ifile)
						sys.exit(1)
					print 'getting histogram for process: %s_%s_%sDown'%(proc,box,syst)
					histoDict['%s_%s_%sDown'%(proc,box,syst)] = tfile.Get('%s_%s_%sDown'%(proc,box,syst))
					if not histoDict['%s_%s_%sDown'%(proc,box,syst)]:
						print "[makeCardsXbb] ERROR : Couldn't find histogram {} in file {}".format('%s_%s_%sDown'%(proc,box,syst), options.ifile)
						sys.exit(1)

		dctpl = open("datacard_Xbb.tpl")
		#dctpl = open("datacardZbb.tpl")

		linel = [];
		for line in dctpl: 
			#print line.strip().split()
			line = line.replace("SIGNALNAME", signal_name).replace("SIGNALMASS", "")
			linel.append(line.strip())

		for i in range(1,numberOfPtBins+1):

			jesErrs = {}
			jerErrs = {}
			bbErrs = {}
			vErrs = {}
			mcstatErrs = {}
			for box in boxes:
				for proc in (sigs+bkgs):
					#print "Taking integral of {}".format('%s_%s'%(proc,box))
					rate = histoDict['%s_%s'%(proc,box)].Integral(1, numberOfMassBins, i, i)
					if rate>0 and len(systs) > 0:
						rateJESUp = histoDict['%s_%s_JESUp'%(proc,box)].Integral(1, numberOfMassBins, i, i)
						rateJESDown = histoDict['%s_%s_JESDown'%(proc,box)].Integral(1, numberOfMassBins, i, i)
						rateJERUp = histoDict['%s_%s_JERUp'%(proc,box)].Integral(1, numberOfMassBins, i, i)
						rateJERDown = histoDict['%s_%s_JERDown'%(proc,box)].Integral(1, numberOfMassBins, i, i)
						jesErrs['%s_%s'%(proc,box)] =  1.0+(abs(rateJESUp-rate)+abs(rateJESDown-rate))/(2.*rate)   
						jerErrs['%s_%s'%(proc,box)] =  1.0+(abs(rateJERUp-rate)+abs(rateJERDown-rate))/(2.*rate)
					else:
						jesErrs['%s_%s'%(proc,box)] =  1.0
						jerErrs['%s_%s'%(proc,box)] =  1.0
						
					vErrs['%s_%s'%(proc,box)] = 1.0+sfs[options.jet_type]["V_ERR"]/sfs[options.jet_type]["V"]
					if box=='pass':
						bbErrs['%s_%s'%(proc,box)] = 1.0+sfs[options.jet_type]["BB_ERR"]/sfs[options.jet_type]["BB"]
					else:
						ratePass = histoDict['%s_%s'%(proc,'pass')].Integral()
						rateFail = histoDict['%s_%s'%(proc,'fail')].Integral()
						if rateFail>0:
							bbErrs['%s_%s'%(proc,box)] = 1.0-sfs[options.jet_type]["BB_ERR"]*(ratePass/rateFail)
						else:
							bbErrs['%s_%s'%(proc,box)] = 1.0
							
						
					for j in range(1,numberOfMassBins):
						mcstatErrs['%s_%s'%(proc,box),i,j] = 1.0
							

			jesString = 'JES lnN'
			jerString = 'JER lnN'
			bbString = 'bbeff lnN'
			vString = 'veff lnN'
			mcStatStrings = {}
			mcStatGroupString = 'mcstat group ='
			qcdGroupString = 'qcd group = qcdeff'
			for box in boxes:
				for proc in sigs+bkgs:
					for j in range(1,numberOfMassBins):
						mcStatStrings['%s_%s'%(proc,box),i,j] = '%s%scat%imcstat%i shape'%(proc,box,i,j)
						
			for box in boxes:
				for proc in sigs+bkgs:
					if proc=='qcd':
						jesString += ' -'
						jerString += ' -'
					else:
						jesString += ' %.3f'%jesErrs['%s_%s'%(proc,box)]
						jerString += ' %.3f'%jerErrs['%s_%s'%(proc,box)]
					if proc in ['qcd','tqq','wqq', 'ttbar']:
						bbString += ' -'
					else:
						bbString += ' %.3f'%bbErrs['%s_%s'%(proc,box)]
					if proc in ['qcd','tqq', 'ttbar']:
						vString += ' -'
					else:
						vString += ' %.3f'%vErrs['%s_%s'%(proc,box)]
					for j in range(1,numberOfMassBins):
						for box1 in boxes:                    
							for proc1 in sigs+bkgs:                            
								if proc1==proc and box1==box:
									mcStatStrings['%s_%s'%(proc1,box1),i,j] += '\t%d'% mcstatErrs['%s_%s'%(proc,box),i,j]
								else:                        
									mcStatStrings['%s_%s'%(proc1,box1),i,j] += '\t-'

			tag = "cat"+str(i)
			os.system("mkdir -pv " + options.odir + "/{}".format(signal_name))
			dctmp = open(options.odir+"{}/card_rhalphabet_{}.txt".format(signal_name, tag), 'w')
			for l in linel:
				if 'JES' in l:
					if not "JES" in systs:
						continue
					newline = jesString
				elif 'JER' in l:
					if not "JER" in systs:
						continue
					newline = jerString
				elif 'bbeff' in l:
					newline = bbString
				elif 'veff' in l:
					newline = vString
				elif 'TQQEFF' in l:
					if (histoDict['tqq_pass'].Integral() + histoDict['tqq_fail'].Integral()) > 0:
						tqqeff = histoDict['tqq_pass'].Integral() / (histoDict['tqq_pass'].Integral() + histoDict['tqq_fail'].Integral())
					else:
						print "[makeCardsXbb] WARNING : tqq histograms have zero integral! Is this right?"
						tqqeff = 0.
					newline = l.replace('TQQEFF','%.4f'%tqqeff)
				else:
					newline = l
				if "CATX" in l:
					newline = newline.replace('CATX',tag)
				dctmp.write(newline + "\n")
			for box in boxes:
				for proc in sigs+bkgs:
					for j in range(1,numberOfMassBins):                    
						# if stat. unc. is greater than 50% 
						if histoDict['%s_%s'%(proc,box)].GetBinContent(j,i) > 0 and histoDict['%s_%s'%(proc,box)].GetBinError(j,i) > 0.5*histoDict['%s_%s'%(proc,box)].GetBinContent(j,i) and proc!='qcd':
						#if histoDict['%s_%s'%(proc,box)].GetBinContent(j,i) > 0 and proc!='qcd':
							massVal = histoDict['%s_%s'%(proc,box)].GetXaxis().GetBinCenter(j)
							ptVal = histoDict['%s_%s'%(proc,box)].GetYaxis().GetBinLowEdge(i) + 0.3*(histoDict['%s_%s'%(proc,box)].GetYaxis().GetBinWidth(i))
							rhoVal = r.TMath.Log(massVal*massVal/ptVal/ptVal)
							if not( options.blind and massVal > BLIND_LO and massVal < BLIND_HI) and not (rhoVal < cuts[options.jet_type]["RHO"][0] or rhoVal > cuts[options.jet_type]["RHO"][1]):
								if do_syst_mcstat:
									dctmp.write(mcStatStrings['%s_%s'%(proc,box),i,j] + "\n")
								#print 'include %s%scat%imcstat%i'%(proc,box,i,j)
								mcStatGroupString += ' %s%scat%imcstat%i'%(proc,box,i,j)
							else:
								#print 'do not include %s%scat%imcstat%i'%(proc,box,i,j)
								pass
						else:
							#print 'do not include %s%scat%imcstat%i'%(proc,box,i,j)
							pass
							
			for im in range(numberOfMassBins):
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
	parser.add_option("--lumi", dest="lumi", type=float, default = 30,help="luminosity", metavar="lumi")
	parser.add_option('-i','--ifile', dest='ifile', default = 'hist_1DZbb.root',help='file with histogram inputs', metavar='ifile')
	parser.add_option('-o','--odir', dest='odir', default = 'cards/',help='directory to write cards', metavar='odir')
	parser.add_option('--pseudo', action='store_true', dest='pseudo', default =False,help='signal comparison', metavar='isData')
	parser.add_option('--blind', action='store_true', dest='blind', default =False,help='blind signal region', metavar='blind')
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