#!/usr/bin/env python
import ROOT as r,sys,math,os
from ROOT import *
from multiprocessing import Process
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array
#r.gSystem.Load("~/Dropbox/RazorAnalyzer/python/lib/libRazorRun2.so")
r.gSystem.Load(os.getenv('CMSSW_BASE')+'/lib/'+os.getenv('SCRAM_ARCH')+'/libHiggsAnalysisCombinedLimit.so')


# including other directories
sys.path.insert(0, '../.')
sys.path.insert(0, '.')
from rhalphabet_builder import RhalphabetBuilder, LoadHistograms

MASS_BINS = 35
MASS_LO = 40
MASS_HI = 285
BLIND_LO = 110
BLIND_HI = 131

def main(options,args):
    # Load the input histograms
    #   - 2D histograms of pass and fail mass,pT distributions
    #   - for each MC sample and the data
    f = r.TFile.Open(options.ifile)
    (pass_hists,fail_hists) = LoadHistograms(f, options.pseudo, options.blind, options.useQCD, mass_range=[MASS_LO, MASS_HI], blind_range=[BLIND_LO, BLIND_HI]);
    f.Close()

    # Build the workspacees
    RhalphabetBuilder(pass_hists, fail_hists, options.odir, mass_nbins=MASS_BINS, mass_lo=MASS_LO, mass_hi=MASS_HI, blind_lo=BLIND_LO, blind_hi=BLIND_HI, mass_fit=options.mass_fit, freeze_poly=options.freeze_poly)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option('-i','--ifile', dest='ifile', default = 'hist_1DZbb.root',help='file with histogram inputs', metavar='ifile')
    parser.add_option('-o','--odir', dest='odir', default = './',help='directory to write plots', metavar='odir')
    parser.add_option('--pseudo', action='store_true', dest='pseudo', default =False,help='use MC', metavar='pseudo')
    parser.add_option('--blind', action='store_true', dest='blind', default =False,help='blind signal region', metavar='blind')
    parser.add_option('--use-qcd', action='store_true', dest='useQCD', default =False,help='use real QCD MC', metavar='useQCD')
    parser.add_option('--mass_fit', action='store_true', dest='mass_fit', default =False,help='mass fit or rho', metavar='mass_fit')
    parser.add_option('--freeze_poly', action='store_true', dest='freeze_poly', default =False,help='freeze pol values', metavar='freeze')
    
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
