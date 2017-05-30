#!/usr/bin/env python
import ROOT
import sys
import math
import os
from ROOT import *
from multiprocessing import Process
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array

#r.gSystem.Load("~/Dropbox/RazorAnalyzer/python/lib/libRazorRun2.so")
gSystem.Load(os.getenv('CMSSW_BASE') + '/lib/' + os.getenv('SCRAM_ARCH') + '/libHiggsAnalysisCombinedLimit.so')

# including other directories
# sys.path.insert(0, '../.')
sys.path.insert(0, os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/fitting/"))
sys.path.insert(0, '.')
from tools import *
from hist import *
from rhalphabet_builder import RhalphabetBuilder, LoadHistograms, GetSF

MASS_BINS = 70
MASS_LO = 40
MASS_HI = 600
BLIND_LO = 0
BLIND_HI = 0
RHO_LO = -5.5
RHO_HI = -2


def main(options, args):
    ifile = options.ifile
    odir = options.odir

    # Load the input histograms
    # 	- 2D histograms of pass and fail mass,pT distributions
    # 	- for each MC sample and the data
    f = TFile.Open(ifile)
    #(hpass, hfail) = loadHistograms(f, options.pseudo, options.blind, options.useQCD, options.scale, options.r)
    (pass_hists,fail_hists) = LoadHistograms(f, options.pseudo, options.blind, options.useQCD, scale=options.scale, r_signal=options.r, mass_range=[MASS_LO, MASS_HI], blind_range=[0.,0.], rho_range=[RHO_LO,RHO_HI])
    #f.Close()

    # Build the workspacees
    #dazsleRhalphabetBuilder(hpass, hfail, f, odir, options.NR, options.NP)

    rhalphabuilder = RhalphabetBuilder(pass_hists, fail_hists, f, options.odir, nr=options.NR, np=options.NP, mass_nbins=MASS_BINS, mass_lo=MASS_LO, mass_hi=MASS_HI, blind_lo=BLIND_LO, blind_hi=BLIND_HI, rho_lo=RHO_LO, rho_hi=RHO_HI, blind=options.blind, mass_fit=options.massfit, freeze_poly=options.freeze, quiet=True)
    rhalphabuilder.run()

##-------------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option('-i', '--ifile', dest='ifile', default='hist_1DZbb.root', help='file with histogram inputs',
                      metavar='ifile')
    parser.add_option('-o', '--odir', dest='odir', default='./', help='directory to write plots', metavar='odir')
    parser.add_option('--pseudo', action='store_true', dest='pseudo', default=False, help='use MC', metavar='pseudo')
    parser.add_option('--blind', action='store_true', dest='blind', default=False, help='blind signal region',
                      metavar='blind')
    parser.add_option('--useQCD', action='store_true', dest='useQCD', default=False, help='use real QCD MC',
                      metavar='useQCD')
    parser.add_option('--massfit', action='store_true', dest='massfit', default=False, help='mass fit or rho',
                      metavar='massfit')
    parser.add_option('--freeze', action='store_true', dest='freeze', default=False, help='freeze pol values',
                      metavar='freeze')
    parser.add_option('--scale', dest='scale', default=1, type='float',
                      help='scale factor to scale MC (assuming only using a fraction of the data)')
    parser.add_option('--nr', dest='NR', default=2, type='int', help='order of rho (or mass) polynomial')
    parser.add_option('--np', dest='NP', default=2, type='int', help='order of pt polynomial')
    parser.add_option('-r', dest='r', default=0, type='float', help='signal strength for MC pseudodataset')
    (options, args) = parser.parse_args()

    import tdrstyle

    tdrstyle.setTDRStyle()
    gStyle.SetPadTopMargin(0.10)
    gStyle.SetPadLeftMargin(0.16)
    gStyle.SetPadRightMargin(0.10)
    gStyle.SetPalette(1)
    gStyle.SetPaintTextFormat("1.1f")
    gStyle.SetOptFit(0000)
    gROOT.SetBatch()

    main(options, args)
##-------------------------------------------------------------------------------------
