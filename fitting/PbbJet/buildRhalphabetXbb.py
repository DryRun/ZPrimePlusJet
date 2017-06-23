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
from rhalphabet_builder import RhalphabetBuilder

gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libDAZSLEPhiBBPlusJet.so"))
from DAZSLE.PhiBBPlusJet.cuts import cuts
from DAZSLE.PhiBBPlusJet.sfs import sfs
import DAZSLE.PhiBBPlusJet.analysis_configuration as config

# Scale factors for MC
# - Double b-tag SF from muon control region
# - Vector(bb) sf... what is this again?
def get_sf(process, cat, jet_type, f, fLoose=None, removeUnmatched=False, iPt=-1):
    SF = 1.
    print process, cat
    if 'hbb' in process or 'zqq' in process or 'Pbb' in process or 'Sbb' in process:
        if 'pass' in cat:
            SF *= sfs[jet_type]["BB"]
            if 'zqq' in process:
                print sfs[jet_type]["BB"]
        else:
            passInt = f.Get(process + '_pass').Integral()
            failInt = f.Get(process + '_fail').Integral()
            if failInt > 0:
                SF *= (1. + (1. - sfs[jet_type]["BB"]) * passInt / failInt)
                if 'zqq' in process:
                    print (1. + (1. - sfs[jet_type]["BB"]) * passInt / failInt)
    if 'wqq' in process or 'zqq' in process or 'hbb' in process or 'Pbb' in process or 'Sbb' in process:
        SF *= sfs[jet_type]["V"]
        if 'zqq' in process:
            print sfs[jet_type]["V"]
    matchingString = ''
    if removeUnmatched and ('wqq' in process or 'zqq' in process):
        matchingString = '_matched'
    if fLoose is not None and ('wqq' in process or 'zqq' in process) and 'pass' in cat:
        if iPt > -1:
            nbinsX = f.Get(process + '_pass' + matchingString).GetXaxis().GetNbins()
            passInt = f.Get(process + '_pass' + matchingString).Integral(1, nbinsX, int(iPt), int(iPt))
            passIntLoose = fLoose.Get(process + '_pass' + matchingString).Integral(1, nbinsX, int(iPt), int(iPt))
        else:
            passInt = f.Get(process + '_pass' + matchingString).Integral()
            passIntLoose = fLoose.Get(process + '_pass' + matchingString).Integral()
        SF *= passInt / passIntLoose
        if 'zqq' in process:
            print passInt / passIntLoose
    # remove cross section from MH=125 signal templates (template normalized to luminosity*efficiency*acceptance)
    ## if process=='hbb125':
    ##     SF *= 1./48.85*5.824E-01
    ## elif process=='zhbb':
    ##     SF *= 1./(8.839E-01*(1.-3.*0.0335962-0.201030)*5.824E-01+8.839E-01*5.824E-01*0.201030+1.227E-01*5.824E-01*0.201030+1.227E-01*5.824E-01*0.201030)
    ## elif process=='whbb':
    ##     SF *= 1./(5.328E-01*(1.-3.*0.108535)*5.824E-01+8.400E-01*(1.-3.*0.108535)*5.824E-01)
    ## elif process=='vbfhbb':
    ##     SF *= 1./(3.782*5.824E-01)
    ## elif process=='tthbb':
    ##     SF *= 1./(5.071E-01*5.824E-01)

    # if 'zqq' in process:
    #    print SF
    #    sys.exit()
    return SF


# Load histograms here, rather than rhalphabet_builder.LoadHistograms. I think it's unreasonable to expect that different analyses will share this function!
# input_file_path = path of TFile containing the pass/fail msd-pT histograms
# mass_range      = range of mSD to use
# rho_range       = range of rho to use
# pseudo          = Substitute pseudodata constructed from MC for real data
# useQCD = For QCD pass, use MC prediction instead of fail * (pass int / fail int)
# fLoose = 
def load_histograms(input_file, mass_range, rho_range, jet_type=None, pseudo=False, useQCD=False, loose_file_path=None, scale=1., r_signal=0., blind_range=None):
    pass_hists = {}
    fail_hists = {}
    #f.ls()

    # backgrounds
    pass_hists_bkg = {}
    fail_hists_bkg = {}
    background_names = ["wqq", "zqq", "qcd", "tqq", "hbb"]
    for i, bkg in enumerate(background_names):
        if bkg == 'qcd':
            qcd_fail = input_file.Get('qcd_fail')
            qcd_fail.Scale(1. / scale)
            qcd_fail.SetBinContent(13, 4, (
                qcd_fail.GetBinContent(12, 4) + qcd_fail.GetBinContent(14, 4)) / 2.)  # REMOVE HIGH WEIGHT EVENT BIN
            qcd_fail.SetBinError(13, 4, (
                qcd_fail.GetBinError(12, 4) + qcd_fail.GetBinError(14, 4)) / 2.)  # REMOVE HIGH WEIGHT EVENT BIN
            if useQCD:
                qcd_pass = input_file.Get('qcd_pass').Clone()
                qcd_pass.Scale(1. / scale)
            else:
                qcd_pass_real = input_file.Get('qcd_pass').Clone('qcd_pass_real')
                qcd_pass_real.Scale(1. / scale)
                qcd_pass = qcd_fail.Clone('qcd_pass')
                qcd_pass_real_integral = 0
                qcd_fail_integral = 0
                for i in range(1, qcd_pass_real.GetNbinsX() + 1):
                    for j in range(1, qcd_pass_real.GetNbinsY() + 1):
                        if qcd_pass_real.GetXaxis().GetBinCenter(i) > mass_range[
                            0] and qcd_pass_real.GetXaxis().GetBinCenter(
                                i) < mass_range[1]:
                            qcd_pass_real_integral += qcd_pass_real.GetBinContent(i, j)
                            qcd_fail_integral += qcd_fail.GetBinContent(i, j)
                if qcd_fail_integral > 0:
                    qcd_pass.Scale(qcd_pass_real_integral / qcd_fail_integral)  # qcd_pass = qcd_fail * eff(pass)/eff(fail)
                else:
                    print "[LoadHistograms] ERROR : qcd_fail_integral = 0" 
                    print "[LoadHistograms] ERROR : file = " + input_file.GetPath()
                    sys.exit(1)
            pass_hists_bkg["qcd"] = qcd_pass
            fail_hists_bkg["qcd"] = qcd_fail
            print 'qcd pass integral', qcd_pass.Integral()
            print 'qcd fail integral', qcd_fail.Integral()
        elif loose_file_path and (bkg == 'wqq' or bkg == 'zqq'):
            input_file_loose = TFile(loose_file_path, "READ")
            hpass_tmp = input_file_loose.Get(bkg + '_pass').Clone()
            hfail_tmp = input_file.Get(bkg + '_fail').Clone()
            hpass_tmp.Scale(1. / scale)
            hfail_tmp.Scale(1. / scale)
            hpass_tmp.Scale(get_sf(bkg, 'pass', jet_type, input_file, input_file_loose))
            hfail_tmp.Scale(get_sf(bkg, 'fail', jet_type, input_file))
            pass_hists_bkg[bkg] = hpass_tmp
            pass_hists_bkg[bkg].SetDirectory(0)
            fail_hists_bkg[bkg] = hfail_tmp
            fail_hists_bkg[bkg].SetDirectory(0)
            input_file_loose.Close()
        else:
            print "Attempting to get {}".format(bkg + '_pass')
            hpass_tmp = input_file.Get(bkg + '_pass').Clone()
            hfail_tmp = input_file.Get(bkg + '_fail').Clone()
            hpass_tmp.Scale(1. / scale)
            hfail_tmp.Scale(1. / scale)
            hpass_tmp.Scale(get_sf(bkg, 'pass', jet_type, input_file))
            hfail_tmp.Scale(get_sf(bkg, 'fail', jet_type, input_file))
            pass_hists_bkg[bkg] = hpass_tmp
            fail_hists_bkg[bkg] = hfail_tmp

    # signals
    pass_hists_sig = {}
    fail_hists_sig = {}

    for signal_name in config.signal_names:
        print "[debug] Getting " + signal_name + "_pass"
        if not input_file.Get(signal_name + "_pass"):
            print "[LoadHistograms] ERROR : Couldn't load histogram " + signal_name + "_pass" + " from file " + input_file.GetPath()
            sys.exit(1)
        passhist = input_file.Get(signal_name + "_pass").Clone()
        failhist = input_file.Get(signal_name + "_fail").Clone()
        for hist in [passhist, failhist]:
            for i in range(0, hist.GetNbinsX() + 2):
                for j in range(0, hist.GetNbinsY() + 2):
                    if hist.GetBinContent(i, j) <= 0:
                        hist.SetBinContent(i, j, 0)
        failhist.Scale(1. / scale)
        passhist.Scale(1. / scale)
        failhist.Scale(get_sf(signal_name, 'fail', jet_type, input_file))
        passhist.Scale(get_sf(signal_name, 'pass', jet_type, input_file))

        pass_hists_sig[signal_name] = passhist
        fail_hists_sig[signal_name] = failhist
        #signal_names.append(signal_name)

    if pseudo:
        for i, bkg in enumerate(background_names):
            print "Pseudo-making data: adding " + bkg
            print "\t",
            print pass_hists_bkg[bkg]
            if i == 0:
                pass_hists["data_obs"] = pass_hists_bkg[bkg].Clone('data_obs_pass')
                fail_hists["data_obs"] = fail_hists_bkg[bkg].Clone('data_obs_fail')
            else:
                pass_hists["data_obs"].Add(pass_hists_bkg[bkg])
                fail_hists["data_obs"].Add(fail_hists_bkg[bkg])

        for i, signal_name in enumerate(config.signal_names):
            pass_hists["data_obs"].Add(pass_hists_sig[signal_name], r_signal)
            fail_hists["data_obs"].Add(fail_hists_sig[signal_name], r_signal)
    else:
        pass_hists["data_obs"] = input_file.Get('data_obs_pass')
        fail_hists["data_obs"] = input_file.Get('data_obs_fail')

    pass_hists.update(pass_hists_bkg)
    pass_hists.update(pass_hists_sig)
    fail_hists.update(fail_hists_bkg)
    fail_hists.update(fail_hists_sig)

    print pass_hists
    print fail_hists
    for histogram in (pass_hists.values() + fail_hists.values()):
        for i in range(1,histogram.GetNbinsX()+1):
            for j in range(1,histogram.GetNbinsY()+1):
                massVal = histogram.GetXaxis().GetBinCenter(i)
                ptVal = histogram.GetYaxis().GetBinLowEdge(j) + histogram.GetYaxis().GetBinWidth(j) * 0.3
                rhoVal = r.TMath.Log(massVal * massVal / ptVal / ptVal)
                if blind_range: 
                    if histogram.GetXaxis().GetBinCenter(i) > blind_range[0] and histogram.GetXaxis().GetBinCenter(i) < blind_range[1]:
                        print "blinding signal region for %s, mass bin [%i,%i] " % (
                        histogram.GetName(), histogram.GetXaxis().GetBinLowEdge(i), histogram.GetXaxis().GetBinUpEdge(i))
                        histogram.SetBinContent(i, j, 0.)
                if rhoVal < rho_range[0] or rhoVal > rho_range[1]:
                    print "removing rho = %.2f for %s, pt bin [%i, %i], mass bin [%i,%i]" % (
                        rhoVal, histogram.GetName(), histogram.GetYaxis().GetBinLowEdge(j),
                        histogram.GetYaxis().GetBinUpEdge(j),
                        histogram.GetXaxis().GetBinLowEdge(i), histogram.GetXaxis().GetBinUpEdge(i))
                    histogram.SetBinContent(i, j, 0.)
        histogram.SetDirectory(0)

        # print "lengths = ", len(pass_hists), len(fail_hists)
    # print pass_hists;
    # print fail_hists;
    return (pass_hists, fail_hists)


def main(options, args):
    input_file = TFile(options.ifile, "READ")
    odir = options.odir

    # Load the input histograms
    # 	- 2D histograms of pass and fail mass,pT distributions
    # 	- for each MC sample and the data
    (pass_hists,fail_hists) = load_histograms(input_file, cuts[options.jet_type]["MSD"], cuts[options.jet_type]["RHO"], useQCD=options.useQCD, jet_type=options.jet_type, scale=options.scale, r_signal=options.r, pseudo=options.pseudo)

    # Build the workspacees
    #dazsleRhalphabetBuilder(hpass, hfail, f, odir, options.NR, options.NP)

    rhalphabuilder = RhalphabetBuilder(pass_hists, fail_hists, input_file, options.odir, nr=options.NR, np=options.NP, mass_nbins=cuts[options.jet_type]["MASS_BINS"], mass_lo=cuts[options.jet_type]["MSD"][0], mass_hi=cuts[options.jet_type]["MSD"][1], rho_lo=cuts[options.jet_type]["RHO"][0], rho_hi=cuts[options.jet_type]["RHO"][1], mass_fit=options.massfit, freeze_poly=options.freeze, quiet=True)
    rhalphabuilder.run()

    input_file.Close()

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
    parser.add_option('--jet_type', dest='jet_type', help='AK8 or CA15')
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
