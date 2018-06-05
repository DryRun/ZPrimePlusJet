#!/usr/bin/env python
import ROOT
import sys
import math
import os
import re
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
from histogram_container import HistogramContainer
from DAZSLE.ZPrimePlusJet.rhalphabet_builder import RhalphabetBuilder

gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libDAZSLEPhiBBPlusJet.so"))
import DAZSLE.ZPrimePlusJet.xbb_config as config

re_sbb = re.compile("Sbb(?P<mass>\d+)")

# Scale factors for MC
# - Double b-tag SF from muon control region
# - Vector(bb) sf... what is this again?
# - n2ddt_fail: SF is different for N2 fail
def GetSF(process, cat, jet_type, f, fLoose=None, removeUnmatched=False, iPt=-1, region=""):
    if region == "SR":
        return GetSF_SR(process, cat, jet_type, f, fLoose=fLoose, removeUnmatched=removeUnmatched, iPt=-1)
    elif region == "muCR":
        return GetSF_SR(process, cat, jet_type, f, fLoose=fLoose, removeUnmatched=removeUnmatched, iPt=-1)
    elif region == "N2SR":
        return GetSF_N2SR(process, cat, jet_type, f, fLoose=fLoose, removeUnmatched=removeUnmatched, iPt=-1)
    elif region == "N2CR":
        return GetSF_N2CR(process, cat, jet_type, f, fLoose=fLoose, removeUnmatched=removeUnmatched, iPt=-1)
    else:
        print "[GetSF] ERROR : Unknown region {}".format(region)

def GetSF_SR(process, cat, jet_type, f, fLoose=None, removeUnmatched=False, iPt=-1):
    SF = 1.

    # bb SF, for MC process with real bb
    if process in ["zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125", "tqq"] or "Pbb" in process or "Sbb" in process:
        if 'pass' in cat:
            SF *= config.analysis_parameters[jet_type]["BB_SF"]
            if 'zqq' in process:
                print config.analysis_parameters[jet_type]["BB_SF"]
        else:
            passInt = f.Get(process + '_pass').Integral()
            failInt = f.Get(process + '_fail').Integral()
            if failInt > 0:
                SF *= (1. + (1. - config.analysis_parameters[jet_type]["BB_SF"]) * passInt / failInt)
                if 'zqq' in process:
                    print (1. + (1. - config.analysis_parameters[jet_type]["BB_SF"]) * passInt / failInt)

    # V SF, for data-MC agreement for substructure cut
    if process in ["wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"] or "Pbb" in process or "Sbb" in process:
        SF *= config.analysis_parameters[jet_type]["V_SF"]
    return SF

# N2SR: Need to apply N2 pass/fail SF, and doubleB pass only.
def GetSF_N2SR(process, cat, jet_type, f, fLoose=None, removeUnmatched=False, iPt=-1):
    SF = 1.

    # bb SF, for MC process with real bb
    if process in ["zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125", "tqq"] or "Pbb" in process or "Sbb" in process:
        SF *= config.analysis_parameters[jet_type]["BB_SF"]
        if 'zqq' in process:
            print config.analysis_parameters[jet_type]["BB_SF"]

    # V SF, for data-MC agreement for substructure cut
    if process in ["wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"] or "Pbb" in process or "Sbb" in process:
        if "pass" in cat:
            SF *= config.analysis_parameters[jet_type]["V_SF"]
        elif "fail" in cat:
            SF *= config.analysis_parameters[jet_type]["VFAIL_SF"]
        else:
            print "[GetSFN2SR] ERROR : cat {} doesn't have pass or fail in it!".format(cat)
    return SF

# N2CR: Need to apply N2/fail SF
def GetSF_N2SR(process, cat, jet_type, f, fLoose=None, removeUnmatched=False, iPt=-1):
    SF = 1.

    # bb SF, for MC process with real bb
    if process in ["zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125", "tqq"] or "Pbb" in process or "Sbb" in process:
        if 'pass' in cat:
            SF *= config.analysis_parameters[jet_type]["BB_SF"]
            if 'zqq' in process:
                print config.analysis_parameters[jet_type]["BB_SF"]
        else:
            passInt = f.Get(process + '_pass').Integral()
            failInt = f.Get(process + '_fail').Integral()
            if failInt > 0:
                SF *= (1. + (1. - config.analysis_parameters[jet_type]["BB_SF"]) * passInt / failInt)
                if 'zqq' in process:
                    print (1. + (1. - config.analysis_parameters[jet_type]["BB_SF"]) * passInt / failInt)

    # V SF, for data-MC agreement for substructure cut
    if process in ["wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"] or "Pbb" in process or "Sbb" in process:
        SF *= config.analysis_parameters[jet_type]["VFAIL_SF"]
    return SF

# muCR: I think this is the same as SR. 
def GetSF_muCR(process, cat, jet_type, f, fLoose=None, removeUnmatched=False, iPt=-1):
    return GetSF_SR(process, cat, jet_type, f, fLoose=fLoose, removeUnmatched=removeUnmatched, iPt=iPt)


# Load histograms here, rather than rhalphabet_builder.LoadHistograms. I think it's unreasonable to expect that different analyses will share this function.
# input_file_path = path of TFile containing the pass/fail msd-pT histograms
# mass_range      = range of mSD to use
# rho_range       = range of rho to use
# pseudo          = Substitute pseudodata constructed from MC for real data
# useQCD = For QCD pass, use MC prediction instead of fail * (pass int / fail int)
# fLoose = 
def LoadHistograms(input_file, interpolation_file, mass_range, rho_range, jet_type=None, pseudo=False, useQCD=False, decidata=False, loose_file_path=None, scale=1., r_signal=0., blind_range=None, do_shift=True, region="", min_mass=-1.):
    pass_hists = {}
    fail_hists = {}
    #f.ls()

    # backgrounds
    pass_hists_bkg = {}
    fail_hists_bkg = {}
    background_names = ["wqq", "zqq", "qcd", "tqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"]
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
            hpass_tmp.Scale(GetSF(bkg, 'pass', jet_type, input_file, input_file_loose, region=region))
            hfail_tmp.Scale(GetSF(bkg, 'fail', jet_type, input_file, region=region))
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
            hpass_tmp.Scale(GetSF(bkg, 'pass', jet_type, input_file, region=region))
            hfail_tmp.Scale(GetSF(bkg, 'fail', jet_type, input_file, region=region))
            pass_hists_bkg[bkg] = hpass_tmp
            fail_hists_bkg[bkg] = hfail_tmp

        # If decidata, scale MC histograms by 0.1
        if decidata:
            pass_hists_bkg[bkg].Scale(0.1)
            fail_hists_bkg[bkg].Scale(0.1)

    # signals
    pass_hists_sig = {}
    fail_hists_sig = {}

    for signal_name in config.limit_signal_names[jet_type]:
        print "[debug] Getting " + signal_name + "_pass"
        if signal_name in config.simulated_signal_names:
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
            failhist.Scale(GetSF(signal_name, 'fail', jet_type, input_file, region=region))
            passhist.Scale(GetSF(signal_name, 'pass', jet_type, input_file, region=region))

            pass_hists_sig[signal_name] = passhist
            fail_hists_sig[signal_name] = failhist
            #signal_names.append(signal_name)
        elif signal_name in config.interpolated_signal_names:
            if not interpolation_file.Get(signal_name + "_pass"):
                print "[LoadHistograms] ERROR : Couldn't load histogram " + signal_name + "_pass" + " from file " + input_file.GetPath()
                sys.exit(1)
            passhist = interpolation_file.Get(signal_name + "_pass").Clone()
            failhist = interpolation_file.Get(signal_name + "_fail").Clone()
            for hist in [passhist, failhist]:
                for i in range(0, hist.GetNbinsX() + 2):
                    for j in range(0, hist.GetNbinsY() + 2):
                        if hist.GetBinContent(i, j) <= 0:
                            hist.SetBinContent(i, j, 0)
            failhist.Scale(1. / scale)
            passhist.Scale(1. / scale)
            failhist.Scale(GetSF(signal_name, 'fail', jet_type, interpolation_file, region=region))
            passhist.Scale(GetSF(signal_name, 'pass', jet_type, interpolation_file, region=region))

            pass_hists_sig[signal_name] = passhist
            fail_hists_sig[signal_name] = failhist

        if decidata:
            pass_hists_sig[signal_name].Scale(0.1)
            fail_hists_sig[signal_name].Scale(0.1)
    pass_hists.update(pass_hists_bkg)
    pass_hists.update(pass_hists_sig)
    fail_hists.update(fail_hists_bkg)
    fail_hists.update(fail_hists_sig)

    # Systematics
    all_systematics = []
    # - Histogram-based
    pass_hists_syst = {}
    fail_hists_syst = {}
    for syst in ['JES', 'JER', 'Trigger','PU']:
        all_systematics.append(syst)
        for direction in ["Up", "Down"]:
            syst_dir = syst + direction
            pass_hists_syst[syst_dir] = {}
            fail_hists_syst[syst_dir] = {}
            for process in ["tqq", "wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"] + config.limit_signal_names[jet_type]:
                if process in config.interpolated_signal_names:
                    this_pass_hist = interpolation_file.Get(process + "_pass_" + syst_dir)
                    this_fail_hist = interpolation_file.Get(process + "_fail_" + syst_dir)
                else:
                    this_pass_hist = input_file.Get(process + "_pass_" + syst_dir)
                    this_fail_hist = input_file.Get(process + "_fail_" + syst_dir)
                    if not this_pass_hist:
                        print "[LoadHistograms] ERROR : Histogram {} not found in file {}".format(process + "_pass_" + syst_dir, input_file.GetPath())
                # Apply manual scale factor (reminder: this is for adjusting to using a fraction of the full dataset)
                this_pass_hist.Scale(1. / scale)
                this_fail_hist.Scale(1. / scale)

                # Apply other SFs
                if process in config.interpolated_signal_names:
                    this_pass_hist.Scale(GetSF(signal_name, 'pass', jet_type, interpolation_file, region=region))
                    this_fail_hist.Scale(GetSF(signal_name, 'fail', jet_type, interpolation_file, region=region))
                else:
                    this_pass_hist.Scale(GetSF(signal_name, 'pass', jet_type, input_file, region=region))
                    this_fail_hist.Scale(GetSF(signal_name, 'fail', jet_type, input_file, region=region))
                pass_hists_syst[syst_dir][process] = this_pass_hist
                fail_hists_syst[syst_dir][process] = this_fail_hist

                if decidata:
                    pass_hists_syst[syst_dir][process].Scale(0.1)
                    fail_hists_syst[syst_dir][process].Scale(0.1)

    # mcstat systematic
    # - This produces one histogram for up, and one histogram for down, with all bins varied coherently.
    # - RhalphabetBuilder is responsible for "uncorrelating" the bins, i.e. make a separate uncertainty for each bin.
    all_systematics.append("mcstat")
    for direction in ["Up", "Down"]:
        pass_hists_syst["mcstat{}".format(direction)] = {}
        fail_hists_syst["mcstat{}".format(direction)] = {}
        for process in ["tqq", "wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"] + config.limit_signal_names[jet_type]:
            if process in config.interpolated_signal_names:
                this_pass_hist = interpolation_file.Get(process + "_pass").Clone()
                this_fail_hist = interpolation_file.Get(process + "_fail").Clone()
            else:
                this_pass_hist = input_file.Get(process + "_pass").Clone()
                this_fail_hist = input_file.Get(process + "_fail").Clone()
            this_pass_hist.SetName(process + "_pass_mcstat{}".format(direction))
            this_fail_hist.SetName(process + "_fail_mcstat{}".format(direction))
            for xbin in xrange(0, this_pass_hist.GetNbinsX() + 2):
                for ybin in xrange(0, this_pass_hist.GetNbinsY() + 2):
                    if direction == "Up":
                        alpha = 1.
                    else:
                        alpha = -1.
                    this_pass_hist.SetBinContent(xbin, ybin, this_pass_hist.GetBinContent(xbin, ybin) + (alpha * this_pass_hist.GetBinError(xbin, ybin)))
                    this_fail_hist.SetBinContent(xbin, ybin, this_fail_hist.GetBinContent(xbin, ybin) + (alpha * this_fail_hist.GetBinError(xbin, ybin)))
            pass_hists_syst["mcstat{}".format(direction)][process] = this_pass_hist
            fail_hists_syst["mcstat{}".format(direction)][process] = this_fail_hist

            if decidata:
                pass_hists_syst["mcstat{}".format(direction)][process].Scale(0.1)
                fail_hists_syst["mcstat{}".format(direction)][process].Scale(0.1)

    # Scale/shift systematics
    # - The central value scale/shift is also done here. 
    if do_shift:
        all_systematics.append("scalept")
        all_systematics.append("smear")

        mass_shift = config.analysis_parameters[jet_type]["MASS_SF"]
        mass_shift_unc = config.analysis_parameters[jet_type]["MASS_SF_ERR"] * 10.# (10 sigma shift) 
        res_shift = config.analysis_parameters[jet_type]["RES_SF"]
        res_shift_unc = config.analysis_parameters[jet_type]["RES_SF_ERR"] * 2.  # (2 sigma shift)

        pass_hists_syst["scaleptUp"] = {}
        pass_hists_syst["scaleptDown"] = {}
        pass_hists_syst["smearUp"] = {}
        pass_hists_syst["smearDown"] = {}
        fail_hists_syst["scaleptUp"] = {}
        fail_hists_syst["scaleptDown"] = {}
        fail_hists_syst["smearUp"] = {}
        fail_hists_syst["smearDown"] = {}

        for process in ["wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"] + config.limit_signal_names[jet_type]:
            if process == 'wqq':
                mass = 80.385
            elif process == 'zqq':
                mass = 91.1876
            elif process in ["hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"]:
                mass = 125
            elif 'Sbb' in process:
                re_match = re_sbb.search(process)
                mass = int(re_match.group("mass"))
            else:
                print "Unrecognized process in scale/shift: {}".format(process)
                sys.exit(1)

            # Start from the existing central value histograms
            # - Note: this code assumes that the histograms are already V-matched, and that we're ignoring the unmatched part
            this_pass_hist = pass_hists[process]
            pass_hists_syst["scaleptUp"][process] = this_pass_hist.Clone()
            pass_hists_syst["scaleptUp"][process].SetName(this_pass_hist.GetName() + "_scaleptUp")
            pass_hists_syst["scaleptDown"][process] = this_pass_hist.Clone()
            pass_hists_syst["scaleptDown"][process].SetName(this_pass_hist.GetName() + "_scaleptDown")
            pass_hists_syst["smearUp"][process] = this_pass_hist.Clone()
            pass_hists_syst["smearUp"][process].SetName(this_pass_hist.GetName() + "_smearUp")
            pass_hists_syst["smearDown"][process] = this_pass_hist.Clone()
            pass_hists_syst["smearDown"][process].SetName(this_pass_hist.GetName() + "_smearDown")

            for ptbin in xrange(1, this_pass_hist.GetNbinsY() + 1):
                this_hist_ptbin = this_pass_hist.ProjectionX(this_pass_hist.GetName() + "_ptbin{}".format(ptbin), ptbin, ptbin)
                hist_container = HistogramContainer([mass], [this_hist_ptbin])
                shift_val = mass - mass * mass_shift
                tmp_shifted_h = hist_container.shift(this_hist_ptbin, shift_val)
                # get new central value and new smeared value
                smear_val = res_shift - 1
                tmp_smeared_h = hist_container.smear(tmp_shifted_h[0], smear_val)
                if smear_val <= 0: 
                    hist_new_central = tmp_smeared_h[1]
                else:
                    hist_new_central = tmp_smeared_h[0]

                # get shift up/down
                hist_syst_shift = hist_container.shift(hist_new_central, mass * mass_shift_unc)
                # get res up/down
                hist_syst_smear = hist_container.smear(hist_new_central, res_shift_unc)

                # Set new bin contents for this pt bin
                for massbin in xrange(1, this_pass_hist.GetNbinsX() + 1):
                    pass_hists[process].SetBinContent(massbin, ptbin, hist_new_central.GetBinContent(massbin, ptbin))
                    pass_hists[process].SetBinError(massbin, ptbin, hist_new_central.GetBinError(massbin, ptbin))

                    pass_hists_syst["scaleptUp"][process].SetBinContent(massbin, ptbin, hist_syst_shift[0].GetBinContent(massbin, ptbin))
                    pass_hists_syst["scaleptUp"][process].SetBinError(massbin, ptbin, hist_syst_shift[0].GetBinError(massbin, ptbin))
                    pass_hists_syst["scaleptDown"][process].SetBinContent(massbin, ptbin, hist_syst_shift[1].GetBinContent(massbin, ptbin))
                    pass_hists_syst["scaleptDown"][process].SetBinError(massbin, ptbin, hist_syst_shift[1].GetBinError(massbin, ptbin))

                    pass_hists_syst["smearUp"][process].SetBinContent(massbin, ptbin, hist_syst_smear[0].GetBinContent(massbin, ptbin))
                    pass_hists_syst["smearUp"][process].SetBinError(massbin, ptbin, hist_syst_smear[0].GetBinError(massbin, ptbin))
                    pass_hists_syst["smearDown"][process].SetBinContent(massbin, ptbin, hist_syst_smear[1].GetBinContent(massbin, ptbin))
                    pass_hists_syst["smearDown"][process].SetBinError(massbin, ptbin, hist_syst_smear[1].GetBinError(massbin, ptbin))
                # End loop: massbin
            # End loop: ptbin

            this_fail_hist = fail_hists[process]
            fail_hists_syst["scaleptUp"][process] = this_fail_hist.Clone()
            fail_hists_syst["scaleptUp"][process].SetName(this_fail_hist.GetName() + "_scaleptUp")
            fail_hists_syst["scaleptDown"][process] = this_fail_hist.Clone()
            fail_hists_syst["scaleptDown"][process].SetName(this_fail_hist.GetName() + "_scaleptDown")
            fail_hists_syst["smearUp"][process] = this_fail_hist.Clone()
            fail_hists_syst["smearUp"][process].SetName(this_fail_hist.GetName() + "_smearUp")
            fail_hists_syst["smearDown"][process] = this_fail_hist.Clone()
            fail_hists_syst["smearDown"][process].SetName(this_fail_hist.GetName() + "_smearDown")

            for ptbin in xrange(1, this_fail_hist.GetNbinsY() + 1):
                this_hist_ptbin = this_fail_hist.ProjectionX(this_fail_hist.GetName() + "_ptbin{}".format(ptbin), ptbin, ptbin)
                hist_container = HistogramContainer([mass], [this_hist_ptbin])
                shift_val = mass - mass * mass_shift
                tmp_shifted_h = hist_container.shift(this_hist_ptbin, shift_val)
                # get new central value and new smeared value
                smear_val = res_shift - 1
                tmp_smeared_h = hist_container.smear(tmp_shifted_h[0], smear_val)
                if smear_val <= 0: 
                    hist_new_central = tmp_smeared_h[1]
                else:
                    hist_new_central = tmp_smeared_h[0]

                # get shift up/down
                hist_syst_shift = hist_container.shift(hist_new_central, mass * mass_shift_unc)
                # get res up/down
                hist_syst_smear = hist_container.smear(hist_new_central, res_shift_unc)

                # Set new bin contents for this pt bin
                for massbin in xrange(1, this_fail_hist.GetNbinsX() + 1):
                    fail_hists[process].SetBinContent(massbin, ptbin, hist_new_central.GetBinContent(massbin, ptbin))
                    fail_hists[process].SetBinError(massbin, ptbin, hist_new_central.GetBinError(massbin, ptbin))

                    fail_hists_syst["scaleptUp"][process].SetBinContent(massbin, ptbin, hist_syst_shift[0].GetBinContent(massbin, ptbin))
                    fail_hists_syst["scaleptUp"][process].SetBinError(massbin, ptbin, hist_syst_shift[0].GetBinError(massbin, ptbin))
                    fail_hists_syst["scaleptDown"][process].SetBinContent(massbin, ptbin, hist_syst_shift[1].GetBinContent(massbin, ptbin))
                    fail_hists_syst["scaleptDown"][process].SetBinError(massbin, ptbin, hist_syst_shift[1].GetBinError(massbin, ptbin))

                    fail_hists_syst["smearUp"][process].SetBinContent(massbin, ptbin, hist_syst_smear[0].GetBinContent(massbin, ptbin))
                    fail_hists_syst["smearUp"][process].SetBinError(massbin, ptbin, hist_syst_smear[0].GetBinError(massbin, ptbin))
                    fail_hists_syst["smearDown"][process].SetBinContent(massbin, ptbin, hist_syst_smear[1].GetBinContent(massbin, ptbin))
                    fail_hists_syst["smearDown"][process].SetBinError(massbin, ptbin, hist_syst_smear[1].GetBinError(massbin, ptbin))
                # End loop: massbin
            # End loop: ptbin
        # End loop: process
    # End if do_shift

    # Load data
    # - Either real data, or pseudodata=sum of background MC, or data PS 10
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

        for i, signal_name in enumerate(config.limit_signal_names[jet_type]):
            pass_hists["data_obs"].Add(pass_hists_sig[signal_name], r_signal)
            fail_hists["data_obs"].Add(fail_hists_sig[signal_name], r_signal)
    elif decidata:
        pass_hists["data_obs"] = input_file.Get('data_obs_ps10_pass')
        if not pass_hists["data_obs"]:
            print "[LoadHistograms] ERROR : Couldn't find histogram data_obs_ps10_pass in file {}".format(input_file.GetPath())
        pass_hists["data_obs"].SetName("data_obs_pass")
        fail_hists["data_obs"] = input_file.Get('data_obs_ps10_fail')
        fail_hists["data_obs"].SetName("data_obs_fail")
    else:
        pass_hists["data_obs"] = input_file.Get('data_obs_pass')
        fail_hists["data_obs"] = input_file.Get('data_obs_fail')

    print pass_hists
    print fail_hists

    # Postprocessing:
    # - Blind signal bins, if requested
    # - Set bins outside the desired rho range to zero
    # - Set bins outside the mSD range to zero
    # - SetDirectory(0)
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
                        histogram.SetBinError(i, j, 0.)
                if rhoVal < rho_range[0] or rhoVal > rho_range[1]:
                    #print "removing rho = %.2f for %s, pt bin [%i, %i], mass bin [%i,%i]" % (
                    #    rhoVal, histogram.GetName(), histogram.GetYaxis().GetBinLowEdge(j),
                    #    histogram.GetYaxis().GetBinUpEdge(j),
                    #    histogram.GetXaxis().GetBinLowEdge(i), histogram.GetXaxis().GetBinUpEdge(i))
                    if "data_obs" in histogram.GetName():
                        print "[debug] ZERO: Setting hist {}, bin mSD={}/rho={}/pt={} to zero because of rho cut".format(histogram.GetName(), massVal, rhoVal, ptVal)
                    histogram.SetBinContent(i, j, 0.)
                    histogram.SetBinError(i, j, 0.)
                if massVal < mass_range[0] or massVal > mass_range[1]:
                    if "data_obs" in histogram.GetName():
                        print "[debug] ZERO: Setting hist {}, bin mSD={}/rho={}/pt={} to zero because of mSD cut".format(histogram.GetName(), massVal, rhoVal, ptVal)
                    histogram.SetBinContent(i, j, 0.)
                    histogram.SetBinError(i, j, 0.)
                if min_mass >= 0.:
                    if massVal < min_mass:
                        if "data_obs" in histogram.GetName():
                            print "[debug] ZERO: Setting hist {}, bin mSD={}/rho={}/pt={} to zero because of min_mass cut".format(histogram.GetName(), massVal, rhoVal, ptVal)

                        histogram.SetBinContent(i, j, 0.)
                        histogram.SetBinError(i, j, 0.)
        histogram.SetDirectory(0)
    for syst in all_systematics:
        for direction in ["Up", "Down"]:
            syst_dir = "{}{}".format(syst, direction)
            for process in pass_hists_syst[syst_dir].keys():
                for histogram in [pass_hists_syst[syst_dir][process], fail_hists_syst[syst_dir][process]]:
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
                                    histogram.SetBinError(i, j, 0.)
                            if rhoVal < rho_range[0] or rhoVal > rho_range[1]:
                                #print "removing rho = %.2f for %s, pt bin [%i, %i], mass bin [%i,%i]" % (
                                #    rhoVal, histogram.GetName(), histogram.GetYaxis().GetBinLowEdge(j),
                                #    histogram.GetYaxis().GetBinUpEdge(j),
                                #    histogram.GetXaxis().GetBinLowEdge(i), histogram.GetXaxis().GetBinUpEdge(i))
                                histogram.SetBinContent(i, j, 0.)
                                histogram.SetBinError(i, j, 0.)
                            if massVal < mass_range[0] or massVal > mass_range[1]:
                                histogram.SetBinContent(i, j, 0.)
                                histogram.SetBinError(i, j, 0.)
                            if min_mass >= 0.:
                                if massVal < min_mass:
                                    histogram.SetBinContent(i, j, 0.)
                                    histogram.SetBinError(i, j, 0.)
                    histogram.SetDirectory(0)


        # print "lengths = ", len(pass_hists), len(fail_hists)
    # print pass_hists;
    # print fail_hists;
    return (pass_hists, fail_hists, pass_hists_syst, fail_hists_syst, all_systematics)


def main(options, args):
    input_file = TFile(config.get_histogram_file(options.region, options.jet_type), "READ")
    interpolation_file = TFile(config.get_interpolation_file(options.region, options.jet_type), "READ")
    #odir =  "{}/Xbb_inputs/{}_{}/".format(config.paths["LimitSetting"], options.region, options.jet_type)
    odir = config.get_datacard_directory("Sbb125", options.jet_type, qcd=options.pseudo, decidata=options.decidata, region=options.region) + "/../"

    # Load the input histograms
    # 	- 2D histograms of pass and fail mass,pT distributions
    # 	- for each MC sample and the data
    if options.jet_type == "CA15" and options.region == "N2SR":
        print "[debug] Setting mass range to [68, 600] for N2SR"
        mass_range = [68., 600.]
    else:
        mass_range = config.analysis_parameters[options.jet_type]["MSD"]
    (pass_hists,fail_hists, pass_hists_syst, fail_hists_syst, all_systematics) = LoadHistograms(input_file, interpolation_file, mass_range, config.analysis_parameters[options.jet_type]["RHO"], useQCD=options.useQCD, jet_type=options.jet_type, scale=options.scale, r_signal=options.r, pseudo=options.pseudo, do_shift=True, decidata=options.decidata, region=options.region)

    # MSD bins: restricted=mass_nbins=config.analysis_parameters[options.jet_type]["MASS_BINS"], mass_lo=config.analysis_parameters[options.jet_type]["MSD"][0], mass_hi=config.analysis_parameters[options.jet_type]["MSD"][1], full = 80, 40, 600
    rhalphabuilder = RhalphabetBuilder(pass_hists, fail_hists, input_file, odir, nr=config.analysis_parameters[options.jet_type]["MAX_NRHO"], np=config.analysis_parameters[options.jet_type]["MAX_NPT"], mass_nbins=80, mass_lo=40, mass_hi=600, rho_lo=config.analysis_parameters[options.jet_type]["RHO"][0], rho_hi=config.analysis_parameters[options.jet_type]["RHO"][1], mass_fit=options.massfit, freeze_poly=options.freeze, quiet=True, signal_names=config.limit_signal_names[options.jet_type])
    rhalphabuilder.add_systematics(all_systematics, pass_hists_syst, fail_hists_syst)
    rhalphabuilder.run()
    if options.addHptShape:
        rhalphabuilder.addHptShape()    
    if options.prefit:
        # List of poly coeffs to freeze
        fix_pars_rhalphabet = {}
        for irho in xrange(0, config.analysis_parameters[options.jet_type]["MAX_NRHO"]):
            for ipt in xrange(0, config.analysis_parameters[options.jet_type]["MAX_NPT"]):
                if irho ==0 and ipt == 0:
                    continue
                fix_pars_rhalphabet["r{}p{}".format(irho, ipt)] = 0.
        cats = config.analysis_parameters[options.jet_type]["FIT_PT_BINS"]
        #if options.region == "N2SR":
        #    cats = [1,2,3,4,5,6]
        rhalphabuilder.prefit(fix_pars_rhalphabet=fix_pars_rhalphabet, category_indices=cats)

    # Copy outputs to subdirectories
    for signal_name in config.limit_signal_names[options.jet_type]:
        os.system("mkdir -pv {}".format(config.get_datacard_directory(signal_name, options.jet_type, qcd=options.pseudo, decidata=options.decidata, region=options.region)))
        #print "cp {}/base.root {}".format(odir, config.get_datacard_directory(signal_name, options.jet_type, qcd=options.pseudo, decidata=options.decidata, region=options.region))
        #os.system("cp {}/base.root {}".format(odir, config.get_datacard_directory(signal_name, options.jet_type, qcd=options.pseudo, decidata=options.decidata, region=options.region)))
        #os.system("cp {}/rhalphabase.root {}".format(odir, config.get_datacard_directory(signal_name, options.jet_type, qcd=options.pseudo, decidata=options.decidata, region=options.region)))
        os.system("mv {}/base.root {}/base.root.old".format(config.get_datacard_directory(signal_name, options.jet_type, qcd=options.pseudo, decidata=options.decidata, region=options.region), config.get_datacard_directory(signal_name, options.jet_type, qcd=options.pseudo, decidata=options.decidata, region=options.region)))
        os.system("mv {}/rhalphabase.root {}/rhalphabase.root.old".format(config.get_datacard_directory(signal_name, options.jet_type, qcd=options.pseudo, decidata=options.decidata, region=options.region), config.get_datacard_directory(signal_name, options.jet_type, qcd=options.pseudo, decidata=options.decidata, region=options.region)))
        print "ln -s {}/base.root {}/base.root".format(config.get_datacard_directory(signal_name, options.jet_type, qcd=options.pseudo, decidata=options.decidata, region=options.region), config.get_datacard_directory(signal_name, options.jet_type, qcd=options.pseudo, decidata=options.decidata, region=options.region))
        os.system("ln -s {}/base.root {}/base.root".format(odir, config.get_datacard_directory(signal_name, options.jet_type, qcd=options.pseudo, decidata=options.decidata, region=options.region)))
        os.system("ln -s {}/rhalphabase.root {}/rhalphabase.root".format(odir, config.get_datacard_directory(signal_name, options.jet_type, qcd=options.pseudo, decidata=options.decidata, region=options.region)))

    input_file.Close()
    interpolation_file.Close()

##-------------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = OptionParser()
    #parser.add_option('-i', '--ifile', dest='ifile', default='hist_1DZbb.root', help='file with histogram inputs', metavar='ifile')
    #parser.add_option('-o', '--odir', dest='odir', default='./', help='directory to write plots', metavar='odir')
    parser.add_option('--pseudo', action='store_true', dest='pseudo', default=False, help='use MC', metavar='pseudo')
    parser.add_option('--decidata', action='store_true', dest='decidata', default=False, help='Use 1/10 data', metavar='decidata')
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
    #parser.add_option('--nr', dest='NR', type='int', help='order of rho (or mass) polynomial')
    #parser.add_option('--np', dest='NP', type='int', help='order of pt polynomial')
    parser.add_option('--jet_type', dest='jet_type', help='AK8 or CA15')
    parser.add_option('--prefit', action='store_true', dest='prefit', default =False,help='do prefit', metavar='prefit')
    parser.add_option('-r', dest='r', default=0, type='float', help='signal strength for MC pseudodataset')
    parser.add_option('--addHptShape',action='store_true',dest='addHptShape',default =False,help='add H pt shape unc', metavar='addHptShape')
    parser.add_option('--region', type=str, default="SR", help="SR, N2SR, or N2CR", dest="region")
    (options, args) = parser.parse_args()

    main(options, args)
    print "Done."
##-------------------------------------------------------------------------------------
