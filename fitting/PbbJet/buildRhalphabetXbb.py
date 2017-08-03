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
from DAZSLE.ZPrimePlusJet.rhalphabet_builder import RhalphabetBuilder

gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libDAZSLEPhiBBPlusJet.so"))
import DAZSLE.ZPrimePlusJet.xbb_config as config

# Scale factors for MC
# - Double b-tag SF from muon control region
# - Vector(bb) sf... what is this again?
def GetSF(process, cat, jet_type, f, fLoose=None, removeUnmatched=False, iPt=-1):
    SF = 1.

    # bb SF, for MC process with real bb
    if 'hbb' in process or 'zqq' in process or 'Pbb' in process or 'Sbb' in process:
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
    if 'wqq' in process or 'zqq' in process or 'hbb' in process or 'Pbb' in process or 'Sbb' in process:
        SF *= config.analysis_parameters[jet_type]["V_SF"]
    return SF


# Load histograms here, rather than rhalphabet_builder.LoadHistograms. I think it's unreasonable to expect that different analyses will share this function!
# input_file_path = path of TFile containing the pass/fail msd-pT histograms
# mass_range      = range of mSD to use
# rho_range       = range of rho to use
# pseudo          = Substitute pseudodata constructed from MC for real data
# useQCD = For QCD pass, use MC prediction instead of fail * (pass int / fail int)
# fLoose = 
def LoadHistograms(input_file, interpolation_file, mass_range, rho_range, jet_type=None, pseudo=False, useQCD=False, loose_file_path=None, scale=1., r_signal=0., blind_range=None, do_shift=True):
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
            hpass_tmp.Scale(GetSF(bkg, 'pass', jet_type, input_file, input_file_loose))
            hfail_tmp.Scale(GetSF(bkg, 'fail', jet_type, input_file))
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
            hpass_tmp.Scale(GetSF(bkg, 'pass', jet_type, input_file))
            hfail_tmp.Scale(GetSF(bkg, 'fail', jet_type, input_file))
            pass_hists_bkg[bkg] = hpass_tmp
            fail_hists_bkg[bkg] = hfail_tmp

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
            failhist.Scale(GetSF(signal_name, 'fail', jet_type, input_file))
            passhist.Scale(GetSF(signal_name, 'pass', jet_type, input_file))

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
            failhist.Scale(GetSF(signal_name, 'fail', jet_type, interpolation_file))
            passhist.Scale(GetSF(signal_name, 'pass', jet_type, interpolation_file))

            pass_hists_sig[signal_name] = passhist
            fail_hists_sig[signal_name] = failhist

    # Systematics
    all_systematics = []
    # - Histogram-based
    pass_hists_syst = {}
    fail_hists_syst = {}
    for syst in ['JES', 'JER', 'trigger','Pu']:
        all_systematics.append(syst)
        for direction in ["Up", "Down"]:
            syst_dir = syst + direction
            pass_hists_syst[syst_dir] = {}
            fail_hists_syst[syst_dir] = {}
            for process in ["tqq", "wqq", "zqq", "hbb"] + config.limit_signal_names[jet_type]:
                if process in config.interpolated_signal_names:
                    this_pass_hist = interpolation_file.Get(process + "_pass_" + syst_dir)
                    this_fail_hist = interpolation_file.Get(process + "_fail_" + syst_dir)
                else:
                    this_pass_hist = input_file.Get(process + "_pass_" + syst_dir)
                    this_fail_hist = input_file.Get(process + "_fail_" + syst_dir)
                # Apply manual scale factor (reminder: this is for adjusting to using a fraction of the full dataset)
                pass_hists_syst[syst_dir][process].Scale(1. / scale)
                fail_hists_syst[syst_dir][process].Scale(1. / scale)

                # Apply other SFs
                if process in config.interpolated_signal_names:
                    pass_hists_syst[syst_dir][process].Scale(GetSF(signal_name, 'pass', jet_type, interpolation_file))
                    fail_hists_syst[syst_dir][process].Scale(GetSF(signal_name, 'fail', jet_type, interpolation_file))
                else:
                    pass_hists_syst[syst_dir][process].Scale(GetSF(signal_name, 'pass', jet_type, input_file))
                    fail_hists_syst[syst_dir][process].Scale(GetSF(signal_name, 'fail', jet_type, input_file))

    # mcstat systematic
    # - This produces one histogram for up, and one histogram for down, with all bins varied coherently.
    # - RhalphabetBuilder is responsible for "uncorrelating" the bins, i.e. make a separate uncertainty for each bin.
    all_systematics.append("mcstat")
    for direction in ["Up", "Down"]:
        pass_hists_syst["mcstat{}".format(direction)] = {}
        fail_hists_syst["mcstat{}".format(direction)] = {}
        for process in ["tqq", "wqq", "zqq", "hbb"] + config.limit_signal_names[jet_type]:
            if process in config.interpolated_signal_names:
                this_pass_hist = interpolation_file.Get(process + "_pass")
                this_fail_hist = interpolation_file.Get(process + "_fail")
            else:
                this_pass_hist = input_file.Get(process + "_pass")
                this_fail_hist = input_file.Get(process + "_fail")
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

    # Scale/shift systematics
    # - The central value scale/shift is also done here. 
    if do_shift:
        all_systematics.append("ptscale")
        all_systematics.append("smear")
        m_data = 82.657
        m_data_err = 0.313
        m_mc = 82.548
        m_mc_err = 0.191
        s_data = 8.701
        s_data_err = 0.433
        s_mc = 8.027
        s_mc_err = 0.607
        mass_shift = m_data / m_mc
        mass_shift_unc = math.sqrt((m_data_err / m_data)**2 + (m_mc_err / m_mc)**2) * 10.  # (10 sigma shift)
        res_shift = s_data / s_mc
        res_shift_unc = math.sqrt((s_data_err / s_data) * (s_data_err / s_data) + (s_mc_err / s_mc) * (
            s_mc_err / s_mc)) * 2.  # (2 sigma shift)

        pass_hists_syst["ptscaleUp"] = {}
        pass_hists_syst["ptscaleDown"] = {}
        pass_hists_syst["smearUp"] = {}
        pass_hists_syst["smearDown"] = {}
        fail_hists_syst["ptscaleUp"] = {}
        fail_hists_syst["ptscaleDown"] = {}
        fail_hists_syst["smearUp"] = {}
        fail_hists_syst["smearDown"] = {}

        for process in ["wqq", "zqq", "hbb"] + config.limit_signal_names[jet_type]:
            if process == 'wqq':
                mass = 80.385
            elif process == 'zqq':
                mass = 91.1876
            elif 'hbb' in process:
                mass = 125
            elif 'Sbb' in process:
                re_match = re_sbb.search(process)
                mass = int(re_match.group("mass"))

            # Start from the existing central value histograms
            # - Note: this code assumes that the histograms are already V-matched, and that we're ignoring the unmatched part
            this_pass_hist = pass_hists[process]
            pass_hists_syst["ptscaleUp"][process] = this_pass_hist.Clone()
            pass_hists_syst["ptscaleUp"][process].SetName(this_pass_hist.GetName() + "_ptscaleUp")
            pass_hists_syst["ptscaleDown"][process] = this_pass_hist.Clone()
            pass_hists_syst["ptscaleDown"][process].SetName(this_pass_hist.GetName() + "_ptscaleDown")
            pass_hists_syst["smearUp"][process] = this_pass_hist.Clone()
            pass_hists_syst["smearUp"][process].SetName(this_pass_hist.GetName() + "_smearUp")
            pass_hists_syst["smearDown"][process] = this_pass_hist.Clone()
            pass_hists_syst["smearDown"][process].SetName(this_pass_hist.GetName() + "_smearDown")

            for ptbin in xrange(1, this_pass_hist.GetNbinsY() + 1):
                this_hist_ptbin = this_pass_hist.ProjectionX(this_pass_hist.GetName() + "_ptbin{}".format(ptbin), ptbin, ptbin)
                hist_container = hist([mass], [this_hist_ptbin])
                shift_val = mass - mass * mass_shift
                tmp_shifted_h = hist_container.shift(tmph_mass_matched, shift_val)
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

                    pass_hists_syst["ptscaleUp"][process].SetBinContent(massbin, ptbin, hist_syst_scale[0].GetBinContent(massbin, ptbin))
                    pass_hists_syst["ptscaleUp"][process].SetBinError(massbin, ptbin, hist_syst_scale[0].GetBinError(massbin, ptbin))
                    pass_hists_syst["ptscaleDown"][process].SetBinContent(massbin, ptbin, hist_syst_scale[1].GetBinContent(massbin, ptbin))
                    pass_hists_syst["ptscaleDown"][process].SetBinError(massbin, ptbin, hist_syst_scale[1].GetBinError(massbin, ptbin))

                    pass_hists_syst["smearUp"][process].SetBinContent(massbin, ptbin, hist_syst_smear[0].GetBinContent(massbin, ptbin))
                    pass_hists_syst["smearUp"][process].SetBinError(massbin, ptbin, hist_syst_smear[0].GetBinError(massbin, ptbin))
                    pass_hists_syst["smearDown"][process].SetBinContent(massbin, ptbin, hist_syst_smear[1].GetBinContent(massbin, ptbin))
                    pass_hists_syst["smearDown"][process].SetBinError(massbin, ptbin, hist_syst_smear[1].GetBinError(massbin, ptbin))
                # End loop: massbin
            # End loop: ptbin

            this_fail_hist = fail_hists[process]
            fail_hists_syst["ptscaleUp"][process] = this_fail_hist.Clone()
            fail_hists_syst["ptscaleUp"][process].SetName(this_fail_hist.GetName() + "_ptscaleUp")
            fail_hists_syst["ptscaleDown"][process] = this_fail_hist.Clone()
            fail_hists_syst["ptscaleDown"][process].SetName(this_fail_hist.GetName() + "_ptscaleDown")
            fail_hists_syst["smearUp"][process] = this_fail_hist.Clone()
            fail_hists_syst["smearUp"][process].SetName(this_fail_hist.GetName() + "_smearUp")
            fail_hists_syst["smearDown"][process] = this_fail_hist.Clone()
            fail_hists_syst["smearDown"][process].SetName(this_fail_hist.GetName() + "_smearDown")

            for ptbin in xrange(1, this_fail_hist.GetNbinsY() + 1):
                this_hist_ptbin = this_fail_hist.ProjectionX(this_fail_hist.GetName() + "_ptbin{}".format(ptbin), ptbin, ptbin)
                hist_container = hist([mass], [this_hist_ptbin])
                shift_val = mass - mass * mass_shift
                tmp_shifted_h = hist_container.shift(tmph_mass_matched, shift_val)
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

                    fail_hists_syst["ptscaleUp"][process].SetBinContent(massbin, ptbin, hist_syst_scale[0].GetBinContent(massbin, ptbin))
                    fail_hists_syst["ptscaleUp"][process].SetBinError(massbin, ptbin, hist_syst_scale[0].GetBinError(massbin, ptbin))
                    fail_hists_syst["ptscaleDown"][process].SetBinContent(massbin, ptbin, hist_syst_scale[1].GetBinContent(massbin, ptbin))
                    fail_hists_syst["ptscaleDown"][process].SetBinError(massbin, ptbin, hist_syst_scale[1].GetBinError(massbin, ptbin))

                    fail_hists_syst["smearUp"][process].SetBinContent(massbin, ptbin, hist_syst_smear[0].GetBinContent(massbin, ptbin))
                    fail_hists_syst["smearUp"][process].SetBinError(massbin, ptbin, hist_syst_smear[0].GetBinError(massbin, ptbin))
                    fail_hists_syst["smearDown"][process].SetBinContent(massbin, ptbin, hist_syst_smear[1].GetBinContent(massbin, ptbin))
                    fail_hists_syst["smearDown"][process].SetBinError(massbin, ptbin, hist_syst_smear[1].GetBinError(massbin, ptbin))
                # End loop: massbin
            # End loop: ptbin
        # End loop: process
    # End if do_shift

    # Load data
    # - Either real data, or pseudodata=sum of background MC
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
    else:
        pass_hists["data_obs"] = input_file.Get('data_obs_pass')
        fail_hists["data_obs"] = input_file.Get('data_obs_fail')

    pass_hists.update(pass_hists_bkg)
    pass_hists.update(pass_hists_sig)
    fail_hists.update(fail_hists_bkg)
    fail_hists.update(fail_hists_sig)

    print pass_hists
    print fail_hists

    # Postprocessing:
    # - Blind signal bins, if requested
    # - Set bins outside the desired rho range to zero
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
                if rhoVal < rho_range[0] or rhoVal > rho_range[1]:
                    #print "removing rho = %.2f for %s, pt bin [%i, %i], mass bin [%i,%i]" % (
                    #    rhoVal, histogram.GetName(), histogram.GetYaxis().GetBinLowEdge(j),
                    #    histogram.GetYaxis().GetBinUpEdge(j),
                    #    histogram.GetXaxis().GetBinLowEdge(i), histogram.GetXaxis().GetBinUpEdge(i))
                    histogram.SetBinContent(i, j, 0.)
        histogram.SetDirectory(0)
    for syst in all_systematics:
        for process in pass_hists_syst[syst].keys():
            for histogram in [pass_hists_syst[syst][process], fail_hists_syst[syst][process]]:
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
                            #print "removing rho = %.2f for %s, pt bin [%i, %i], mass bin [%i,%i]" % (
                            #    rhoVal, histogram.GetName(), histogram.GetYaxis().GetBinLowEdge(j),
                            #    histogram.GetYaxis().GetBinUpEdge(j),
                            #    histogram.GetXaxis().GetBinLowEdge(i), histogram.GetXaxis().GetBinUpEdge(i))
                            histogram.SetBinContent(i, j, 0.)
                histogram.SetDirectory(0)


        # print "lengths = ", len(pass_hists), len(fail_hists)
    # print pass_hists;
    # print fail_hists;
    return (pass_hists, fail_hists, pass_hists_syst, fail_hists_syst)


def main(options, args):
    input_file = TFile(config.get_histogram_file("SR", options.jet_type), "READ")
    interpolation_file = TFile(config.get_interpolation_file(options.jet_type), "READ")
    odir = config.paths["LimitSetting"] + "/Xbb_inputs/"

    # Load the input histograms
    # 	- 2D histograms of pass and fail mass,pT distributions
    # 	- for each MC sample and the data
    (pass_hists,fail_hists, pass_hists_syst, fail_hists_syst) = LoadHistograms(input_file, interpolation_file, config.analysis_parameters[options.jet_type]["MSD"], config.analysis_parameters[options.jet_type]["RHO"], useQCD=options.useQCD, jet_type=options.jet_type, scale=options.scale, r_signal=options.r, pseudo=options.pseudo, do_shift=True)

    rhalphabuilder = RhalphabetBuilder(pass_hists, fail_hists, input_file, odir, nr=options.NR, np=options.NP, mass_nbins=config.analysis_parameters[options.jet_type]["MASS_BINS"], mass_lo=config.analysis_parameters[options.jet_type]["MSD"][0], mass_hi=config.analysis_parameters[options.jet_type]["MSD"][1], rho_lo=config.analysis_parameters[options.jet_type]["RHO"][0], rho_hi=config.analysis_parameters[options.jet_type]["RHO"][1], mass_fit=options.massfit, freeze_poly=options.freeze, quiet=True, signal_names=config.limit_signal_names[options.jet_type], interpolation_file=interpolation_file)
    rhalphabuilder.add_systematics(all_systematics, pass_hists_syst, fail_hists_syst)
    rhalphabuilder.run()
    if options.addHptShape:
        rhalphabuilder.addHptShape()    
    if options.prefit:
        rhalphabuilder.prefit()
    elif options.loadfit is not None:
        rhalphabuilder.loadfit(options.loadfit)

    # Copy outputs to subdirectories
    for signal_name in config.limit_signal_names[options.jet_type]:
        os.system("cp {}/base.root {}".format(odir, config.get_datacard_directory(signal_name, options.jet_type, qcd=options.pseudo)))
        os.system("cp {}/rhalphabase.root {}".format(odir, config.get_datacard_directory(signal_name, options.jet_type, qcd=options.pseudo)))

    input_file.Close()
    interpolation_file.Close()

##-------------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = OptionParser()
    #parser.add_option('-i', '--ifile', dest='ifile', default='hist_1DZbb.root', help='file with histogram inputs', metavar='ifile')
    #parser.add_option('-o', '--odir', dest='odir', default='./', help='directory to write plots', metavar='odir')
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
    parser.add_option('--np', dest='NP', default=1, type='int', help='order of pt polynomial')
    parser.add_option('--jet_type', dest='jet_type', help='AK8 or CA15')
    parser.add_option('--prefit', action='store_true', dest='prefit', default =False,help='do prefit', metavar='prefit')
    parser.add_option('-r', dest='r', default=0, type='float', help='signal strength for MC pseudodataset')
    (options, args) = parser.parse_args()

    main(options, args)
##-------------------------------------------------------------------------------------
