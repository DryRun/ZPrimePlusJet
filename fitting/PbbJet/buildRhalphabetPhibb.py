#!/usr/bin/env python
import ROOT as r, sys, math, os
from ROOT import TFile, TTree, TChain, gPad, gDirectory
from multiprocessing import Process
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array

#r.gSystem.Load("~/Dropbox/RazorAnalyzer/python/lib/libRazorRun2.so")
r.gSystem.Load(os.getenv('CMSSW_BASE') + '/lib/' + os.getenv('SCRAM_ARCH') + '/libHiggsAnalysisCombinedLimit.so')

# including other directories
# sys.path.insert(0, '../.')
from tools import *
from hist import *
from rhalphabet_builder_Phibb import RhalphabetBuilder, LoadHistograms, GetSF

#MASS_BINS = 23
#MASS_LO = 40
#MASS_HI = 201
MASS_BINS = 80 
MASS_LO = 40
MASS_HI = 600
BLIND_LO = 110
BLIND_HI = 131
#RHO_LO = -6
#RHO_HI = -2.1

BB_SF = 0.91
BB_SF_ERR = 0.03
V_SF = 0.993
V_SF_ERR = 0.043
def GetSF(process, cut, cat, f, fLoose=None, removeUnmatched=False, iPt=-1):    
    SF = 1
    if 'DMSbb' in process:
        if '50' in process:
            SF *= 0.8 * 1.574e-02
        elif '100' in process:
            SF *= 0.8 * 1.526e-02
        elif '125' in process:
            SF *= 0.8 * 1.486e-02
        elif '200' in process:
            SF *= 0.8 * 1.359e-02
        elif '300' in process:
            SF *= 0.8 * 1.251e-02
        elif '350' in process:
            SF *= 0.8 * 1.275e-02
        elif '400' in process:
            SF *= 0.8 * 1.144e-02
        elif '500' in process:
            SF *= 0.8 * 7.274e-03

    if 'hqq' in process or 'zqq' in process or 'Pbb' in process or 'Sbb' in process:
        if 'pass' in cat:
            SF *= BB_SF
            if 'zqq' in process:
                print BB_SF
        else:
            passInt = f.Get(process + '_' + cut + '_pass').Integral()
            failInt = f.Get(process + '_' + cut + '_fail').Integral()
            if failInt > 0:
                SF *= (1. + (1. - BB_SF) * passInt / failInt)
                if 'zqq' in process:
                    print (1. + (1. - BB_SF) * passInt / failInt)
    if 'wqq' in process or 'zqq' in process or 'hqq' in process or 'Pbb' in process or 'Sbb' in process:
        SF *= V_SF
        if 'zqq' in process:
            print V_SF
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
        SF *= passInt/passIntLoose
        if 'zqq' in process:
            print passInt/passIntLoose
    # remove cross section from MH=125 signal templates (template normalized to luminosity*efficiency*acceptance)
    ## if process=='hqq125':
    ##     SF *= 1./48.85*5.824E-01
    ## elif process=='zhqq':
    ##     SF *= 1./(8.839E-01*(1.-3.*0.0335962-0.201030)*5.824E-01+8.839E-01*5.824E-01*0.201030+1.227E-01*5.824E-01*0.201030+1.227E-01*5.824E-01*0.201030)
    ## elif process=='whqq':
    ##     SF *= 1./(5.328E-01*(1.-3.*0.108535)*5.824E-01+8.400E-01*(1.-3.*0.108535)*5.824E-01)
    ## elif process=='vbfhqq':
    ##     SF *= 1./(3.782*5.824E-01)
    ## elif process=='tthqq':
    ##     SF *= 1./(5.071E-01*5.824E-01)
    
    #if 'zqq' in process:
    #    print SF
    #    sys.exit()
    return SF

def LoadHistograms(f, pseudo, blind, useQCD, scale, r_signal, mass_range, blind_range, rho_range, fLoose=None, cut='p9'):

    pass_hists = {}
    fail_hists = {}
    f.ls()
    # backgrounds
    pass_hists_bkg = {}
    fail_hists_bkg = {}
    background_names = ["wqq", "zqq", "qcd", "tqq", "hqq125", "whqq125", "zhqq125", "vbfhqq125", "tthqq125"]
    for i, bkg in enumerate(background_names):
        if bkg=='qcd':
            qcd_fail = f.Get('qcd_' + cut +'_fail')
            qcd_fail.Scale(1. / scale)
            qcd_fail.SetBinContent(13, 4, (
            qcd_fail.GetBinContent(12, 4) + qcd_fail.GetBinContent(14, 4)) / 2.)  # REMOVE HIGH WEIGHT EVENT BIN
            qcd_fail.SetBinError(13, 4, (
            qcd_fail.GetBinError(12, 4) + qcd_fail.GetBinError(14, 4)) / 2.)  # REMOVE HIGH WEIGHT EVENT BIN
            if useQCD:
                qcd_pass = f.Get('qcd_' + cut +'_pass').Clone()
                qcd_pass.Scale(1. / scale)
            else:
                qcd_pass_real = f.Get('qcd_' + cut +'_pass').Clone('qcd_pass_real')
                qcd_pass_real.Scale(1. / scale)
                qcd_pass = qcd_fail.Clone('qcd_pass')
                qcd_pass_real_integral = 0
                qcd_fail_integral = 0
                for i in range(1, qcd_pass_real.GetNbinsX() + 1):
                    for j in range(1, qcd_pass_real.GetNbinsY() + 1):
                        if qcd_pass_real.GetXaxis().GetBinCenter(i) > mass_range[0] and qcd_pass_real.GetXaxis().GetBinCenter(
                                i) < mass_range[1]:
                            qcd_pass_real_integral += qcd_pass_real.GetBinContent(i, j)
                            qcd_fail_integral += qcd_fail.GetBinContent(i, j)
                qcd_pass.Scale(qcd_pass_real_integral / qcd_fail_integral)  # qcd_pass = qcd_fail * eff(pass)/eff(fail)
            pass_hists_bkg["qcd"] = qcd_pass
            fail_hists_bkg["qcd"] = qcd_fail
            print 'qcd pass integral', qcd_pass.Integral()
            print 'qcd fail integral', qcd_fail.Integral()
        elif (fLoose is not None) and (bkg=='wqq' or bkg=='zqq'):
            hpass_tmp = fLoose.Get(bkg + '_' + cut + '_pass').Clone()
            hfail_tmp = f.Get(bkg + '_' + cut + '_fail').Clone()
            hpass_tmp.Scale(1. / scale)
            hfail_tmp.Scale(1. / scale)
            hpass_tmp.Scale(GetSF(bkg, cut, 'pass', f, fLoose))
            hfail_tmp.Scale(GetSF(bkg, cut, 'fail', f))
            pass_hists_bkg[bkg] = hpass_tmp
            fail_hists_bkg[bkg] = hfail_tmp            
        else: 
            print "Attempting to get {}".format(bkg + '_' + cut + '_pass')
            hpass_tmp = f.Get(bkg + '_' + cut + '_pass').Clone()
            hfail_tmp = f.Get(bkg + '_' + cut + '_fail').Clone()
            hpass_tmp.Scale(1. / scale)
            hfail_tmp.Scale(1. / scale)
            hpass_tmp.Scale(GetSF(bkg, cut, 'pass', f))
            hfail_tmp.Scale(GetSF(bkg, cut, 'fail', f))
            pass_hists_bkg[bkg] = hpass_tmp
            fail_hists_bkg[bkg] = hfail_tmp

    # signals
    pass_hists_sig = {}
    fail_hists_sig = {}
    # for Pbb
    #masses = [50, 75, 125, 100, 150, 250, 300]
    #sigs = ['Pbb_']
    #signal_names = []
    # for Hbb
    masses = [50,100,125,200,300,350,400,500]
    #sigs = ["hqq", "zhqq", "whqq", "vbfhqq", "tthqq"]
    sigs = ["DMSbb"]
    signal_names = []
    for mass in masses:
        for sig in sigs:
            print "[debug] Getting " + sig + str(mass) + '_' + cut + "_pass"                
            passhist = f.Get(sig + str(mass) + "_" + cut + "_pass").Clone()
            failhist = f.Get(sig + str(mass) + "_" + cut + "_fail").Clone()
            for hist in [passhist, failhist]:
                for i in range(0, hist.GetNbinsX() + 2):
                    for j in range(0, hist.GetNbinsY() + 2):
                        if hist.GetBinContent(i, j) <= 0:
                            hist.SetBinContent(i, j, 0)
            failhist.Scale(1. / scale)
            passhist.Scale(1. / scale)
            failhist.Scale(GetSF(sig + str(mass), cut, 'fail', f))
            passhist.Scale(GetSF(sig + str(mass), cut, 'pass', f))
            pass_hists_sig[sig + str(mass)] = passhist
            fail_hists_sig[sig + str(mass)] = failhist
            signal_names.append(sig + str(mass))

    # Systematics
    all_systematics = []
    pass_hists_syst = {}
    fail_hists_syst = {}
    for syst in ['JES', 'JER', 'trigger', 'mcstat','Pu']:
        all_systematics.append(syst)
        for direction in ["Up", "Down"]:
            syst_dir = syst + direction
            pass_hists_syst[syst_dir] = {}
            fail_hists_syst[syst_dir] = {}
            for process in ["tqq", "wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"] + signal_names:
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
                    this_pass_hist.Scale(GetSF(signal_name, 'pass', jet_type, interpolation_file))
                    this_fail_hist.Scale(GetSF(signal_name, 'fail', jet_type, interpolation_file))
                else:
                    this_pass_hist.Scale(GetSF(signal_name, 'pass', jet_type, input_file))
                    this_fail_hist.Scale(GetSF(signal_name, 'fail', jet_type, input_file))
                pass_hists_syst[syst_dir][process] = this_pass_hist
                fail_hists_syst[syst_dir][process] = this_fail_hist
    # mcstat systematic
    # - This produces one histogram for up, and one histogram for down, with all bins varied coherently.
    # - RhalphabetBuilder is responsible for "uncorrelating" the bins, i.e. make a separate uncertainty for each bin.
    all_systematics.append("mcstat")
    for direction in ["Up", "Down"]:
        pass_hists_syst["mcstat{}".format(direction)] = {}
        fail_hists_syst["mcstat{}".format(direction)] = {}
        for process in ["tqq", "wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"] + signal_names:
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

        for process in ["wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"] + signal_names:
            if process == 'wqq':
                mass = 80.385
            elif process == 'zqq':
                mass = 91.1876
            elif process in ["hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"]:
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

                    pass_hists_syst["ptscaleUp"][process].SetBinContent(massbin, ptbin, hist_syst_shift[0].GetBinContent(massbin, ptbin))
                    pass_hists_syst["ptscaleUp"][process].SetBinError(massbin, ptbin, hist_syst_shift[0].GetBinError(massbin, ptbin))
                    pass_hists_syst["ptscaleDown"][process].SetBinContent(massbin, ptbin, hist_syst_shift[1].GetBinContent(massbin, ptbin))
                    pass_hists_syst["ptscaleDown"][process].SetBinError(massbin, ptbin, hist_syst_shift[1].GetBinError(massbin, ptbin))

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

                    fail_hists_syst["ptscaleUp"][process].SetBinContent(massbin, ptbin, hist_syst_shift[0].GetBinContent(massbin, ptbin))
                    fail_hists_syst["ptscaleUp"][process].SetBinError(massbin, ptbin, hist_syst_shift[0].GetBinError(massbin, ptbin))
                    fail_hists_syst["ptscaleDown"][process].SetBinContent(massbin, ptbin, hist_syst_shift[1].GetBinContent(massbin, ptbin))
                    fail_hists_syst["ptscaleDown"][process].SetBinError(massbin, ptbin, hist_syst_shift[1].GetBinError(massbin, ptbin))

                    fail_hists_syst["smearUp"][process].SetBinContent(massbin, ptbin, hist_syst_smear[0].GetBinContent(massbin, ptbin))
                    fail_hists_syst["smearUp"][process].SetBinError(massbin, ptbin, hist_syst_smear[0].GetBinError(massbin, ptbin))
                    fail_hists_syst["smearDown"][process].SetBinContent(massbin, ptbin, hist_syst_smear[1].GetBinContent(massbin, ptbin))
                    fail_hists_syst["smearDown"][process].SetBinError(massbin, ptbin, hist_syst_smear[1].GetBinError(massbin, ptbin))
                # End loop: massbin
            # End loop: ptbin
        # End loop: process
    # End if do_shift

    if pseudo:        
        for i, bkg in enumerate(background_names):
            if i==0:
                pass_hists["data_obs"] = pass_hists_bkg[bkg].Clone('data_obs_pass')
                fail_hists["data_obs"] = fail_hists_bkg[bkg].Clone('data_obs_fail')
            else:                
                pass_hists["data_obs"].Add(pass_hists_bkg[bkg])
                fail_hists["data_obs"].Add(fail_hists_bkg[bkg])

        for i, signal in enumerate(signal_names):
            pass_hists["data_obs"].Add(pass_hists_sig[signal],r_signal)
            fail_hists["data_obs"].Add(fail_hists_sig[signal],r_signal)
    else:
        pass_hists["data_obs"] = f.Get('data_obs_pass')
        fail_hists["data_obs"] = f.Get('data_obs_fail')
    pass_hists.update(pass_hists_bkg)
    pass_hists.update(pass_hists_sig)
    fail_hists.update(fail_hists_bkg)
    fail_hists.update(fail_hists_sig)
    
    for histogram in (pass_hists.values() + fail_hists.values()):
        for i in range(1,histogram.GetNbinsX()+1):
            for j in range(1,histogram.GetNbinsY()+1):
                massVal = histogram.GetXaxis().GetBinCenter(i)
                ptVal = histogram.GetYaxis().GetBinLowEdge(j) + histogram.GetYaxis().GetBinWidth(j) * 0.3
                rhoVal = r.TMath.Log(massVal * massVal / ptVal / ptVal)
                if blind and histogram.GetXaxis().GetBinCenter(i) > blind_range[0] and histogram.GetXaxis().GetBinCenter(i) < blind_range[1]:
                    print "blinding signal region for %s, mass bin [%i,%i] "%(histogram.GetName(),histogram.GetXaxis().GetBinLowEdge(i),histogram.GetXaxis().GetBinUpEdge(i))
                    histogram.SetBinContent(i,j,0.)
                if rhoVal < rho_range[0] or rhoVal > rho_range[1]:
                    print "removing rho = %.2f for %s, pt bin [%i, %i], mass bin [%i,%i]" % (
                    rhoVal, histogram.GetName(), histogram.GetYaxis().GetBinLowEdge(j), histogram.GetYaxis().GetBinUpEdge(j),
                    histogram.GetXaxis().GetBinLowEdge(i), histogram.GetXaxis().GetBinUpEdge(i))
                    histogram.SetBinContent(i, j, 0.)
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
    return (pass_hists, fail_hists, pass_hists_syst, fail_hists_syst, all_systematics)


def main(options, args):
    ifile = options.ifile
    odir = options.odir

    # Load the input histograms
    # 	- 2D histograms of pass and fail mass,pT distributions
    # 	- for each MC sample and the data
    f = r.TFile.Open(ifile)    
    fLoose = None
    if options.ifile_loose is not None:
        fLoose = r.TFile.Open(options.ifile_loose)
    #(hpass, hfail) = loadHistograms(f, options.pseudo, options.blind, options.useQCD, options.scale, options.r)
    (pass_hists, fail_hists,pass_hists_syst, fail_hists_syst, all_systematics) = LoadHistograms(f, options.pseudo, options.blind, options.useQCD, scale=options.scale, r_signal=options.r, mass_range=[MASS_LO, MASS_HI], blind_range=[BLIND_LO, BLIND_HI], rho_range=[options.lrho, options.hrho], fLoose=fLoose, cuts = options.cuts)
    #f.Close()

    # Build the workspacees
    #dazsleRhalphabetBuilder(hpass, hfail, f, odir, options.NR, options.NP)

    rhalphabuilder = RhalphabetBuilder(pass_hists, fail_hists, f, options.odir, nr=options.NR, np=options.NP, mass_nbins=MASS_BINS, mass_lo=MASS_LO, mass_hi=MASS_HI, blind_lo=BLIND_LO, blind_hi=BLIND_HI, rho_lo=options.lrho, rho_hi=options.hrho, blind=options.blind, mass_fit=options.massfit, freeze_poly=options.freeze, remove_unmatched=options.removeUnmatched, input_file_loose=fLoose, cuts=options.cuts)
    rhalphabuilder.add_systematics(all_systematics, pass_hists_syst, fail_hists_syst)
    rhalphabuilder.run()
    if options.prefit:
        rhalphabuilder.prefit()
    elif options.loadfit is not None:
        rhalphabuilder.loadfit(options.loadfit)
        

##-------------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option('-i', '--ifile', dest='ifile', default='hist_1DZbb.root', help='file with histogram inputs',
                      metavar='ifile')
    parser.add_option('--ifile-loose', dest='ifile_loose', default=None, help='second file with histogram inputs (looser b-tag cut to take W/Z/H templates)',
                      metavar='ifile_loose')
    parser.add_option('-o', '--odir', dest='odir', default='./', help='directory to write plots', metavar='odir')
    parser.add_option('--pseudo', action='store_true', dest='pseudo', default=False, help='use MC', metavar='pseudo')
    parser.add_option('--blind', action='store_true', dest='blind', default=False, help='blind signal region',
                      metavar='blind')
    parser.add_option('--use-qcd', action='store_true', dest='useQCD', default=False, help='use real QCD MC',
                      metavar='useQCD')
    parser.add_option('--massfit', action='store_true', dest='massfit', default=False, help='mass fit or rho',
                      metavar='massfit')
    parser.add_option('--freeze', action='store_true', dest='freeze', default=False, help='freeze pol values',
                      metavar='freeze')
    parser.add_option('--scale', dest='scale', default=1, type='float',
                      help='scale factor to scale MC (assuming only using a fraction of the data)')
    parser.add_option('--nr', dest='NR', default=2, type='int', help='order of rho (or mass) polynomial')
    parser.add_option('--np', dest='NP', default=1, type='int', help='order of pt polynomial')
    parser.add_option('-r', dest='r', default=0, type='float', help='signal strength for MC pseudodataset')
    parser.add_option('--remove-unmatched', action='store_true', dest='removeUnmatched', default =False,help='remove unmatched', metavar='removeUnmatched')
    parser.add_option('--prefit', action='store_true', dest='prefit', default =False,help='do prefit', metavar='prefit')
    parser.add_option('--loadfit', dest='loadfit', default=None, help='load qcd polynomial parameters from alternative rhalphabase.root',metavar='loadfit')
    parser.add_option('--lrho', dest='lrho', default=-6.0, type= 'float', help='low value rho cut')
    parser.add_option('--hrho', dest='hrho', default=-2.1, type='float', help=' high value rho cut')

    parser.add_option('-c', '--cut', dest='cut', default='p9', type='string', help='double b-tag cut value')

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

    main(options, args)
##-------------------------------------------------------------------------------------
