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
import tools as tools

##############################################################################
##############################################################################
#### B E G I N N I N G   O F   C L A S S
##############################################################################
##############################################################################

class RhalphabetBuilder(): 
    def __init__( self, pass_hists, fail_hists, odir, mass_nbins=35, mass_lo=40, mass_hi=285, blind_lo=110, blind_hi=131, mass_fit=False, freeze_poly=False): 
        self._pass_hists = pass_hists;
        self._fail_hists = fail_hists;
        self._mass_fit = mass_fit
        self._freeze = freeze_poly

        self._output_path = "{}/base.root".format(odir)
        self._rhalphabet_output_path = "{}/rhalphabase.root".format(odir)

        self._mass_nbins = mass_nbins
        self._mass_lo    = mass_lo
        self._mass_hi    = mass_hi
        self._mass_blind_lo = blind_lo
        self._mass_blind_hi = blind_hi
        # self._mass_nbins = pass_hists[0].GetXaxis().GetNbins();
        # self._mass_lo    = pass_hists[0].GetXaxis().GetBinLowEdge( 1 );
        # self._mass_hi    = pass_hists[0].GetXaxis().GetBinUpEdge( self._mass_nbins );

        print "number of mass bins and lo/hi: ", self._mass_nbins, self._mass_lo, self._mass_hi;

        #polynomial order for fit
        self._poly_degree_pt = 2 #1 = linear ; 2 is quadratic
        self._poly_degree_rho = 2
        #self._poly_degree_rhoP =1;

        self._nptbins = pass_hists["data_obs"].GetYaxis().GetNbins();
        self._pt_lo = pass_hists["data_obs"].GetYaxis().GetBinLowEdge( 1 );
        self._pt_hi = pass_hists["data_obs"].GetYaxis().GetBinUpEdge( self._nptbins );

        # define RooRealVars
        self._lMSD    = r.RooRealVar("x","x",self._mass_lo,self._mass_hi)
        self._lMSD.setRange('Low',self._mass_lo,self._mass_blind_lo)
        self._lMSD.setRange('Blind',self._mass_blind_lo,self._mass_blind_hi)
        self._lMSD.setRange('High',self._mass_blind_hi,self._mass_hi)
        #self._lMSD.setBins(self._mass_nbins)       
        self._lPt     = r.RooRealVar("pt","pt",self._pt_lo,self._pt_hi);
        self._lPt.setBins(self._nptbins)
        self._lRho    = r.RooFormulaVar("rho","log(x*x/pt/pt)",r.RooArgList(self._lMSD,self._lPt))

        self._lEff    = r.RooRealVar("veff"      ,"veff"      ,0.5 ,0.,1.0)

        self._lEffQCD = r.RooRealVar("qcdeff"    ,"qcdeff"   ,0.01,0.,10.)
        qcd_pass_integral = 0
        qcd_fail_integral = 0
        for i in range(1,fail_hists["qcd"].GetNbinsX()+1):
            for j in range(1,fail_hists["qcd"].GetNbinsY()+1):
                if fail_hists["qcd"].GetXaxis().GetBinCenter(i) > self._mass_lo and fail_hists["qcd"].GetXaxis().GetBinCenter(i) < self._mass_hi:
                    qcd_fail_integral += fail_hists["qcd"].GetBinContent(i,j)
                    qcd_pass_integral += pass_hists["qcd"].GetBinContent(i,j)
        if qcd_fail_integral>0:
            qcdeff = qcd_pass_integral / qcd_fail_integral
            self._lEffQCD.setVal(qcdeff)
        print "qcdeff = %f"%qcdeff
        self._lDM     = r.RooRealVar("dm","dm", 0.,-10,10)
        self._lShift  = r.RooFormulaVar("shift",self._lMSD.GetName()+"-dm",r.RooArgList(self._lMSD,self._lDM)) 

        self._all_vars = [];
        self._all_shapes = [];
        self._all_data = [];
        self._all_pars = [];

        self._background_names = ["wqq", "zqq", "qcd", "tqq"]
        self._signal_names = []
        for mass in [50,75,125,100,150,250,300]:
            self._signal_names.append("Pbb_" + str(mass))

        self.LoopOverPtBins();

    def LoopOverPtBins(self):

        print "number of pt bins = ", self._nptbins;
        for pt_bin in range(1,self._nptbins+1):
        # for pt_bin in range(1,2):
            print "------- pT bin number ",pt_bin      

            # 1d histograms in each pT bin (in the order... data, w, z, qcd, top, signals)
            pass_hists_ptbin = {};
            fail_hists_ptbin = {};
            for name, hist in self._pass_hists.iteritems():
                pass_hists_ptbin[name] = tools.proj("cat",str(pt_bin),hist,self._mass_nbins,self._mass_lo,self._mass_hi)
            for name, hist in self._fail_hists.iteritems():
                fail_hists_ptbin[name] = tools.proj("cat",str(pt_bin),hist,self._mass_nbins,self._mass_lo,self._mass_hi)

            # make RooDataset, RooPdfs, and histograms
            # GetWorkspaceInputs returns: RooDataHist (data), then RooHistPdf of each electroweak
            (data_pass_rdh, data_fail_rdh, pass_rhps, fail_rhps) = self.GetWorkspaceInputs(pass_hists_ptbin, fail_hists_ptbin,"cat"+str(pt_bin))
            #Get approximate pt bin value
            this_pt = self._pass_hists["data_obs"].GetYaxis().GetBinLowEdge(pt_bin)+self._pass_hists["data_obs"].GetYaxis().GetBinWidth(pt_bin)*0.3;
            print "------- this bin pT value ",this_pt
            
            #Make the rhalphabet fit for this pt bin
            (rhalphabet_hist_pass, rhalphabet_hist_fail) = self.MakeRhalphabet(["data_obs", "wqq", "zqq", "tqq"], fail_hists_ptbin, this_pt, "cat"+str(pt_bin))

            # Get signals
            (signal_rdhs_pass, signal_rdhs_fail) = self.GetSignalInputs(pass_hists_ptbin, fail_hists_ptbin, "cat"+str(pt_bin))
            pass_rhps.update(signal_rdhs_pass)
            fail_rhps.update(signal_rdhs_fail)

            # #Write to file
            print "pass_rhps = "
            print pass_rhps
            self.MakeWorkspace(self._output_path, [data_pass_rdh] + pass_rhps.values(), "pass_cat"+str(pt_bin), True)
            self.MakeWorkspace(self._output_path, [data_fail_rdh] + fail_rhps.values(), "fail_cat"+str(pt_bin), True)

        
        for pt_bin in range(1,self._nptbins+1):
            for mass_bin in range(1,self._mass_nbins+1):
                print "qcd_fail_cat%i_Bin%i flatParam" % (pt_bin, mass_bin);

    # iHs = dict of fail histograms
    def MakeRhalphabet(self, samples, fail_histograms, pt, category):
        print "---- [MakeRhalphabet]"   

        rhalph_bkgd_name ="qcd";
        lUnity = r.RooConstVar("unity","unity",1.)
        lZero  = r.RooConstVar("lZero","lZero",0.)

        #Fix the pt (top) and the qcd eff
        self._lPt.setVal(pt)
        self._lEffQCD.setConstant(False)

        polynomial_variables = []
        self.buildPolynomialArray(polynomial_variables,self._poly_degree_pt,self._poly_degree_rho,"p","r",-0.3,0.3)
        print "polynomial_variables=",
        print polynomial_variables

        #Now build the function
        pass_bins = r.RooArgList()
        fail_bins = r.RooArgList()

        for mass_bin in range(1,self._mass_nbins+1):
            self._lMSD.setVal(fail_histograms["data_obs"].GetXaxis().GetBinCenter(mass_bin))
            if self._mass_fit :
                print ("Pt/mass poly")
                roopolyarray = self.buildRooPolyArray(self._lPt.getVal(),self._lMSD.getVal(),lUnity,lZero,polynomial_variables)
            else :
                print ("Pt/Rho poly")
                roopolyarray = self.buildRooPolyRhoArray(self._lPt.getVal(),self._lRho.getVal(),lUnity,lZero,polynomial_variables)
            print "RooPolyArray:"
            roopolyarray.Print()
            fail_bin_content = 0
            for sample in samples:
                if sample == "data_obs":
                    print sample, fail_histograms[sample].GetName(), "add data"
                    print "\t+={}".format(fail_histograms[sample].GetBinContent(mass_bin))
                    fail_bin_content += fail_histograms[sample].GetBinContent(mass_bin) # add data
                else:                    
                    print sample, fail_histograms[sample].GetName(), "subtract W/Z/ttbar"
                    print "\t-={}".format(fail_histograms[sample].GetBinContent(mass_bin))
                    fail_bin_content -= fail_histograms[sample].GetBinContent(mass_bin) # subtract W/Z/ttbar from data
            if fail_bin_content < 0: fail_bin_content = 0

            print rhalph_bkgd_name+"_fail_"+category+"_Bin"+str(mass_bin), fail_bin_content

            #5 sigma range + 10 events
            fail_bin_unc = math.sqrt(fail_bin_content)*5+10
            #Define the failing category
            #fail_bin_var = r.RooRealVar(rhalph_bkgd_name+"_fail_"+category+"_Bin"+str(mass_bin),rhalph_bkgd_name+"_fail_"+category+"_Bin"+str(mass_bin),fail_bin_content,max(fail_bin_content-fail_bin_unc,0),max(fail_bin_content+fail_bin_unc,0))
            fail_bin_var = r.RooRealVar(rhalph_bkgd_name+"_fail_"+category+"_Bin"+str(mass_bin),rhalph_bkgd_name+"_fail_"+category+"_Bin"+str(mass_bin),fail_bin_content,0,max(fail_bin_content+fail_bin_unc,0))

            print "[david debug] fail_bin_var:"
            fail_bin_var.Print()

            #Now define the passing cateogry based on the failing (make sure it can't go negative)
            lArg = r.RooArgList(fail_bin_var,roopolyarray,self._lEffQCD)
            pass_bin_var = r.RooFormulaVar(rhalph_bkgd_name+"_pass_"+category+"_Bin"+str(mass_bin),rhalph_bkgd_name+"_pass_"+category+"_Bin"+str(mass_bin),"@0*max(@1,0)*@2",lArg)
            print "Pass=fail*poly*eff RooFormulaVar:"
            print pass_bin_var.Print()

            # print pass_bin_var.GetName()

            #If the number of events in the failing is small remove the bin from being free in the fit
            if fail_bin_content < 4:
                print "too small number of events", fail_bin_content, "Bin", str(mass_bin)
                fail_bin_var.setConstant(True)
                pass_bin_var = r.RooRealVar(rhalph_bkgd_name+"_pass_"+category+"_Bin"+str(mass_bin),rhalph_bkgd_name+"_pass_"+category+"_Bin"+str(mass_bin),0,0,0)
                pass_bin_var.setConstant(True)

            #Add bins to the array
            pass_bins.add(pass_bin_var)
            fail_bins.add(fail_bin_var)
            self._all_vars.extend([pass_bin_var,fail_bin_var])
            self._all_pars.extend([pass_bin_var,fail_bin_var])
            # print  fail_bin_var.GetName(),"flatParam",lPass#,lPass+"/("+lFail+")*@0"

        #print "Printing pass_bins:"
        #for i in xrange(pass_bins.getSize()):
        #    pass_bins[i].Print()
        pass_rparh  = r.RooParametricHist(rhalph_bkgd_name+"_pass_"+category,rhalph_bkgd_name+"_pass_"+category,self._lMSD,pass_bins,fail_histograms["data_obs"])
        fail_rparh  = r.RooParametricHist(rhalph_bkgd_name+"_fail_"+category,rhalph_bkgd_name+"_fail_"+category,self._lMSD,fail_bins,fail_histograms["data_obs"])
        print "Print pass and fail RooParametricHists"
        pass_rparh.Print()
        fail_rparh.Print()
        pass_norm = r.RooAddition(rhalph_bkgd_name+"_pass_"+category+"_norm",rhalph_bkgd_name+"_pass_"+category+"_norm",pass_bins)
        fail_norm = r.RooAddition(rhalph_bkgd_name+"_fail_"+category+"_norm",rhalph_bkgd_name+"_fail_"+category+"_norm",fail_bins)
        print "Printing NPass and NFail variables:"
        pass_norm.Print()
        fail_norm.Print()        
        self._all_shapes.extend([pass_rparh,fail_rparh,pass_norm,fail_norm])

        #Now write the wrokspace with the rooparamhist
        pass_workspace = r.RooWorkspace("w_pass_"+str(category))
        fail_workspace = r.RooWorkspace("w_fail_"+str(category))
        getattr(pass_workspace,'import')(pass_rparh,r.RooFit.RecycleConflictNodes())
        getattr(pass_workspace,'import')(pass_norm,r.RooFit.RecycleConflictNodes())
        getattr(fail_workspace,'import')(fail_rparh,r.RooFit.RecycleConflictNodes())
        getattr(fail_workspace,'import')(fail_norm,r.RooFit.RecycleConflictNodes())
        print "Printing rhalphabet workspace:"
        pass_workspace.Print()
        if category.find("1") > -1:
            pass_workspace.writeToFile(self._rhalphabet_output_path)
        else:
            pass_workspace.writeToFile(self._rhalphabet_output_path,False)
        fail_workspace.writeToFile(self._rhalphabet_output_path,False)
        return [pass_rparh,fail_rparh]

    def buildRooPolyArray(self,iPt,iMass,iQCD,iZero,iVars):

        # print "---- [buildRooPolyArray]"  
        # print len(iVars);

        lPt  = r.RooConstVar("Var_Pt_" +str(iPt)+"_"+str(iMass), "Var_Pt_" +str(iPt)+"_"+str(iMass),(iPt))
        lMass = r.RooConstVar("Var_Mass_"+str(iPt)+"_"+str(iMass), "Var_Mass_"+str(iPt)+"_"+str(iMass),(iMass))
        lMassArray = r.RooArgList()
        lNCount=0
        for pRVar in range(0,self._poly_degree_rho+1):
            lTmpArray = r.RooArgList()
            for pVar in range(0,self._poly_degree_pt+1):
                if lNCount == 0:
                    lTmpArray.add(iQCD) # for the very first constant (e.g. p0r0), just set that to 1
                else:
                    lTmpArray.add(iVars[lNCount])
                lNCount=lNCount+1
            pLabel="Var_Pol_Bin_"+str(round(iPt,2))+"_"+str(round(iMass,3))+"_"+str(pRVar)
            pPol = r.RooPolyVar(pLabel,pLabel,lPt,lTmpArray)
            lMassArray.add(pPol)
            self._all_vars.append(pPol)

        lLabel="Var_MassPol_Bin_"+str(round(iPt,2))+"_"+str(round(iMass,3))
        lMassPol = r.RooPolyVar(lLabel,lLabel,lMass,lMassArray)
        self._all_vars.extend([lPt,lMass,lMassPol])
        return lMassPol

    def buildRooPolyRhoArray(self,iPt,iRho,iQCD,iZero,iVars):

        # print "---- [buildRooPolyArray]"      

        lPt  = r.RooConstVar("Var_Pt_" +str(iPt)+"_"+str(iRho), "Var_Pt_" +str(iPt)+"_"+str(iRho),(iPt))
        lRho = r.RooConstVar("Var_Rho_"+str(iPt)+"_"+str(iRho), "Var_Rho_"+str(iPt)+"_"+str(iRho),(iRho))
        lRhoArray = r.RooArgList()
        lNCount=0
        for pRVar in range(0,self._poly_degree_rho+1):
            lTmpArray = r.RooArgList()
            for pVar in range(0,self._poly_degree_pt+1):
                if lNCount == 0: 
                    lTmpArray.add(iQCD); # for the very first constant (e.g. p0r0), just set that to 1
                else: 
                    print "lNCount = " + str(lNCount)
                    lTmpArray.add(iVars[lNCount])
                lNCount=lNCount+1
            pLabel="Var_Pol_Bin_"+str(round(iPt,2))+"_"+str(round(iRho,3))+"_"+str(pRVar)
            pPol = r.RooPolyVar(pLabel,pLabel,lPt,lTmpArray)
            print "pPol:"
            print pPol.Print()
            lRhoArray.add(pPol);
            self._all_vars.append(pPol)

        lLabel="Var_RhoPol_Bin_"+str(round(iPt,2))+"_"+str(round(iRho,3))
        lRhoPol = r.RooPolyVar(lLabel,lLabel,lRho,lRhoArray)
        self._all_vars.extend([lPt,lRho,lRhoPol])
        return lRhoPol


    def buildPolynomialArray(self, iVars,iNVar0,iNVar1,iLabel0,iLabel1,iXMin0,iXMax0):

        print "---- [buildPolynomialArray]"
        ## form of polynomial
        ## (p0r0 + p1r0 * pT + p2r0 * pT^2 + ...) + 
        ## (p0r1 + p1r1 * pT + p2r1 * pT^2 + ...) * rho + 
        ## (p0r2 + p1r2 * pT + p2r2 * pT^2 + ...) * rho^2 + ...
        '''
        r0p0    =    0, pXMin,pXMax
        r1p0    =   -3.7215e-03 +/-  1.71e-08
        r2p0    =    2.4063e-06 +/-  2.76e-11
            r0p1    =   -2.1088e-01 +/-  2.72e-06I  
            r1p1    =    3.6847e-05 +/-  4.66e-09
            r2p1    =   -3.8415e-07 +/-  7.23e-12
            r0p2    =   -8.5276e-02 +/-  6.90e-07
            r1p2    =    2.2058e-04 +/-  1.10e-09
            r2p2    =   -2.2425e-07 +/-  1.64e-12
        '''
        value = [ 0.,
            -3.7215e-03,
            2.4063e-06,
            -2.1088e-01, 
            3.6847e-05, 
            -3.8415e-07, 
            -8.5276e-02, 
            2.2058e-04,
            -2.2425e-07]
        error = [iXMax0,
            1.71e-08,
            2.76e-11,
            2.72e-06,
            4.66e-09,
            7.23e-12,
            6.90e-07,
            1.10e-09,
            1.64e-12]
        
        for i0 in range(iNVar0+1):
            for i1 in range(iNVar1+1):
                pVar = iLabel1+str(i1)+iLabel0+str(i0);     
                if self._freeze :
                
                    start = value [i0*3+i1]
                    pXMin = value [i0*3+i1]-error[i0*3+i1]
                    pXMax = value [i0*3+i1]+error[i0*3+i1]
                    
                else: 
                    start = 0.0
                    pXMin = iXMin0
                    pXMax = iXMax0
                
                pRooVar = r.RooRealVar(pVar,pVar,0.0,pXMin,pXMax)
                print("========  here i0 %s i1 %s"%(i0,i1))
                print pVar
                print(" is : %s  +/- %s"%(value[i0*3+i1],error[i0*3+i1]))        
                iVars.append(pRooVar)
                self._all_vars.append(pRooVar)

    def GetWorkspaceInputs(self, pass_histograms, fail_histograms,iBin):

        roocategories = r.RooCategory("sample","sample") 
        roocategories.defineType("pass",1) 
        roocategories.defineType("fail",0) 
        data_rdh_pass = r.RooDataHist("data_obs_pass_"+iBin,"data_obs_pass_"+iBin,r.RooArgList(self._lMSD),pass_histograms["data_obs"])
        data_rdh_fail = r.RooDataHist("data_obs_fail_"+iBin,"data_obs_fail_"+iBin,r.RooArgList(self._lMSD), fail_histograms["data_obs"])
        data_rdh_comb  = r.RooDataHist("comb_data_obs","comb_data_obs",r.RooArgList(self._lMSD),r.RooFit.Index(roocategories),r.RooFit.Import("pass",data_rdh_pass),r.RooFit.Import("fail",data_rdh_fail)) 

        roofit_shapes = {}
        for sample in ["wqq", "zqq", "qcd", "tqq"]:
            roofit_shapes[sample] = self.GetRoofitHistObjects(pass_histograms[sample], fail_histograms[sample], sample, iBin)

        total_pdf_pass = r.RooAddPdf("tot_pass"+iBin,"tot_pass"+iBin,r.RooArgList(roofit_shapes["qcd"]["pass_epdf"]))
        total_pdf_fail = r.RooAddPdf("tot_fail"+iBin,"tot_fail"+iBin,r.RooArgList(roofit_shapes["qcd"]["fail_epdf"]))
        ewk_pdf_pass = r.RooAddPdf("ewk_pass"+iBin,"ewk_pass"+iBin,r.RooArgList(roofit_shapes["wqq"]["pass_epdf"],roofit_shapes["zqq"]["pass_epdf"], roofit_shapes["tqq"]["pass_epdf"]))
        ewk_pdf_fail = r.RooAddPdf("ewk_fail"+iBin,"ewk_fail"+iBin,r.RooArgList(roofit_shapes["wqq"]["fail_epdf"],roofit_shapes["zqq"]["fail_epdf"], roofit_shapes["tqq"]["fail_epdf"]))

        total_simulpdf  = r.RooSimultaneous("tot","tot",roocategories) 
        total_simulpdf.addPdf(total_pdf_pass,"pass") 
        total_simulpdf.addPdf(total_pdf_fail,"fail")     
        self._all_data.extend([data_rdh_pass,data_rdh_fail])
        self._all_shapes.extend([total_pdf_pass,total_pdf_fail,ewk_pdf_pass,ewk_pdf_fail])

        ## find out which to make global
        ## RooDataHist (data), then RooHistPdf of each electroweak
        # Previous return values 2 and 3 (RooAbsPdf (qcd,ewk)) removed by David on 19/1/2017, because they don't seem to be used. 
        return [
            data_rdh_pass, 
            data_rdh_fail, 
            #{"qcd":total_pdf_pass, "ewk":ewk_pdf_pass},
            #{"qcd":total_pdf_fail, "ewk":ewk_pdf_fail},
            {"wqq":roofit_shapes["wqq"]["pass_rdh"], "zqq":roofit_shapes["zqq"]["pass_rdh"], "tqq":roofit_shapes["tqq"]["pass_rdh"]}, 
            {"wqq":roofit_shapes["wqq"]["fail_rdh"], "zqq":roofit_shapes["zqq"]["fail_rdh"], "tqq":roofit_shapes["tqq"]["fail_rdh"]}, 
        ]

    # Get (RooHistPdf, RooExtendPdf, RooDataHist) for a pair of pass/fail histograms
    # - The RooExtendPdfs are coupled via their normalizations, N*eff or N*(1-eff). 
    def GetRoofitHistObjects(self, hist_pass, hist_fail, label="w", category="_cat0"):
        # normalization
        total_norm   = r.RooRealVar(label+"norm"+category, label+"norm"+category, (hist_pass.Integral()+hist_fail.Integral()), 0., 5.*(hist_pass.Integral()+hist_fail.Integral()))
        pass_norm  = r.RooFormulaVar(label+"fpass"+category, label+"norm"+category+"*(veff)", r.RooArgList(total_norm,self._lEff))
        fail_norm  = r.RooFormulaVar(label+"fqail"+category, label+"norm"+category+"*(1-veff)", r.RooArgList(total_norm,self._lEff))

        # shapes
        pass_rdh  = r.RooDataHist(label+"_pass_"+category, label+"_pass_"+category, r.RooArgList(self._lMSD),hist_pass)
        fail_rdh  = r.RooDataHist(label+"_fail_"+category, label+"_fail_"+category, r.RooArgList(self._lMSD),hist_fail) 
        pass_rhp      = r.RooHistPdf (label+"passh"+category, label+"passh"+category, r.RooArgList(self._lShift), r.RooArgList(self._lMSD), pass_rdh, 0)
        fail_rhp      = r.RooHistPdf (label+"failh"+category, label+"failh"+category, r.RooArgList(self._lShift), r.RooArgList(self._lMSD), fail_rdh, 0)

        # extended likelihood from normalization and shape above
        pass_epdf     = r.RooExtendPdf(label+"_passe_" +category, label+"pe" +category, pass_rhp, pass_norm)
        fail_epdf     = r.RooExtendPdf(label+"_faile_" +category, label+"fe" +category, fail_rhp, fail_norm)

        #lHist   = [pass_pdf,fail_rhp,pass_epdf,fail_epdf,pass_rdh,fail_rdh]
        return_dict = {
            "pass_rdh":pass_rdh, 
            "fail_rdh":fail_rdh,
            "pass_pdf":pass_rhp,
            "fail_pdf":fail_rhp, 
            "pass_epdf":pass_epdf,
            "fail_epdf":fail_epdf
        }
        self._all_vars.extend([total_norm, pass_norm, fail_norm])
        self._all_shapes.extend(return_dict.values())
        return return_dict

    def GetSignalInputs(self,iHP,iHF,iBin):
        # get signals
        lPSigs  = {}
        lFSigs  = {}
        for signal_name in self._signal_names:
            roofit_shapes = self.GetRoofitHistObjects(iHP[signal_name], iHF[signal_name],signal_name,iBin)
            lPSigs[signal_name] = roofit_shapes["pass_rdh"]
            lFSigs[signal_name] = roofit_shapes["fail_rdh"]
        return (lPSigs,lFSigs)


    #def MakeWorkspace(self,iOutput,iDatas,iFuncs,iVars,iCat="cat0",iShift=True):
    def MakeWorkspace(self, output_path, import_objects, category="cat0", shift=True):
        print "Making workspace " + "w_" + str(category)
        workspace = r.RooWorkspace("w_"+str(category))
        for import_object in import_objects:
            print "Importing {}".format(import_object.GetName())
            getattr(workspace,'import')(import_object, r.RooFit.RecycleConflictNodes())

        #for name, pFunc in iFuncs.iteritems():
        #    print "Importing {}".format(name)
        #    pFunc.Print()
        #    getattr(workspace,'import')(pFunc, r.RooFit.RecycleConflictNodes())
        #    # if iShift and pFunc.GetName().find("qq") > -1:
        #    #   (pFUp, pFDown)  = tools.shift(iVars[0],pFunc,5.)
        #    #   (pSFUp,pSFDown) = tools.smear(iVars[0],pFunc,0.05)
        #    #   getattr(workspace,'import')(pFUp,  r.RooFit.RecycleConflictNodes())
        #    #   getattr(workspace,'import')(pFDown,r.RooFit.RecycleConflictNodes())
        #    #   getattr(workspace,'import')(pSFUp,  r.RooFit.RecycleConflictNodes())
        #    #   getattr(workspace,'import')(pSFDown,r.RooFit.RecycleConflictNodes())

        #for pData in iDatas:
        #    pData.Print()
        #    getattr(workspace,'import')(pData, r.RooFit.RecycleConflictNodes()) # , r.RooFit.Rename("data_obs")

        if category.find("pass_cat1") == -1:
            workspace.writeToFile(output_path,False)
        else:
            workspace.writeToFile(output_path) 
        workspace.Print()
        # workspace.writeToFile(output_path)   

##############################################################################
##############################################################################
#### E N D   O F   C L A S S
##############################################################################
##############################################################################

##-------------------------------------------------------------------------------------
def LoadHistograms(f, pseudo, blind, useQCD, mass_range, blind_range):
    pass_hists = {}
    fail_hists = {}
    f.ls()

    pass_hists_bkg = {}
    fail_hists_bkg = {}
    background_names = ["wqq", "zqq", "qcd", "tqq"]
    for i, bkg in enumerate(background_names):
        if bkg=='qcd':
            qcd_fail = f.Get('qcd_fail')
            if useQCD:
                qcd_pass = f.Get('qcd_pass')
            else:
                qcd_pass_real = f.Get('qcd_pass').Clone('qcd_pass_real')
                qcd_pass = qcd_fail.Clone('qcd_pass')
                qcd_pass_real_integral = 0
                qcd_fail_integral = 0
                for i in range(1,qcd_pass_real.GetNbinsX()+1):
                    for j in range(1,qcd_pass_real.GetNbinsY()+1):
                        if qcd_pass_real.GetXaxis().GetBinCenter(i) > mass_range[0] and qcd_pass_real.GetXaxis().GetBinCenter(i) < mass_range[1]:
                            qcd_pass_real_integral += qcd_pass_real.GetBinContent(i,j)
                            qcd_fail_integral += qcd_fail.GetBinContent(i,j)                   
                qcd_pass.Scale(qcd_pass_real_integral/qcd_fail_integral) # qcd_pass = qcd_fail * eff(pass)/eff(fail)
            pass_hists_bkg["qcd"] = qcd_pass
            fail_hists_bkg["qcd"] = qcd_fail
            print 'qcd pass integral', qcd_pass.Integral()
            print 'qcd fail integral', qcd_fail.Integral()
        else:            
            pass_hists_bkg[bkg] = f.Get(bkg+'_pass')
            fail_hists_bkg[bkg] = f.Get(bkg+'_fail')
        
    if pseudo:        
        for i, bkg in enumerate(background_names):
            if i==0:
                pass_hists["data_obs"] = pass_hists_bkg[bkg].Clone('data_obs_pass')
                fail_hists["data_obs"] = fail_hists_bkg[bkg].Clone('data_obs_fail')
            else:                
                pass_hists["data_obs"].Add(pass_hists_bkg[bkg])
                fail_hists["data_obs"].Add(fail_hists_bkg[bkg])        
    else:
        pass_hists["data_obs"] = f.Get('data_obs_pass')
        fail_hists["data_obs"] = f.Get('data_obs_fail')

    #signals
    pass_hists_sig = {}
    fail_hists_sig = {}
    masses=[50,75,125,100,150,250,300]
    signal_names = []
    for mass in masses:
        print "Signal histogram name = Pbb_" +str(mass)+"_pass"
        passhist = f.Get("Pbb" + "_" +str(mass)+"_pass").Clone()
        failhist = f.Get("Pbb" + "_" +str(mass)+"_fail").Clone()
        for hist in [passhist, failhist]:
            for i in range(0,hist.GetNbinsX()+2):
                for j in range(0,hist.GetNbinsY()+2):
                    if hist.GetBinContent(i,j) <= 0:
                        hist.SetBinContent(i,j,0)
        pass_hists_sig["Pbb_" + str(mass)] = passhist
        fail_hists_sig["Pbb_" + str(mass)] = failhist
        signal_names.append("Pbb_" + str(mass))
        #pass_hists_sig.append(f.Get(sig+str(mass)+"_pass"))

    pass_hists.update(pass_hists_bkg)
    pass_hists.update(pass_hists_sig)
    fail_hists.update(fail_hists_bkg)
    fail_hists.update(fail_hists_sig)
    
    for histogram in (pass_hists.values() + fail_hists.values()):
        if blind:            
            for i in range(1,histogram.GetNbinsX()+1):
                for j in range(1,histogram.GetNbinsY()+1):
                    if histogram.GetXaxis().GetBinCenter(i) > blind_range[0] and histogram.GetXaxis().GetBinCenter(i) < blind_range[1]:
                        print "blinding signal region for %s, mass bin [%i,%i] "%(histogram.GetName(),histogram.GetXaxis().GetBinLowEdge(i),histogram.GetXaxis().GetBinUpEdge(i))
                        histogram.SetBinContent(i,j,0.)
                        print histogram.GetBinContent(i,j)
        histogram.SetDirectory(0)  

    # print "lengths = ", len(pass_hists), len(fail_hists)
    # print pass_hists;
    # print fail_hists;
    return (pass_hists,fail_hists)
    # return (fail_hists,pass_hists)

##-------------------------------------------------------------------------------------

def main(options,args):
    
    ifile = options.ifile
    odir = options.odir

    # Load the input histograms
    #   - 2D histograms of pass and fail mass,pT distributions
    #   - for each MC sample and the data
    f = r.TFile.Open(ifile)
    (pass_hists,fail_hists) = LoadHistograms(f,options.pseudo,options.blind,options.useQCD);

    # Build the workspacees
    RhalphabetBuilder(pass_hists, fail_hists, odir, mass_fit=options.mass_fit, freeze_poly=args.freeze_poly)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option('-i','--ifile', dest='ifile', default = 'hist_1DZbb.root',help='file with histogram inputs', metavar='ifile')
    parser.add_option('-o','--odir', dest='odir', default = './',help='directory to write plots', metavar='odir')
    parser.add_option('--pseudo', action='store_true', dest='pseudo', default =False,help='use MC', metavar='pseudo')
    parser.add_option('--blind', action='store_true', dest='blind', default =False,help='blind signal region', metavar='blind')
    parser.add_option('--use-qcd', action='store_true', dest='useQCD', default =False,help='use real QCD MC', metavar='useQCD')
    parser.add_option('--mass_fit', action='store_true', dest='mass_fit', default =False,help='Do mass fit instead of rho fit', metavar='mass_fit')
    parser.add_option('--freeze_poly', action='store_true', dest='freeze_poly', default =False,help='freeze_poly pol values', metavar='freeze_poly')
    
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

