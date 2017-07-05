import ROOT as r
import sys
import math
import array
import os

def uncorrelate(hists, sysName, suppressLevel=None):
  """Replaces each histogram whose name contains 'sysName' with many copies that represent uncorrelated bin-by-bin systematics.
  suppressLevel: if provided, new histograms will only be created for bins that differ from nominal by a fractional amount greater than suppressLevel."""
  #get all histograms that match the input string
  toUncorrelate = [name for name in hists if sysName in name]
  print "Treating the following distributions as uncorrelated for",sysName,": "
  for name in toUncorrelate: print name
  
  #get names of individual systematics
  systNames = []
  for name in toUncorrelate:
    systName = name.replace("Up","").replace("Down","")
    if systName not in systNames:
      systNames.append(systName)

  for name in systNames:
    print("Uncorrelating "+name)
    #get histogram with central values
    centerName = name.split("_")[:-1]
    centerName = '_'.join(centerName)
    #get up and down variants
    upName = name+'Up'
    downName = name+'Down'
    uncName = name.split("_")[-1]
    print("Central values taken from "+centerName)
    #for each bin create a new histogram in which that bin is up/down and the rest are centered
    for b in range(1,hists[centerName].GetNbinsX()+1):
      newHistUpName = centerName+"_"+uncName+str(b)+"Up"
      newHistDownName = centerName+"_"+uncName+str(b)+"Down"

      #check level of agreement with the nominal
      if suppressLevel is not None:
        #get percent difference from nominal
        if hists[centerName].GetBinContent(b) > 0:
          percDifferenceUp = abs(hists[upName].GetBinContent(b)-hists[centerName].GetBinContent(b))/hists[centerName].GetBinContent(b)
          percDifferenceDown = abs(hists[downName].GetBinContent(b)-hists[centerName].GetBinContent(b))/hists[centerName].GetBinContent(b)
          percDifference = max(percDifferenceUp, percDifferenceDown)
          if percDifference <= suppressLevel: 
            #print "Suppressing nuisance in bin",b,"(agrees at",percDifference,"level)"
            continue
        elif hists[upName].GetBinContent(b) == hists[centerName].GetBinContent(b) and hists[downName].GetBinContent(b) == hists[centerName].GetBinContent(b): 
            #print "Suppressing nuisance in bin",b,"because there is no change from the nominal"
            continue

      #new up histogram
      hists[newHistUpName] = hists[centerName].Clone(newHistUpName)
      hists[newHistUpName].SetDirectory(0)
      hists[newHistUpName].SetBinContent(b, hists[upName].GetBinContent(b)) #new hist has the unperturbed value in every bin except one
      hists[newHistUpName].SetBinError(b, hists[upName].GetBinError(b))

      #new down histogram
      hists[newHistDownName] = hists[centerName].Clone(newHistDownName)
      hists[newHistDownName].SetDirectory(0)
      hists[newHistDownName].SetBinContent(b, hists[downName].GetBinContent(b)) #new hist has the unperturbed value in every bin except one
      hists[newHistDownName].SetBinError(b, hists[downName].GetBinError(b))

    #remove the original histogram
    del hists[upName]
    del hists[downName]

def shift(iVar,iDataHist,iShift=0.):
    '''
      Creates up/down shifted RooDataHists
      iDataHist=input RooDataHist
      iShift=shift value
      Shifting is performed by creating a new x variable, x'=x-+xshift.
    '''
    lInt    = iDataHist.createHistogram("x").Integral()
    lDM     = r.RooRealVar   ("Xdm","Xdm", 0.,-10,10)
    lShift  = r.RooFormulaVar("Xshift",iVar.GetName()+"-Xdm",r.RooArgList(iVar,lDM))  
    lSPdf   = r.RooHistPdf(iDataHist.GetName()+"P",iDataHist.GetName()+"P", r.RooArgList(lShift),r.RooArgList(iVar),iDataHist,0)
    lDM.setVal(iShift)
    lHUp   = lSPdf.createHistogram("x")
    lHUp.Scale(lInt)
    lUp    = r.RooDataHist(iDataHist.GetName()+"_scaleUp",iDataHist.GetName()+"_scaleUp", r.RooArgList(iVar),lHUp)
    lDM.setVal(-iShift)
    lHDown = lSPdf.createHistogram("x")
    lHDown.Scale(lInt)
    lDown  = r.RooDataHist(iDataHist.GetName()+"_scaleDown",iDataHist.GetName()+"_scaleDown", r.RooArgList(iVar),lHDown)
    return (lUp,lDown)

def smear(iVar,iDataHist,iScale=0.):
    '''
      Creates up/down smeared RooDataHists
      iDataHist=input RooDataHist
      iScale=smear value
      Smearing is performed by creating a new x variable, x'=(x-xmean)/(1+-scale)+xmean.
    '''
    lDM     = r.RooRealVar("Xshift","Xshift", 1.,0.,2.)
    lVar    = iDataHist.createHistogram("x").GetMean()
    lInt    = iDataHist.createHistogram("x").Integral()
    lShift  = r.RooFormulaVar("Xsmear","("+iVar.GetName()+"-"+str(lVar)+")/Xshift+"+str(lVar),r.RooArgList(iVar,lDM))  
    lHPdf   = r.RooHistPdf(iDataHist.GetName()+"S",iDataHist.GetName()+"S", r.RooArgList(lShift),r.RooArgList(iVar),iDataHist,0)
    lDM.setVal(1.+iScale)
    lHUp = lHPdf.createHistogram("x")
    lHUp.Scale(lInt)
    lUp = r.RooDataHist(iDataHist.GetName()+"_smearUp",iDataHist.GetName()+"_smearUp", r.RooArgList(iVar),lHUp)    
    lDM.setVal(1.-iScale)
    lHDown = lHPdf.createHistogram("x")
    lHDown.Scale(lInt)
    lDown  = r.RooDataHist(iDataHist.GetName()+"_smearDown",iDataHist.GetName()+"_smearDown", r.RooArgList(iVar),lHDown)
    return [lUp,lDown]    
