import os
import sys
from array import array
import re
import ROOT
from ROOT import *
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/slc6_amd64_gcc491/libDAZSLEPhiBBPlusJet.so"))
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so"))
import DAZSLE.PhiBBPlusJet.analysis_configuration as config
import DAZSLE.PhiBBPlusJet.style as style
seaborn = Root.SeabornInterface()
seaborn.Initialize()

def make_plots(base_file_path, combine_file_path):
	# Get data from base
	base_file = TFile(base_file_path, "READ")
	