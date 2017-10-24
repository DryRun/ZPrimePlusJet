import os
import sys
from ROOT import *
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libDAZSLEPhiBBPlusJet.so"))
import DAZSLE.ZPrimePlusJet.xbb_config as config

ftest_directory = "/uscms/home/dryu/DAZSLE/data/Fits/ftest/"

# Submit fits via condor
def run_fits(signal_name, jet_type, region, nrhos=[1,2,3,4,5,6], npts=[1,2], pseudodata=False, decidata=False, ntoys=500, signal_mu=0., freeze_nuisances=None):
	datacard_directory = config.get_datacard_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region)
	if region == "SR":
		card_name = "card_rhalphabet_muonCR.txt" # This is the SR+muCR datacard
	elif region == "N2CR":
		card_name = "card_rhalphabet_N2CR.txt" # This is the N2CR-only datacard

	for nrho in nrhos:
		for npt in npts:
			# Freeze higher orders
			freeze_poly_string = ""
			for irho in xrange(config.analysis_parameters[jet_type]["MAX_NRHO"]+1):
				for ipt in xrange(config.analysis_parameters[jet_type]["MAX_NPT"]+1):
					if irho > nrho or ipt > npt:
						freeze_poly_string += "r{}p{},".format(irho, ipt)
			freeze_poly_string = freeze_poly_string[:-1]
			if freeze_nuisances:
				freeze_nuisances = freeze_nuisances + "," + freeze_poly_string
			else:
				freeze_nuisances = freeze_poly_string

			job_name = "ftestfits_{}_{}_{}_r{}p{}".format(signal_name, jet_type, region, nrho, npt)
			fit_directory = config.get_ftest_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region, nrho=nrho, npt=npt)
			os.system("mkdir -pv {}".format(fit_directory))
			cwd = os.getcwd()
			os.chdir(fit_directory)

			# Run script
			run_script_path = "{}/run_ftest_fits.sh".format(fit_directory)
			run_script = open(run_script_path, "w")
			run_script.write("#!/bin/bash\n")
			run_script.write("cp {} tmpcard.txt\n".format(card_name))
			run_script.write("sed -i 's@{}/@@g' tmpcard.txt\n".format(datacard_directory)) # Combine cards sometimes pick up the absolute path to the workspace, which doesn't work on condor

			command1 = "combine -M GoodnessOfFit {}  --rMax 20 --rMin -20 --fixedSignalStrength 0 --algorithm saturated -n {} ".format("tmpcard.txt", job_name + "_centralfit")
			if freeze_nuisances:
				command1 += " --freezeNuisances {}".format(freeze_nuisances)
			run_script.write(command1 + "\n")

			command2 = "combine -M GenerateOnly {} --rMax 20 --rMin -20 --toysFrequentist -t {} --expectSignal {} --saveToys -n {}".format("tmpcard.txt", ntoys, signal_mu, job_name)
			run_script.write(command2 + "\n")

			command3 = "combine -M GoodnessOfFit {} --rMax 20 --rMin -20 -t {} --fixedSignalStrength 0 --toysFile higgsCombine{}.GenerateOnly.mH120.123456.root --algorithm saturated -n {}".format("tmpcard.txt", ntoys, job_name, job_name + "_toysfit")
			if freeze_nuisances:
				command3 +=  " --freezeNuisances {}".format(freeze_nuisances)
			run_script.write(command3 + "\n")
			run_script.close()

			files_to_transfer = ["{}/{}".format(datacard_directory, card_name), "{}/base.root".format(datacard_directory), "{}/rhalphabase.root".format(datacard_directory)]
			csub_command = "csub {} --cmssw --no_retar -F {}".format(run_script_path, ",".join(files_to_transfer))
			csub_script_path = "{}/csub_submit.sh".format(fit_directory)
			csub_script = open(csub_script_path, "w")
			csub_script.write("#!/bin/bash\n")
			csub_script.write(csub_command + "\n")
			csub_script.close()

			os.system("source {}".format(csub_script_path))
			#print "cat {}".format(csub_script_path)

def get_fit_likelihoods(signal_name, jet_type, region, pseudodata=False, decidata=False, nrho=2, npt=1):
	fit_directory = config.get_ftest_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region, nrho=nrho, npt=npt)



if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="One stop shop for F tests")
	action_group = parser.add_mutually_exclusive_group() 
	action_group.add_argument('--fits', action="store_true", help="Run fits")
	action_group.add_argument('--pval', action="store_true", help="Compute F statistics and p-values")
	action_group.add_argument('--plots', action="store_true", help="Plot fit results")

	parser.add_argument("--nrhos", type=str, help="nrho values (comma-separated)")
	parser.add_argument("--npts", type=str, help="npt values (comma-separated)")
	parser.add_argument("--signals", type=str, default="Sbb125", help="Signal names to run")
	parser.add_argument("--jet_types", type=str, default="AK8,CA15", help="Jet types to run")
	parser.add_argument("--region", type=str, default="SR", help="SR or N2CR")
	parser.add_argument("--pseudodata", action="store_true", help="Run over QCD pseudodata rather than data")
	parser.add_argument("--decidata", action="store_true", help="Run over 1/10 data")
	parser.add_argument("--mu", type=float, default=0., help="Signal mu for toy generation")
	parser.add_argument("--ntoys", type=int, default=500, help="Number of toys")
	parser.add_argument("--freeze_nuisances", type=str, help="Pass to --freezeNuisances argument of combine")
	args = parser.parse_args()

	jet_types = args.jet_types.split(",")
	signals = args.signals.split(",")
	nrhos = [int(x) for x in args.nrhos.split(",")]
	npts = [int(x) for x in args.npts.split(",")]

	if args.fits:
		for signal in signals:
			for jet_type in jet_types:
				run_fits(signal, jet_type, args.region, nrhos=nrhos, npts=npts, pseudodata=args.pseudodata, decidata=args.decidata, ntoys=500, signal_mu=args.mu, freeze_nuisances=args.freeze_nuisances)
	print "Done."
