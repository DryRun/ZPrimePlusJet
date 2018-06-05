import os
import sys
from ROOT import *
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libDAZSLEPhiBBPlusJet.so"))
import DAZSLE.ZPrimePlusJet.xbb_config as config

def csub_combine(signal_name, jet_type, method="Asymptotic", region="SR", decidata=False, pseudodata=False, card_name=None, verbose=0, nrho=2, npt=1, prefit=False, freeze_nuisances=None, freeze_groups=None):
	if not card_name:
		# Infer card name automatically
		if region == "SR":
			card_name = "card_rhalphabet_muonCR.txt" # This is the SR+muCR datacard
		elif region == "N2SR":
			card_name = "card_rhalphabet_nomuonCR.txt" # This is the SR+muCR datacard
		elif region == "N2CR":
			card_name = "card_rhalphabet_N2CR.txt" # This is the N2CR-only datacard
		elif region == "muCR":
			card_name = "datacard_muonCR.txt" # This is the muCR-only datacard

	datacard_directory = config.get_datacard_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region)
	working_directory = config.get_limit_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region)
	os.system("mkdir -pv {}".format(working_directory))
	cwd = os.getcwd()
	os.chdir(working_directory)
	input_files = ["{}/base.root".format(datacard_directory), "{}/rhalphabase.root".format(datacard_directory), "{}/{}".format(datacard_directory, card_name)]
	if region == "SR" or region == "muCR":
		input_files.append("{}/workspace_muonCR.root".format(datacard_directory))

	# Run script
	job_name = "combine_{}_{}_{}_{}_r{}p{}".format(signal_name, jet_type, method, region, nrho, npt)
	run_script_path = "{}/run_combine_{}.sh".format(working_directory, job_name)
	run_script = open(run_script_path, 'w')
	run_script.write("#!/bin/bash\n")
	run_script.write("cp {} tmpcard.txt\n".format(card_name))
	run_script.write("mkdir -pv plots_{}\n".format(job_name))
	run_script.write("sed -i 's@{}/@@g' tmpcard.txt\n".format(datacard_directory)) # Combine cards sometimes pick up the absolute path to the workspace, which doesn't work on condor

	# Prefit
	if prefit:
		prefit_fix_pars_str = ""
		prefit_zero_pars_str = ""
		for irho in xrange(config.analysis_parameters[jet_type]["MAX_NRHO"]+1):
			for ipt in xrange(config.analysis_parameters[jet_type]["MAX_NPT"]+1):
				if irho > nrho or ipt > npt:
					prefit_fix_pars_str += "r{}p{}:0,".format(irho, ipt)
				else:
					prefit_zero_pars_str += "r{}p{},".format(irho, ipt)
		prefit_fix_pars_str = prefit_fix_pars_str[:-1] # Chop trailing comma
		prefit_zero_pars_str = prefit_zero_pars_str[:-1] # Chop trailing comma
		cats = config.analysis_parameters[jet_type]["FIT_PT_BINS"]
		#if jet_type == "CA15" and region == "N2SR":
		#	cats = [1,2,3,4,5,6]
		run_script.write("python $CMSSW_BASE/python/DAZSLE/ZPrimePlusJet/prefit_workspace.py --base_path base.root --rhalphabase_path rhalphabase.root --fix_pars_rhalphabet {} --zero_pars_rhalphabet {}  --signals {} --cats {} --no_backup_original\n".format(prefit_fix_pars_str, prefit_zero_pars_str, signal_name, ",".join([str(x) for x in cats])))

	combine_command = "combine -M {} -v {} tmpcard.txt -n {}".format(method, verbose, job_name) # -t -1 --toysFreq
	if method == "Asymptotic":
		combine_command += " --saveWorkspace"
	elif method == "MaxLikelihoodFit":
		combine_command += " --saveNormalizations --plot --saveShapes --saveWorkspace --out plots_{}".format(job_name)
	frozen_nuisance_list = []		
	for irho in xrange(config.analysis_parameters[jet_type]["MAX_NRHO"]+1):
		for ipt in xrange(config.analysis_parameters[jet_type]["MAX_NPT"]+1):
			if irho > nrho or ipt > npt:
				frozen_nuisance_list.append("r{}p{}".format(irho, ipt))
	if freeze_nuisances:
		frozen_nuisance_list.extend(freeze_nuisances)
	freeze_string = " --freezeNuisances " + ",".join(frozen_nuisance_list)
	combine_command += freeze_string

	if freeze_groups:
		combine_command += " --freezeNuisanceGroups " + ",".join(freeze_groups)
	combine_command += " 2>&1 | tee log_{}.txt".format(job_name)
	run_script.write(combine_command + "\n")

	# Prepare return objects
	run_script.write("tar -czvf plots_{}.tar.gz plots_{}\n".format(job_name, job_name))
	#run_script.write("mv base.root base_{}.root\n".format(job_name))
	run_script.write("rm -f base.root\n")
	run_script.write("rm -f rhalphabase.root\n")
	run_script.write("mv {} {}.{}\n".format(card_name, card_name, job_name))
	run_script.close()
	print combine_command

	# csub script
	csub_script_path = "{}/csub_{}.sh".format(working_directory, job_name)
	print csub_script_path
	csub_script = open(csub_script_path, 'w')
	csub_script.write("#!/bin/bash\n")
	csub_command = "csub {} -F {} --cmssw --no_retar ".format(run_script_path, ",".join(input_files))
	csub_script.write(csub_command + "\n")
	csub_script.close()

	# run
	os.system("source {}".format(csub_script_path))

	os.chdir(cwd)

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Run combine on condor")
	parser.add_argument("--regions", type=str, default="SR", help="SR, muCR, or N2CR")
	parser.add_argument("--methods", type=str, default="Asymptotic", help="Asymptotic or MaxLikelihoodFit")
	parser.add_argument("--jet_types", type=str, default="AK8,CA15")
	parser.add_argument("--decidata", action="store_true", help="Use 1/10 data")
	parser.add_argument("--pseudodata", action="store_true", help="Use MC pseudodata")
	parser.add_argument("--card", type=str, help="Manually specify card name")
	parser.add_argument("--signals", type=str, default="all", help="all or comma-separated list of signal names")
	parser.add_argument("--verbose", type=int, default="0", help="Combine verbosity")
	parser.add_argument("--nrho", type=int, help="Polynomial degree in rho")
	parser.add_argument("--npt", type=int, help="Polynomial degree in pt")
	parser.add_argument("--freeze_nuisances", type=str, help="Freeze nuisance parameters")
	parser.add_argument("--freeze_groups", type=str, help="Freeze nuisance parameters")
	#parser.add_argument("--prefit", action="store_true", help="Prefit the workspace") # Prefit should always be true, otherwise you get nonzero poly parameters!
	args = parser.parse_args()

	if args.signals == "all":
		signal_names = {
			"AK8":["Sbb{}".format(x) for x in range(50,425,25)] + ["PSbb{}".format(x) for x in range(50,425,25)],
			"CA15":["Sbb{}".format(x) for x in range(50,525,25)] + ["PSbb{}".format(x) for x in range(50,525,25)],
		}
	else:
		signal_names = {
			"AK8":args.signals.split(","),
			"CA15":args.signals.split(","),
		}

	if args.freeze_nuisances:
		freeze_nuisances = args.freeze_nuisances.split(",")
	else:
		freeze_nuisances = None
	if args.freeze_groups:
		freeze_groups = args.freeze_groups.split(",")
	else:
		freeze_groups = None
	for region in args.regions.split(","):
		for method in args.methods.split(","):
			for jet_type in args.jet_types.split(","):
				for signal_name in signal_names[jet_type]:
					if args.nrho:
						nrho = args.nrho
					else:
						nrho = config.analysis_parameters[jet_type]["DEFAULT_NRHO"]
					if args.npt:
						npt = args.npt
					else:
						npt = config.analysis_parameters[jet_type]["DEFAULT_NPT"]
					csub_combine(signal_name, jet_type, method=method, region=region, decidata=args.decidata, pseudodata=args.pseudodata, verbose=args.verbose, nrho=nrho, npt=npt, prefit=False, freeze_nuisances=freeze_nuisances, freeze_groups=freeze_groups)

