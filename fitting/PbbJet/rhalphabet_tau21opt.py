import os
import sys
from glob import glob
from joblib import Parallel, delayed

from ROOT import *
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/slc6_amd64_gcc491/libDAZSLEPhiBBPlusJet.so"))
import DAZSLE.PhiBBPlusJet.analysis_configuration as config

step0 = False # makeCardsXbb
step1 = False # buildRhalphabetXbb
step2 = True # combineCards.py and call combine

tau21_values = [0.4, 0.45, 0.5, 0.525, 0.55, 0.575, 0.6, 0.65, 0.7] #   
dcsv_values = [0.7, 0.75, 0.8, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975]
#tau21_values = [0.55]
#dcsv_values = [0.9]
jet_types = ["AK8", "CA15"] # CA15

def run_single(tau_21, dcsv, jet_type):
	print "On {}/{}/{}".format(tau_21, dcsv, jet_type)
	wp_string = "tau21ddt{}_dcsv{}_{}".format(tau_21, dcsv, jet_type)
	print wp_string
	directory = "/uscms_data/d3/dryu/DAZSLE/data/LimitSetting/Xbb_inputs/{}".format(wp_string)

	if step0:
		os.system("mkdir -pv {}/cards_mcstat/".format(directory))
		print "[makeCardsXbb]"
		os.system("python makeCardsXbb.py -o {}/cards_mcstat/ -i {}/histograms_SR_{}.root --jet_type {}".format(directory, directory, wp_string, jet_type))
	elif step1:
		os.system("mkdir -pv {}/cards_mcstat/".format(directory))
		print "[buildRhalphabetXbb]"
		previous_directory = os.getcwd()
		working_directory = directory + "/cards_mcstat/"
		os.chdir(working_directory)
		run_script_path = working_directory + "/run_rhalphabet.sh"
		run_script = open(run_script_path, "w")
		run_script.write("#!/bin/bash\n")
		run_script.write("python $CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/fitting/PbbJet/buildRhalphabetXbb.py -i histograms_SR_{}.root -o . --pseudo --jet_type {}\n".format(wp_string, jet_type))
		run_script.close()
		csub_script = open(working_directory + "/csub.sh", "w")
		csub_script.write("#!/bin/bash\n")
		files_to_transfer = ["{}/histograms_SR_{}.root".format(directory, wp_string)]
		csub_script.write("csub {} --cmssw --no_retar -F {}\n".format(run_script_path, ",".join(files_to_transfer)))
		csub_script.close()
		os.system("source {}/csub.sh".format(working_directory))
		os.chdir(previous_directory)
		#os.system("python buildRhalphabetXbb.py -i {}/histograms_SR_{}.root -o {}/cards_mcstat/ --pseudo --useQCD".format(directory, wp_string, directory))
		#print "[writeMuonCRDatacardXbb]"
		#os.system("python writeMuonCRDatacardXbb.py -i {} -o {}/cards_mcstat/".format(directory, directory))
	elif step2:
		cwd = os.getcwd()
		#print "[debug] Glob pattern = " + directory + "/cards_mcstat/*Sbb*"
		#signal_dirs = glob(directory + "/cards_mcstat/*Sbb*")
		#print "List of signal directories:"
		#print signal_dirs
		#for signal_dir in signal_dirs:
		working_directory = directory + "/cards_mcstat/"
		os.chdir(working_directory)
		run_script = open("{}/run_combine.sh".format(working_directory), "w")
		run_script.write("#!/bin/bash\n")
		datacards = []
		signal_masses = []
		for signal_name in config.simulated_signal_names:
			signal_dir = directory + "/cards_mcstat/" + signal_name
			#print "Working directory: " + signal_dir
			os.chdir(signal_dir)
			#os.system("cp ../*base.root .")
			os.system("combineCards.py card_rhalphabet_cat1.txt card_rhalphabet_cat2.txt card_rhalphabet_cat3.txt card_rhalphabet_cat4.txt  card_rhalphabet_cat5.txt card_rhalphabet_cat6.txt   >  card_rhalphabet_{}.txt".format(signal_name))
			datacards.append("{}/card_rhalphabet_{}.txt".format(signal_dir, signal_name))
			signal_masses.append(str(config.signal_masses[signal_name]))
			os.chdir(working_directory)
		run_script.write("datacards=( " + " ".join([os.path.basename(x) for x in datacards]) + " )\n")
		run_script.write("signal_names=( " + " ".join(config.simulated_signal_names) + " )\n")
		run_script.write("signal_masses=( " + " ".join(signal_masses) + " )\n")
		run_script.write("combine -M Asymptotic -v 0 -n ${signal_names[$1]} -t -1 ${datacards[$1]} 2>&1 | tee combine_log_${signal_names[$1]}.txt\n") # -m ${signal_masses[$1]}
		run_script.close()
		csub_script = open("{}/submit.sh".format(working_directory), "w")
		csub_script.write("#!/bin/bash\n")
		csub_script.write("csub run_combine.sh --cmssw --no_retar -F {},base.root,rhalphabase.root --queue_n {}\n".format(",".join(datacards), len(datacards)))
		csub_script.close()
		os.system("source {}/submit.sh".format(working_directory))
		os.chdir(cwd)

Parallel(n_jobs=8)(delayed(run_single)(tau_21, dcsv, jet_type) for tau_21 in tau21_values for dcsv in dcsv_values for jet_type in jet_types)

print "Done."