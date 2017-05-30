import os
import sys
from glob import glob

step1 = False
step2 = True

tau21_values = [0.6, 0.65, 0.7] # 0.4, 0.45, 0.5, 0.525, 0.55, 0.575, 
dcsv_values = [0.7, 0.75, 0.8, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975]
#tau21_values = [0.55]
#dcsv_values = [0.9]
jet_types = ["AK8", "CA15"]
for tau_21 in tau21_values:
	for dcsv in dcsv_values:
		for jet_type in jet_types:
			print "On {}/{}/{}".format(tau_21, dcsv, jet_type)
			wp_string = "tau21ddt{}_dcsv{}_{}".format(tau_21, dcsv, jet_type)
			print wp_string
			directory = "/uscms_data/d3/dryu/DAZSLE/data/LimitSetting/Xbb_inputs/{}".format(wp_string)

			if step1:
				os.system("mkdir -pv {}/cards_mcstat/".format(directory))
				print "[makeCardsXbb]"
				os.system("python makeCardsXbb.py -o {}/cards_mcstat/ -i {}/histograms_SR_{}.root".format(directory, directory, wp_string))
				print "[buildRhalphabetXbb]"
				previous_directory = os.getcwd()
				working_directory = directory + "/cards_mcstat/"
				os.chdir(working_directory)
				run_script_path = working_directory + "/run_rhalphabet.sh"
				run_script = open(run_script_path, "w")
				run_script.write("#!/bin/bash\n")
				run_script.write("python $CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/fitting/PbbJet/buildRhalphabetXbb.py -i histograms_SR_{}.root -o . --pseudo --useQCD\n".format(wp_string))
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
				print "[debug] Glob pattern = " + directory + "/cards_mcstat/*Sbb*"
				signal_dirs = glob(directory + "/cards_mcstat/*Sbb*")
				print "List of signal directories:"
				print signal_dirs
				for signal_dir in signal_dirs:
					print "Working directory: " + signal_dir
					os.chdir(signal_dir)
					os.system("cp ../*base.root .")
					os.system("combineCards.py card_rhalphabet_cat1.txt card_rhalphabet_cat2.txt card_rhalphabet_cat3.txt card_rhalphabet_cat4.txt  card_rhalphabet_cat5.txt card_rhalphabet_cat6.txt   >  card_rhalphabet.txt")
					# Make and submit condor job
					run_script = open("{}/run_combine.sh".format(signal_dir), "w")
					run_script.write("#!/bin/bash\n")
					run_script.write("combine -M Asymptotic -v 2 -t -1 card_rhalphabet.txt 2>&1 | tee combine_log.txt\n")
					run_script.close()
					csub_script = open("{}/submit.sh".format(signal_dir), "w")
					csub_script.write("#!/bin/bash\n")
					csub_script.write("csub run_combine.sh --cmssw --no_retar -F card_rhalphabet.txt,base.root,rhalphabase.root\n")
					csub_script.close()
					os.system("source {}/submit.sh".format(signal_dir))
					os.chdir(cwd)


print "Done."