import os
import sys
from glob import glob

step1 = False
step2 = True

tau21_values = [0.4, 0.45, 0.5, 0.525, 0.55, 0.575, 0.6, 0.65, 0.7]
for tau_21 in tau21_values:
	wp_string = "tau21ddt_{}".format(tau_21)
	print wp_string
	directory = "~/DAZSLE/data/LimitSetting/Xbb_inputs/{}".format(wp_string)

	if step1:
		os.system("mkdir -pv {}/cards_mcstat/".format(directory))
		print "[makeCardsXbb]"
		os.system("python makeCardsXbb.py -o {}/cards_mcstat/ -i {}/histograms_SR_{}_AK8.root".format(directory, directory, wp_string))
		print "[buildRhalphabetXbb]"
		os.system("python buildRhalphabetXbb.py -i {}/histograms_SR_{}_AK8.root -o {}/cards_mcstat/ --pseudo".format(directory, wp_string, directory))
		#print "[writeMuonCRDatacardXbb]"
		#os.system("python writeMuonCRDatacardXbb.py -i {} -o {}/cards_mcstat/".format(directory, directory))
	if step2:
		cwd = os.getcwd()
		print "[debug] Glob pattern = " + directory + "/cards_mcstat/DMSbb*"
		signal_dirs = glob(directory + "/cards_mcstat/DMSbb*")
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


