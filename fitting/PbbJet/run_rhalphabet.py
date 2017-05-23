import os
import sys
from glob import glob

step1 = True
step2 = True

analyses = ["qcd", "data"]
jet_types = ["AK8", "CA15"]
run_fit = True

for analysis in analyses:
	for jet_type in jet_types:
		directory = "/uscms_data/d3/dryu/DAZSLE/data/LimitSetting/Xbb_inputs/SR_{}".format(jet_type)
		if analysis == "qcd":
			card_dir = directory + "/cards_qcd_mcstat/"
		elif analysis == "data":
			card_dir = directory + "/cards_data_mcstat/"
		if step1:
			os.system("mkdir -pv {}".format(card_dir))
			print "[makeCardsXbb]"
			os.system("python makeCardsXbb.py -o {} -i {}/histograms_SR_{}.root".format(card_dir, directory, jet_type))
			print "[buildRhalphabetXbb]"
			os.system("python buildRhalphabetXbb.py -i {}/histograms_SR_{}.root -o {}/ --pseudo --useQCD".format(directory, jet_type, card_dir))
			#print "[writeMuonCRDatacardXbb]"
			#os.system("python writeMuonCRDatacardXbb.py -i {} -o {}/cards_mcstat/".format(directory, directory))
		if step2:
			cwd = os.getcwd()
			print "[debug] Glob pattern = " + card_dir + "/*Sbb*"
			signal_dirs = glob(card_dir + "/*Sbb*")
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

				if run_fit:
					# Run MaxLikelihoodFit for plots
					os.system("combine -M MaxLikelihoodFit -v2 --saveWithUncertainties --saveShapes --plots --out plots card_rhalphabet.txt 2>&1 | tee fit_log.txt")



