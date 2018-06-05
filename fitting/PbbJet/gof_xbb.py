import os
import sys
import math
import time
import numpy as np
from ROOT import *
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libDAZSLEPhiBBPlusJet.so"))
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
import DAZSLE.ZPrimePlusJet.xbb_config as config

gof_directory = "/uscms/home/dryu/DAZSLE/data/Fits/gof/"

def run_gof(signal_name, jet_type, algorithms=["saturated", "KS", "AD"], region="SR", muonCR=True, decidata=False, pseudodata=False, card_name=None, verbose=0, nrho=2, npt=1, freeze_nuisances=None, fixed_signal_strength=None, ntoys=500, toys_per_job=10, condor=True, freeze_groups=None):
	if not card_name:
		# Infer card name automatically
		if region == "SR":
			if muonCR:
				card_name = "card_rhalphabet_muonCR.txt" # This is the SR+muCR datacard
			else:
				card_name = "card_rhalphabet_nomuonCR.txt" # This is the SR datacard
		elif region == "N2SR":
			card_name = "card_rhalphabet_nomuonCR.txt" # This is the SR datacard. N2SR can't use the muonCR... needs it's own CR, which you haven't made yet
		elif region == "N2CR":
			card_name = "card_rhalphabet_N2CR.txt" # This is the N2CR-only datacard
		elif region == "muCR":
			card_name = "datacard_muonCR.txt" # This is the muCR-only datacard

	datacard_directory = config.get_datacard_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region)
	job_name = "{}_{}_{}_r{}p{}".format(signal_name, jet_type, region, nrho, npt)
	if region == "SR":
		if muonCR:
			job_name += "_yesmuonCR"
		else:
			job_name += "_nomuonCR"
	if fixed_signal_strength != None:
		job_name += "_fixmu{}".format(fixed_signal_strength)
	else:
		job_name += "_floatmu"
	if decidata:
		job_name += "_ps10"
	working_directory = "{}/{}".format(gof_directory, job_name)
	os.system("mkdir -pv " + working_directory)

	cwd = os.getcwd()
	os.chdir(working_directory)
	input_files = ["{}/base.root".format(datacard_directory), "{}/rhalphabase.root".format(datacard_directory), "{}/{}".format(datacard_directory, card_name)]
	if muonCR:
		input_files.append("{}/workspace_muonCR.root".format(datacard_directory))

	# Run script
	run_script_path = "{}/run_combine_{}.sh".format(working_directory, job_name)
	run_script = open(run_script_path, 'w')
	run_script.write("#!/bin/bash\n")
	run_script.write("cp {} tmpcard.txt\n".format(card_name))
	run_script.write("mkdir -pv plots_{}\n".format(job_name))
	run_script.write("sed -i 's@{}/@@g' tmpcard.txt\n".format(datacard_directory)) # Combine cards sometimes pick up the absolute path to the workspace, which doesn't work on condor

	# Prefit
	prefit_fix_pars_str = ""
	for irho in xrange(config.analysis_parameters[jet_type]["MAX_NRHO"]+1):
		for ipt in xrange(config.analysis_parameters[jet_type]["MAX_NPT"]+1):
			if irho > nrho or ipt > npt:
				prefit_fix_pars_str += "r{}p{}:0,".format(irho, ipt)
	prefit_fix_pars_str = prefit_fix_pars_str[:-1] # Chop trailing comma
	cats = config.analysis_parameters[jet_type]["FIT_PT_BINS"]
	#if jet_type == "CA15" and region == "N2SR":
	#	cats = [1,2,3,4,5,6]
	run_script.write("python $CMSSW_BASE/python/DAZSLE/ZPrimePlusJet/prefit_workspace.py --base_path base.root --rhalphabase_path rhalphabase.root --fix_pars_rhalphabet {} --signals {} --cats {} --no_backup_original  2>&1 | tee log_prefit_{}.txt\n".format(prefit_fix_pars_str, signal_name, ",".join([str(x) for x in cats]), job_name))

	# String for freezing nuisance parameters. Both manually specified and extra polnomial terms.
	frozen_nps = []
	for irho in xrange(config.analysis_parameters[jet_type]["MAX_NRHO"]+1):
		for ipt in xrange(config.analysis_parameters[jet_type]["MAX_NPT"]+1):
			if irho > nrho or ipt > npt:
				frozen_nps.append("r{}p{}".format(irho, ipt))
	if freeze_nuisances:
		frozen_nps.extend(freeze_nuisances)
	if len(frozen_nps) >= 1:
		freeze_string = " --freezeNuisances "
		for np in frozen_nps:
			freeze_string += np + ","
		freeze_string = freeze_string[:-1]
	else:
		freeze_string = ""

	if freeze_groups:
		freeze_groups_string = " --freezeNuisanceGroups " + ",".join(freeze_groups)
	else:
		freeze_groups_string = ""

	# String for fixing signal strength
	if fixed_signal_strength != None:
		fixed_signal_strength_string = " --fixedSignalStrength {}".format(fixed_signal_strength)
	else:
		fixed_signal_strength_string = ""

	# Max likelihood fit, for plots
	run_script.write("if [ \"$1\" -eq \"0\" ];\n")
	run_script.write("\tthen\n")

	run_script.write("\tmkdir plots\n")
	command_mlf = "combine -M MaxLikelihoodFit -v 0 tmpcard.txt -n {}_mlfit  --saveNormalizations --plot --saveShapes --saveWorkspace --out plots ".format(job_name + "_maxlikelihoodfit")
	command_mlf += freeze_string
	command_mlf += freeze_groups_string
	run_script.write("\t" + command_mlf + " 2>&1 | tee log_mlfit_{}.txt\n".format(job_name))
	run_script.write("\ttar -czvf plots.tar.gz plots\n")

	# GoF of central fit
	for algorithm in algorithms:
		command_centralgof = "combine -M GoodnessOfFit --algorithm {} --rMax 20 --rMin -20 -v {} tmpcard.txt -n {}".format(algorithm, verbose, job_name + "_" + algorithm)
		command_centralgof += freeze_string
		command_centralgof += freeze_groups_string
		command_centralgof += fixed_signal_strength_string
		run_script.write("\t" + command_centralgof + " 2>&1 | tee log_gof_central_{}_{}.txt\n".format(job_name, algorithm))

	run_script.write("fi\n")

	# GoF toys, for p0 value
	command_gentoys = "combine -M GenerateOnly tmpcard.txt --rMax 20 --rMin -20 --toysFrequentist -t {} --expectSignal 0 --saveToys -n {} --seed $1".format(toys_per_job, job_name + "_gentoys_subjob$1")
	command_gentoys += freeze_string
	command_gentoys += freeze_groups_string
	run_script.write(command_gentoys + " 2>&1\n")
	toys_filename = "higgsCombine{}.GenerateOnly.mH120.$1.root".format(job_name + "_gentoys_subjob$1")

	for algorithm in algorithms:
		command_goftoys = "combine tmpcard.txt -M GoodnessOfFit --algorithm {} --rMax 20 --rMin -20 -t {} --toysFile {} -n {} --seed $1".format(algorithm, toys_per_job, toys_filename, job_name + "_" + algorithm + "_toysfit" + "_subjob$1")
		command_goftoys += freeze_string
		command_goftoys += freeze_groups_string
		command_goftoys += fixed_signal_strength_string
		run_script.write(command_goftoys + "  2>&1 | tee log_gof_toys$1_{}_{}.txt\n".format(job_name, algorithm))


	# Prepare return objects
	# - Rename workspace and datacard files to avoid overwriting
	if condor:
		run_script.write("rm -f base.root\n")
		run_script.write("rm -f rhalphabase.root\n")
		run_script.write("mv {} {}.{}\n".format(card_name, card_name, job_name))
	run_script.close()

	# csub script
	if condor:
		n_jobs = int(math.ceil(1. * ntoys / toys_per_job))
		csub_script_path = "{}/csub_{}.sh".format(working_directory, job_name)
		csub_script = open(csub_script_path, 'w')
		csub_script.write("#!/bin/bash\n")
		csub_command = "csub {} -F {} --cmssw --no_retar -n {}".format(run_script_path, ",".join(input_files), n_jobs)
		csub_script.write(csub_command + "\n")
		csub_script.close()

		# run
		os.system("source {}".format(csub_script_path))
		#print "source {}".format(csub_script_path)
	else:
		# Local run
		local_run_dir = gof_directory + "/local/" + job_name
		os.system("mkdir -pv {}".format(local_run_dir))
		os.system("cp {} {}".format(run_script_path, local_run_dir))
		for filename in input_files:
			os.system("cp {} {}".format(filename, local_run_dir))
		os.chdir(local_run_dir)
		os.system("source {}/{} 0".format(local_run_dir, os.path.basename(run_script_path)))

	os.chdir(cwd)

def get_pval(gof_central, gof_toys):
	return 1. * len([x for x in gof_toys if x > gof_central]) / len(gof_toys)


def pval_plot(central_gof, toy_gofs, plot_name):
		pval = get_pval(central_gof, toy_gofs)

		cname = "c_gof_{}".format(plot_name)
		c = TCanvas(cname, cname, 800, 600)

		# Compute some stats about the toys to determine histogram bounds
		toy_array = np.array(toy_gofs)
		toy_mean = np.mean(toy_array)
		toy_stddev = np.std(toy_array)

		h_toys = TH1D("h_toys_{}".format(cname), "Toy GOFs", 50, 0., toy_mean+4*toy_stddev)
		for gof in toy_gofs:
			h_toys.Fill(gof)

		frame = TH1D("frame_{}".format(cname), "Frame", 100, 0., max(toy_mean+4*toy_stddev, central_gof * 1.05))
		frame.GetXaxis().SetTitle("GOF")
		frame.SetMaximum(h_toys.GetMaximum())
		frame.Draw("axis")
		h_toys.Draw("hist same")

		arrow_central  = TArrow(central_gof, 0.25 * h_toys.GetMaximum(), central_gof,0)
		arrow_central.SetLineColor(kBlue+1)
		arrow_central.SetLineWidth(2)
		arrow_central.Draw("same")

		l = TLegend(0.5, 0.7, 0.85, 0.85)
		l.SetFillColor(0)
		l.SetHeader("p_{{0}}={}".format(round(pval, 2)))
		l.AddEntry(arrow_central, "Observed")
		l.AddEntry(h_toys, "Toys")
		l.Draw()

		c.SaveAs(config.paths["Fits"] + "/gof/figures/" + c.GetName() + ".pdf")


# Get the number of bins in the data histogram
# - Don't count bins where the fail region has zero events.
def get_nbins(signal_name, jet_type, region, pseudodata=False, decidata=False):
	data_file = TFile(config.get_datacard_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region) + "/base.root", "READ")
	nbins = 0
	for cat in xrange(1, 7):
		w = data_file.Get("w_fail_cat{}".format(cat))
		roodatahist = w.data("data_obs_fail_cat{}".format(cat))
		#roohistpdf = RooHistPdf("histpdf_cat{}".format(cat), "histpdf_{}".format(cat), RooArgSet(w.var("x")), rooabsdata)
		hist = roodatahist.createHistogram("h_data_obs_fail_cat{}".format(cat), w.var("x"))
		for bin in xrange(1, hist.GetNbinsX() + 1):
			if hist.GetBinContent(bin) > 0:
				nbins += 1
	return nbins

def pval_main(signal_name, jet_type, region, nrho, npt, pseudodata=False, decidata=False, nbins=-1, muonCR=False, fixed_signal_strength=None, algorithm="saturated"):
	import glob

	job_name = "{}_{}_{}_r{}p{}".format(signal_name, jet_type, region, nrho, npt)
	if region == "SR":
		if muonCR:
			job_name += "_yesmuonCR"
		else:
			job_name += "_nomuonCR"
	if fixed_signal_strength != None:
		job_name += "_fixmu{}".format(fixed_signal_strength)
	else:
		job_name += "_floatmu"
	if decidata:
		job_name += "_ps10"
	working_directory = "{}/{}".format(gof_directory, job_name)
	toy_tree = TChain("limit")
	toy_gofs = []
	for toy_file in glob.glob("{}/higgsCombine*{}_toysfit_subjob*root".format(working_directory, algorithm)): #higgsCombineSbb125_AK8_SR_r2p1_nomuonCR_floatmu_saturated_toysfit_subjob38.GoodnessOfFit.mH120.38.root
		toy_tree.Add(toy_file)
	for entry in xrange(toy_tree.GetEntries()):
		toy_tree.GetEntry(entry)
		toy_gofs.append(toy_tree.GetLeaf("limit").GetValue(0))
	print "Toy GOFs:"
	print toy_gofs

	central_file = TFile("{}/higgsCombine{}_{}.GoodnessOfFit.mH120.root".format(working_directory, job_name, algorithm))
	central_tree = central_file.Get("limit")
	central_tree.GetEntry(0)
	central_gof = central_tree.GetLeaf("limit").GetValue(0)
	print "Central GOF: {}".format(central_gof)
	central_file.Close()

	plot_name = job_name + "_" + algorithm
	pval_plot(central_gof, toy_gofs, plot_name)
	print "GOF p-value [{}] = {}".format(plot_name, get_pval(central_gof, toy_gofs))

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Run combine on condor")
	parser.add_argument("--regions", type=str, default="SR", help="SR, muCR, N2SR, or N2CR")
	parser.add_argument("--nomuonCR", action="store_false", help="Include muonCR in GOF (SR only)")
	parser.add_argument("--jet_types", type=str, default="AK8,CA15")
	parser.add_argument("--decidata", action="store_true", help="Use 1/10 data")
	parser.add_argument("--pseudodata", action="store_true", help="Use MC pseudodata")
	parser.add_argument("--card", type=str, help="Manually specify card name")
	parser.add_argument("--signals", type=str, default="Sbb125", help="all or comma-separated list of signal names")
	parser.add_argument("--fixed_signal_strength", type=float, default=None, help="Specify a fixed signal strength (otherwise, left floating)")
	parser.add_argument("--verbose", type=int, default="0", help="Combine verbosity")
	parser.add_argument("--poly_degrees", type=str, default="", help="nrho1:npt1:nrho2:npt2,...")
	parser.add_argument("--ntoys", type=int, default=500, help="Number of toys")
	parser.add_argument("--toys_per_job", type=int, default=10, help="Toys per subjob")
	parser.add_argument("--run", action="store_true", help="Run jobs on condor")
	parser.add_argument("--condor_run", action="store_true", help="Run jobs on condor")
	parser.add_argument("--pval", action="store_true", help="Get results and make p value figures")
	parser.add_argument("--freeze_nuisances", type=str, help="Freeze nuisance parameters")
	parser.add_argument("--freeze_groups", type=str, help="Freeze nuisance parameters")
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

	poly_degrees_string = {}
	if args.poly_degrees == "":
		poly_degrees_string["AK8"] = "1:1,2:1,3:1,4:1,1:2,2:2,3:2"
		poly_degrees_string["CA15"] = "1:1,2:1,3:1,4:1,5:1,6:1,1:2,2:2,3:2,4:2,5:2,1:3,2:3,3:3,4:3"
	else:
		poly_degrees_string["AK8"] = args.poly_degrees
		poly_degrees_string["CA15"] = args.poly_degrees
	poly_degrees = {"AK8":[], "CA15":[]}
	for jet_type in ["AK8", "CA15"]:
		for poly_degree_string in poly_degrees_string[jet_type].split(","):
			nrho = int(poly_degree_string.split(":")[0])
			npt = int(poly_degree_string.split(":")[1])
			poly_degrees[jet_type].append((nrho, npt))

	if args.freeze_nuisances:
		freeze_nuisances = args.freeze_nuisances.split(",")
	else:
		freeze_nuisances = None
	if args.freeze_groups:
		freeze_groups = args.freeze_groups.split(",")
	else:
		freeze_groups = None

	if args.run or args.condor_run:
		for region in args.regions.split(","):
			for jet_type in args.jet_types.split(","):
				for signal_name in signal_names[jet_type]:
					for poly_degree in poly_degrees[jet_type]:
						run_gof(signal_name, jet_type, region=region, decidata=args.decidata, pseudodata=args.pseudodata, verbose=args.verbose, nrho=poly_degree[0], npt=poly_degree[1], fixed_signal_strength=args.fixed_signal_strength, muonCR=args.nomuonCR, ntoys=args.ntoys, toys_per_job=args.toys_per_job, condor=args.condor_run, freeze_nuisances=freeze_nuisances, freeze_groups=freeze_groups)
	elif args.pval:
		for region in args.regions.split(","):
			for jet_type in args.jet_types.split(","):
				for signal_name in signal_names[jet_type]:
					for poly_degree in poly_degrees[jet_type]:
						pval_main(signal_name, jet_type, region=region, nrho=poly_degree[0], npt=poly_degree[1], pseudodata=args.pseudodata, decidata=args.decidata, muonCR=args.nomuonCR, fixed_signal_strength=args.fixed_signal_strength, algorithm="saturated")