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

ftest_directory = "/uscms/home/dryu/DAZSLE/data/Fits/ftest/"

# Submit fits via condor
def run_fits(signal_name, jet_type, region, poly_degree_pairs, pseudodata=False, decidata=False, ntoys=500, toys_per_job=10, signal_mu=0., freeze_nuisances=None, prefit=False):
	print "[run_fits] INFO : Called with args {}, {}, {}, ...".format(signal_name, jet_type, region)
	print "[run_fits] INFO : poly_degree_pairs = ",
	print poly_degree_pairs
	datacard_directory = config.get_datacard_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region)
	if region == "SR":
		card_name = "card_rhalphabet_nomuonCR.txt" # This is the SR+muCR datacard
	elif region == "N2SR":
		card_name = "card_rhalphabet_nomuonCR.txt" # This is the SR+muCR datacard
	elif region == "N2CR":
		card_name = "card_rhalphabet_N2CR.txt" # This is the N2CR-only datacard


	n_jobs = int(math.ceil(1. * ntoys / toys_per_job))
	hadd_script_paths = []
	for poly_degree_pair in poly_degree_pairs:
		nrho1 = poly_degree_pair[0][0]
		npt1 = poly_degree_pair[0][1]
		nrho2 = poly_degree_pair[1][0]
		npt2 = poly_degree_pair[1][1]

		# Freeze higher orders
		freeze_poly_string1 = ""
		for irho in xrange(config.analysis_parameters[jet_type]["MAX_NRHO"]+1):
			for ipt in xrange(config.analysis_parameters[jet_type]["MAX_NPT"]+1):
				if irho > nrho1 or ipt > npt1:
					freeze_poly_string1 += "r{}p{},".format(irho, ipt)
		freeze_poly_string1 = freeze_poly_string1[:-1]
		if freeze_nuisances:
			freeze_nuisances1 = freeze_nuisances + "," + freeze_poly_string1
		else:
			freeze_nuisances1 = freeze_poly_string1

		freeze_poly_string2 = ""
		for irho in xrange(config.analysis_parameters[jet_type]["MAX_NRHO"]+1):
			for ipt in xrange(config.analysis_parameters[jet_type]["MAX_NPT"]+1):
				if irho > nrho2 or ipt > npt2:
					freeze_poly_string2 += "r{}p{},".format(irho, ipt)
		freeze_poly_string2 = freeze_poly_string2[:-1]
		if freeze_nuisances:
			freeze_nuisances2 = freeze_nuisances + "," + freeze_poly_string2
		else:
			freeze_nuisances2 = freeze_poly_string2



		job_name1 = "ftestfits1_{}_{}_{}_r{}p{}".format(signal_name, jet_type, region, nrho1, npt1)
		job_name2 = "ftestfits2_{}_{}_{}_r{}p{}".format(signal_name, jet_type, region, nrho2, npt2)
		fit_directory = config.get_ftest_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region, nrho1=nrho1, npt1=npt1, nrho2=nrho2, npt2=npt2) + "/condor/"
		os.system("mkdir -pv {}".format(fit_directory))
		cwd = os.getcwd()
		os.chdir(fit_directory)

		# Run script
		run_script_path = "{}/run_ftest_fits.sh".format(fit_directory)
		run_script = open(run_script_path, "w")
		run_script.write("#!/bin/bash\n")
		run_script.write("cp {} tmpcard.txt\n".format(card_name))
		run_script.write("sed -i 's@{}/@@g' tmpcard.txt\n".format(datacard_directory)) # Combine cards sometimes pick up the absolute path to the workspace, which doesn't work on condor

		if prefit:
			prefit_fix_pars_str = ""
			for irho in xrange(config.analysis_parameters[jet_type]["MAX_NRHO"]+1):
				for ipt in xrange(config.analysis_parameters[jet_type]["MAX_NPT"]+1):
					if irho > nrho1 or ipt > npt1:
						prefit_fix_pars_str += "r{}p{}:0,".format(irho, ipt)
			prefit_fix_pars_str = prefit_fix_pars_str[:-1] # Chop trailing comma
			run_script.write("python $CMSSW_BASE/python/DAZSLE/ZPrimePlusJet/prefit_workspace.py --base_path base.root --rhalphabase_path rhalphabase.root --fix_pars_rhalphabet {} --signals {}\n".format(prefit_fix_pars_str, signal_name))


		command_fit1 = "combine -M MaxLikelihoodFit -v 0 -t -1 --toysFreq tmpcard.txt -n {}_mlfit  --saveNormalizations --plot --saveShapes --saveWorkspace --out plots1".format(job_name1)
		if freeze_nuisances1:
			command_fit1 += " --freezeNuisances {}".format(freeze_nuisances1)
		run_script.write("if [ \"$1\" -eq \"0\" ];\n\tthen\n\tmkdir plots1\n\t" + command_fit1 + "\n\ttar -czvf plots1.tar.gz plots1\nfi\n")

		command_fit2 = "combine -M MaxLikelihoodFit -v 0 -t -1 --toysFreq tmpcard.txt -n {}_mlfit  --saveNormalizations --plot --saveShapes --saveWorkspace --out plots2".format(job_name2)
		if freeze_nuisances2:
			command_fit2 += " --freezeNuisances {}".format(freeze_nuisances2)
		run_script.write("if [ \"$1\" -eq \"0\" ];\n\tthen\n\tmkdir plots2\n\t" + command_fit2 + "\n\ttar -czvf plots2.tar.gz plots2\nfi\n")


		command_gof1 = "combine -M GoodnessOfFit {}  --rMax 20 --rMin -20 --fixedSignalStrength 0 --algorithm saturated -n {} ".format("tmpcard.txt", job_name1 + "_centralfit")
		if freeze_nuisances1:
			command_gof1 += " --freezeNuisances {}".format(freeze_nuisances1)
		run_script.write("if [ \"$1\" -eq \"0\" ];\n\tthen\n\t" + command_gof1 + "\nfi\n")

		command_gof2 = "combine -M GoodnessOfFit {}  --rMax 20 --rMin -20 --fixedSignalStrength 0 --algorithm saturated -n {} ".format("tmpcard.txt", job_name2 + "_centralfit")
		if freeze_nuisances2:
			command_gof2 += " --freezeNuisances {}".format(freeze_nuisances2)
		run_script.write("if [ \"$1\" -eq \"0\" ];\n\tthen\n\t" + command_gof2 + "\nfi\n")

		command_gentoys = "combine -M GenerateOnly {} --rMax 20 --rMin -20 --toysFrequentist -t {} --expectSignal {} --saveToys -n {} --seed $1".format("tmpcard.txt", toys_per_job, signal_mu, job_name1 + "_gentoys_subjob$1")
		if freeze_nuisances1:
			command_gentoys += " --freezeNuisances {}".format(freeze_nuisances1)
		run_script.write(command_gentoys + "\n")
		toys_filename = "higgsCombine{}.GenerateOnly.mH120.$1.root".format(job_name1 + "_gentoys_subjob$1")

		command_fittoys1 = "combine -M GoodnessOfFit {} --rMax 20 --rMin -20 -t {} --fixedSignalStrength 0 --toysFile {} --algorithm saturated -n {} --seed $1".format("tmpcard.txt", toys_per_job, toys_filename, job_name1 + "_toysfit" + "_subjob$1")
		if freeze_nuisances1:
			command_fittoys1 +=  " --freezeNuisances {}".format(freeze_nuisances1)
		run_script.write(command_fittoys1 + "\n")

		command_fittoys2 = "combine -M GoodnessOfFit {} --rMax 20 --rMin -20 -t {} --fixedSignalStrength 0 --toysFile {} --algorithm saturated -n {} --seed $1".format("tmpcard.txt", toys_per_job, toys_filename, job_name2 + "_toysfit" + "_subjob$1")
		if freeze_nuisances2:
			command_fittoys2 +=  " --freezeNuisances {}".format(freeze_nuisances2)
		run_script.write(command_fittoys2 + "\n")
		run_script.close()

		files_to_transfer = ["{}/{}".format(datacard_directory, card_name), "{}/base.root".format(datacard_directory), "{}/rhalphabase.root".format(datacard_directory)]
		csub_command = "csub {} --cmssw --no_retar -m 4000 -F {} -n {}".format(run_script_path, ",".join(files_to_transfer), n_jobs)
		csub_script_path = "{}/csub_submit.sh".format(fit_directory)
		csub_script = open(csub_script_path, "w")
		csub_script.write("#!/bin/bash\n")
		csub_script.write(csub_command + "\n")
		csub_script.close()

		os.system("source {}".format(csub_script_path))
		#print "cat {}".format(csub_script_path)
		
		# Postprocessing script
		hadd_script = open("{}/hadd.sh".format(fit_directory), "w")
		hadd_script.write("#!/bin/bash\n")
		hadd_script.write("hadd -f {}/../toyfits1.root {}/higgsCombine{}_toysfit_subjob*.GoodnessOfFit.mH120.*.root\n".format(fit_directory, fit_directory, job_name1))
		hadd_script.write("hadd -f {}/../toyfits2.root {}/higgsCombine{}_toysfit_subjob*.GoodnessOfFit.mH120.*.root\n".format(fit_directory, fit_directory, job_name2))
		hadd_script.write("cp {}/higgsCombine{}_centralfit.GoodnessOfFit.mH120.root  {}/../centralfit1.root\n".format(fit_directory, job_name1, fit_directory))
		hadd_script.write("cp {}/higgsCombine{}_mlfit.MaxLikelihoodFit.mH120.root  {}/../mlfit1.root\n".format(fit_directory, job_name1, fit_directory))
		hadd_script.write("cp {}/higgsCombine{}_centralfit.GoodnessOfFit.mH120.root  {}/../centralfit2.root\n".format(fit_directory, job_name2, fit_directory))
		hadd_script.write("cp {}/higgsCombine{}_mlfit.MaxLikelihoodFit.mH120.root  {}/../mlfit2.root\n".format(fit_directory, job_name2, fit_directory))
		hadd_script.write("cp {}/plots1.tar.gz {}/../plots1.tar.gz \n".format(fit_directory, fit_directory))
		hadd_script.write("cp {}/plots2.tar.gz {}/../plots2.tar.gz \n".format(fit_directory, fit_directory))
		hadd_script.write("cd {}/..\n".format(fit_directory))
		hadd_script.write("tar -xzvf plots1.tar.gz\n")
		hadd_script.write("tar -xzvf plots2.tar.gz\n")
		hadd_script.write("cd -\n")
		hadd_script.close()
		hadd_script_paths.append("{}/hadd.sh".format(fit_directory))
	master_hadd_script = open("{}/ftest/master_hadd_{}_{}_{}_{}.sh".format(config.paths["Fits"], signal_name, jet_type, region, int(time.time())), "w")
	for hadd_script_path in hadd_script_paths:
		master_hadd_script.write("source {}\n".format(hadd_script_path))
	master_hadd_script.close()


# Stuff for extracting and calculating F statistics
def get_fit_likelihoods(signal_name, jet_type, region, nrho1, npt1, nrho2, npt2, pseudodata=False, decidata=False):
	fit_directory = config.get_ftest_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region, nrho1=nrho1, npt1=npt1, nrho2=nrho2, npt2=npt2)
	f_central1 = TFile("{}/centralfit1.root".format(fit_directory), "READ")
	if not f_central1.IsOpen():
		print "[get_fit_likelihoods] ERROR : Couldn't open file {}".format("{}/centralfit1.root".format(fit_directory))
		sys.exit(1)
	t_central1 = f_central1.Get("limit")
	if not t_central1:
		print "[get_fit_likelihoods] ERROR : Couldn't find tree limit in file {}".format(f_central1.GetPath())
		sys.exit(1)
	t_central1.GetEntry(0)
	limit_central1 = t_central1.GetLeaf("limit").GetValue(0)
	f_central1.Close()

	f_central2 = TFile("{}/centralfit2.root".format(fit_directory), "READ")
	t_central2 = f_central2.Get("limit")
	if not t_central2:
		print "[get_fit_likelihoods] ERROR : Couldn't find tree limit in file {}".format(f_central2.GetPath())
		sys.exit(1)
	t_central2.GetEntry(0)
	limit_central2 = t_central2.GetLeaf("limit").GetValue(0)
	f_central2.Close()

	f_toys1 = TFile("{}/toyfits1.root".format(fit_directory), "READ")
	t_toys1 = f_toys1.Get("limit")
	limit_toys1 = []
	for entry in xrange(t_toys1.GetEntriesFast()):
		t_toys1.GetEntry(entry)
		limit_toys1.append(t_toys1.GetLeaf("limit").GetValue(0))
	f_toys1.Close()

	f_toys2 = TFile("{}/toyfits2.root".format(fit_directory), "READ")
	t_toys2 = f_toys2.Get("limit")
	limit_toys2 = []
	for entry in xrange(t_toys2.GetEntriesFast()):
		t_toys2.GetEntry(entry)
		limit_toys2.append(t_toys2.GetLeaf("limit").GetValue(0))
	f_toys2.Close()

	return limit_central1, limit_toys1, limit_central2, limit_toys2

def get_fstat(likelihood1, likelihood2, ndf1, ndf2):
	return (likelihood1 - likelihood2) / (ndf1 - ndf2) / (likelihood2 / ndf2)

def get_fstats(signal_name, jet_type, region, pseudodata=False, decidata=False, nbins=-1, nrho1=0, npt1=0, nrho2=0, npt2=0):
	if nbins < 0:
		print "[get_fstats] ERROR : Must specify number of bins (arg nbins)"
		sys.exit(1)
	l_central1, l_toys1, l_central2, l_toys2 = get_fit_likelihoods(signal_name, jet_type, region, nrho1=nrho1, npt1=npt1, nrho2=nrho2, npt2=npt2, pseudodata=pseudodata, decidata=decidata)

	ndf1 = nbins - nrho1*npt1
	ndf2 = nbins - nrho2*npt2

	fstat_central = get_fstat(l_central1, l_central2, ndf1, ndf2)
	fstat_toys = []
	for itoy in xrange(len(l_toys1)):
		fstat_toys.append(get_fstat(l_toys1[itoy], l_toys2[itoy], ndf1, ndf2))

	return fstat_central, fstat_toys

def get_pval(fstat_central, fstat_toys):
	return 1. * len([x for x in fstat_toys if x > fstat_central]) / len(fstat_toys)

def pval_table(signal_name, jet_type, region, poly_degree_pairs, pseudodata=False, decidata=False, nbins=-1):
	print "[pval_table] INFO : Called with args signal_name={}, jet_type={}, region={}, pseudodata={}, decidata={}, nbins={}, and poly_degree_pairs=...".format(signal_name, jet_type, region, pseudodata, decidata, nbins)
	print poly_degree_pairs

	table_file = open(config.get_ftest_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region, nrho1=2, npt1=1, nrho2=3, npt2=1) + "/../../pvals.tex", "w")
	table_file.write("\\usepackage{collcell,array}\n")
	table_file.write("\\newcommand{\\shiftdown}[1]{\\smash{\\raisebox{-.5\\normalbaselineskip}{#1}}}\n")
	table_file.write("\\newcolumntype{C}{>{\\collectcell\\shiftdown}c<{\\endcollectcell}}\n\n")

	table_file.write("\\begin{table}\n")
	table_file.write("\t\\begin{tabular}{|c|C|C|}\n")
	table_file.write("\t\t\\hline\n")
	table_file.write("\t\t$(n_{\\rho},\\,n_{\\pt})$\t&\t\\multicolumn{1}{c}{$F$}\t&\t\\multicolumn{1}{c}{$p_0$}\t\n")
	table_file.write("\t\t\\hline\n")
	last_fstat = -1e20
	for poly_degree_pair in poly_degree_pairs:
		nrho1 = poly_degree_pair[0][0]
		npt1 = poly_degree_pair[0][1]
		nrho2 = poly_degree_pair[1][0]
		npt2 = poly_degree_pair[1][1]

		fstat_central, fstat_toys = get_fstats(signal_name, jet_type, region, pseudodata=pseudodata, decidata=decidata, nbins=nbins, nrho1=nrho1, npt1=npt1, nrho2=nrho2, npt2=npt2)
		pval = get_pval(fstat_central, fstat_toys)
		row = "\t\t$({},\,{})$\t&\t${}$\t&\t${}$\t\n".format(nrho1, npt1, round(fstat_central, 3), round(pval, 3))
		table_file.write(row)
		table_file.write("\t\t\\hline\n")
		table_file.write("\t\t$({},\,{})$\n".format(nrho2, npt2))
		table_file.write("\t\t\\hline\n\t\t\\hline\n")
		print "{},{}\tF\tp0".format(nrho1, npt1)
		print "\t\t{}\t{}".format(fstat_central, pval)
		print "{},{}\n".format(nrho2, npt2)
	table_file.write("\t\\end{tabular}\n")
	table_file.write("\t\\caption{{$F$ test results for {}/{}}}\n".format(region, jet_type))
	table_file.write("\\end{table}\n")
	table_file.close()

def pval_plots(signal_name, jet_type, region, poly_degree_pairs, pseudodata=False, decidata=False, nbins=-1):
	print "[pval_table] INFO : Called with args signal_name={}, jet_type={}, region={}, pseudodata={}, decidata={}, nbins={}, and poly_degree_pairs=...".format(signal_name, jet_type, region, pseudodata, decidata, nbins)
	print poly_degree_pairs
	for poly_degree_pair in poly_degree_pairs:
		nrho1 = poly_degree_pair[0][0]
		npt1 = poly_degree_pair[0][1]
		nrho2 = poly_degree_pair[1][0]
		npt2 = poly_degree_pair[1][1]

		fstat_central, fstat_toys = get_fstats(signal_name, jet_type, region, pseudodata=pseudodata, decidata=decidata, nbins=nbins, nrho1=nrho1, npt1=npt1, nrho2=nrho2, npt2=npt2)
		pval = get_pval(fstat_central, fstat_toys)

		cname = "c_r{}p{}_vs_r{}p{}_fstat_{}_{}_{}".format(nrho1, npt1, nrho2, npt2, signal_name, jet_type, region)
		if pseudodata:
			cname += "_pseudodata"
		if decidata:
			cname += "_ps10"
		c = TCanvas(cname, cname, 800, 600)

		# Compute some stats about the toys to determine histogram bounds
		toy_array = np.array(fstat_toys)
		toy_mean = np.mean(toy_array)
		toy_stddev = np.std(toy_array)

		h_toys = TH1D("h_toys_{}".format(cname), "Toy F stats", 50, 0., toy_mean+4*toy_stddev)
		for fstat in fstat_toys:
			h_toys.Fill(fstat)

		frame = TH1D("frame_{}".format(cname), "Frame", 100, 0., max(toy_mean+4*toy_stddev, fstat_central * 1.05))
		frame.GetXaxis().SetTitle("F stat")
		frame.SetMaximum(h_toys.GetMaximum())
		frame.Draw("axis")
		h_toys.Draw("hist same")

		arrow_central  = TArrow(fstat_central, 0.25 * h_toys.GetMaximum(), fstat_central,0)
		arrow_central.SetLineColor(kBlue+1)
		arrow_central.SetLineWidth(2)
		arrow_central.Draw("same")

		l = TLegend(0.5, 0.7, 0.85, 0.85)
		l.SetFillColor(0)
		l.SetHeader("p_{{0}}={}".format(round(pval, 2)))
		l.AddEntry(arrow_central, "Observed")
		l.AddEntry(h_toys, "Toys")
		l.Draw()

		c.SaveAs(config.paths["Fits"] + "/ftest/figures/" + c.GetName() + ".pdf")


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


if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="One stop shop for F tests")
	action_group = parser.add_mutually_exclusive_group() 
	action_group.add_argument('--fits', action="store_true", help="Run fits")
	action_group.add_argument('--pval', action="store_true", help="Compute F statistics and p-values")
	action_group.add_argument('--plots', action="store_true", help="Plot fit results")

	parser.add_argument("--poly_degrees", type=str, default="", help="nrho1:npt1:nrho2:npt2,...")
	parser.add_argument("--signals", type=str, default="Sbb125", help="Signal names to run")
	parser.add_argument("--jet_type", type=str, default="AK8", help="AK8 or CA15")
	parser.add_argument("--region", type=str, default="SR", help="SR, N2SR, or N2CR")
	parser.add_argument("--pseudodata", action="store_true", help="Run over QCD pseudodata rather than data")
	parser.add_argument("--decidata", action="store_true", help="Run over 1/10 data")
	parser.add_argument("--mu", type=float, default=0., help="Signal mu for toy generation")
	parser.add_argument("--ntoys", type=int, default=500, help="Number of toys")
	parser.add_argument("--toys_per_job", type=int, default=10, help="Toys per subjob")
	parser.add_argument("--freeze_nuisances", type=str, help="Pass to --freezeNuisances argument of combine")
	parser.add_argument("--prefit", action="store_true", help="Prefit the workspace (recommended)")
	args = parser.parse_args()

	signals = args.signals.split(",")
	if args.poly_degrees == "":
		if args.jet_type == "AK8":
			poly_degrees_string = "1:1:2:1,2:1:3:1,3:1:4:1,2:1:2:2,3:1:3:2"
		elif args.jet_type == "CA15":
			poly_degrees_string = "2:1:3:1,3:1:4:1,4:1:5:1,5:1:6:1,5:1:5:2,5:1:4:2,4:2:5:2,5:2:6:2"
	else:
		poly_degrees_string = args.poly_degrees
	poly_degree_pairs = []
	for poly_degree_string in poly_degrees_string.split(","):
		nrho1 = int(poly_degree_string.split(":")[0])
		npt1 = int(poly_degree_string.split(":")[1])
		nrho2 = int(poly_degree_string.split(":")[2])
		npt2 = int(poly_degree_string.split(":")[3])
		if nrho1*npt1 >= nrho2*npt2:
			print "[ftest_xbb] ERROR : Invalid degree pair (nrho1, npt1)=({}, {}), (nrho2, npt2)=({}, {})".format(nrho1, npt1, nrho2, npt2)
		poly_degree_pairs.append(((nrho1, npt1), (nrho2, npt2)))

	if args.fits:
		for signal in signals:
			run_fits(signal, args.jet_type, args.region, poly_degree_pairs=poly_degree_pairs, pseudodata=args.pseudodata, decidata=args.decidata, ntoys=args.ntoys, signal_mu=args.mu, freeze_nuisances=args.freeze_nuisances, toys_per_job=args.toys_per_job, prefit=args.prefit)

	if args.pval:
		for signal in signals:
			nbins = get_nbins(signal, args.jet_type, args.region, pseudodata=args.pseudodata, decidata=args.decidata)
			pval_table(signal, args.jet_type, args.region, poly_degree_pairs, pseudodata=args.pseudodata, decidata=args.decidata, nbins=nbins)
			pval_plots(signal, args.jet_type, args.region, poly_degree_pairs, pseudodata=args.pseudodata, decidata=args.decidata, nbins=nbins)
	print "Done."

