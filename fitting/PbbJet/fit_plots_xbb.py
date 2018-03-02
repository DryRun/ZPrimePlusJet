import os
import sys
import math
from ROOT import *
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/slc6_amd64_gcc491/libDAZSLEPhiBBPlusJet.so"))
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so"))
import DAZSLE.PhiBBPlusJet.style as style
import DAZSLE.ZPrimePlusJet.xbb_config as config
seaborn = Root.SeabornInterface()
seaborn.Initialize()

def fitb_plots(signal_name, jet_type, region, nrho, npt, pseudodata=False, decidata=False, logy=False):
	print "[fitb_plot] DEBUG : Getting data histogram from {}".format(config.get_datacard_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region) + "/base.root")
	data_file = TFile(config.get_datacard_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region) + "/base.root", "READ")
	# mlfit plots file = e.g. ~/DAZSLE/data/LimitSetting/Xbb_inputs/N2CR_CA15/cards_mcstat/Sbb125/plots_combine_Sbb125_CA15_MaxLikelihoodFit_N2CR_r5p1/mlfitcombine_Sbb125_CA15_MaxLikelihoodFit_N2CR_r5p1.root
	fit_jobname = "{}_{}_MaxLikelihoodFit_{}_r{}p{}".format(signal_name, jet_type, region,nrho, npt)

	# Untar plots dir
	cwd = os.getcwd()
	os.chdir(config.get_limit_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region))
	os.system("tar -xzvf " + config.get_limit_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region) + "/plots_combine_{}.tar.gz".format(fit_jobname))
	os.chdir(cwd)

	print "[fitb_plot] DEBUG : Getting fitted s and b from " + config.get_limit_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region) + "/plots_combine_{}/mlfitcombine_{}.root".format(fit_jobname, fit_jobname)
	fit_file = TFile(config.get_limit_directory(signal_name, jet_type, qcd=pseudodata, decidata=decidata, region=region) + "/plots_combine_{}/mlfitcombine_{}.root".format(fit_jobname, fit_jobname), "READ")
	if not fit_file.IsOpen():
		print "[fitb_plot] ERROR : Couldn't open file {}".format(fit_file.GetPath())
		return -1
	for cat in config.analysis_parameters[jet_type]["FIT_PT_BINS"]:
		for what in ["pass", "fail"]:
			cname = "c_msdfitb_{}_cat{}_{}_{}_{}_{}".format(what, cat, jet_type, region, nrho, npt)
			w_data = data_file.Get("w_{}_cat{}".format(what, cat))
			if not w_data:
				print "[fitb_plot] ERROR : Couldn't load object {} from file {}".format("w_{}_cat{}".format(what, cat), data_file.GetPath())
				return -1
			data_hist = w_data.data("data_obs_{}_cat{}".format(what, cat)).createHistogram("h_data_obs_{}_cat{}".format(what, cat), w_data.var("x"))

			if not fit_file.Get("shapes_fit_b"):
				print "[fitb_plot] ERROR : Couldn't load shapes_fit_b from file {}".format(fit_file.GetPath())
				return -1
			if not fit_file.Get("shapes_fit_b").Get("cat{}_{}_cat{}".format(cat, what, cat)):
				print "[fitb_plot] ERROR : Couldn't load shapes_fit_b/{} from file {}".format("cat{}_{}_cat{}".format(cat, what, cat), fit_file.GetPath())
				return -1
			fit_hist_dir = fit_file.Get("shapes_fit_b").Get("cat{}_{}_cat{}".format(cat, what, cat))
			if not fit_hist_dir:
				print "[fitb_plot] ERROR : Couldn't find directory shapes_fit_b/{} in file {}".format("cat{}_{}_cat{}".format(cat, what, cat), fit_file.GetPath())
			background_hists = {}
			total_background_hist = None
			for key in fit_hist_dir.GetListOfKeys():
				if not "TH1" in key.GetClassName():
					continue
				if key.GetName() == "total_background":
					total_background_hist = key.ReadObj()
				elif key.GetName() in ["total", "total_signal", "total_covar"] or "Sbb" in key.GetName():
					continue
				elif key.GetName() in ["hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"]:
					if "hbb" in background_hists:
						if background_hists["hbb"].GetNbinsX() != key.ReadObj().GetNbinsX():
							print "WARNING : Existing hbb hist has {} bins, while {} has {}".format(background_hists["hbb"].GetNbinsX(), key.GetName(), key.ReadObj().GetNbinsX())
							sys.exit(1)
						background_hists["hbb"].Add(key.ReadObj())
					else:
						background_hists["hbb"] = key.ReadObj().Clone()
						background_hists["hbb"].SetName("hbb_" + cname)
				else:
					background_hists[key.GetName()] = key.ReadObj()

			# Scale background histograms by bin width -- the combine outputs are PDFs normalized to bin width 1
			for name, hist in background_hists.iteritems():
				hist.Scale(7)
			total_background_hist.Scale(7)

			if pseudodata:
				cname += "_pseudodata"
			if decidata:
				cname += "_ps10"
			if logy:
				cname += "_logy"
			c = TCanvas(cname, cname, 1000, 1000)
			top = TPad("top_{}".format(cname), "top_{}".format(cname), 0., 0.5, 1., 1.)
			top.SetBottomMargin(0.03)
			top.Draw()
			top.cd()

			frame = TH1D("frame_top_{}".format(cname), "frame", 100, 25, 500.)
			if logy:
				frame.SetMinimum(1.)
				frame.SetMaximum(data_hist.GetMaximum() * 10.)
			else:
				frame.SetMinimum(0.)
				frame.SetMaximum(data_hist.GetMaximum() * 1.3)
			frame.GetXaxis().SetLabelSize(0)
			frame.GetXaxis().SetTitleSize(0)
			frame.GetYaxis().SetTitle("Events / bin")
			frame.Draw()

			background_list = sorted(background_hists.keys(), key=lambda x: background_hists[x].Integral())
			background_stack = THStack("bkgd_stack_{}".format(cname), "bkgd_stack")
			for background in background_list:
				print "[debug] Background {} has bins {}, {}, {}".format(background, background_hists[background].GetNbinsX(), background_hists[background].GetXaxis().GetXmin(), background_hists[background].GetXaxis().GetXmax())
				background_hists[background].SetFillStyle(1001)
				background_hists[background].SetFillColor(style.background_colors[background])
				background_hists[background].SetLineColor(1)
				#print "[debug] Adding {} to the stack now".format(background)
				background_stack.Add(background_hists[background])
				#print "[debug]\t...done."
			#print "[debug] Drawing background stack..."
			background_stack.Draw("hist same")
			#print "[debug]\t...done."

			#print "[debug] Drawing total background histogram"
			total_background_hist.SetLineColor(1)
			total_background_hist.SetLineStyle(1)
			total_background_hist.SetLineWidth(2)
			total_background_hist.Draw("hist same")

			#print "[debug] Drawing data histogram"
			data_hist.SetMarkerStyle(20)
			data_hist.SetMarkerColor(1)
			data_hist.SetMarkerSize(1)
			data_hist.Draw("p same")

			#print "[debug] Making legend"
			l = TLegend(0.6, 0.4, 0.85, 0.85)
			l.SetFillColor(0)
			l.SetBorderSize(0)
			l.AddEntry(data_hist, "Data 2016", "p")
			l.AddEntry(total_background_hist, "Total Bkgd", "l")
			for background in reversed(background_list):
				l.AddEntry(background_hists[background], style.background_legends[background], "fl")
			l.Draw()

			c.cd()
			bottom = TPad("bottom_{}".format(cname), "bottom_{}".format(cname), 0., 0.0, 1., 0.5)
			bottom.SetTopMargin(0.03)
			bottom.SetBottomMargin(0.2)
			bottom.Draw()
			bottom.cd()

			frame_ratio = TH1D("frame_ratio_{}".format(cname), "frame_ratio", 100, 25, 500.)
			frame_ratio.GetXaxis().SetTitle("m_{SD} [GeV]")
			frame_ratio.GetYaxis().SetTitle("#frac{data - bkgd}{#sqrt{bkgd}}")
			frame_ratio.SetMinimum(-4.)
			frame_ratio.SetMaximum(4.)
			frame_ratio.Draw()

			print "[debug] Total data / bkgd = {} / {} = {}".format(data_hist.Integral(), total_background_hist.Integral(), data_hist.Integral() / total_background_hist.Integral())
			ratio = data_hist.Clone()
			ratio.Reset()
			for bin in xrange(1, data_hist.GetNbinsX() + 1):
				data = data_hist.GetBinContent(bin)
				bkgd = total_background_hist.GetBinContent(bin)
				if bkgd > 0:
					pull = (data - bkgd) / math.sqrt(bkgd)
				else:
					pull = 0.
				ratio.SetBinContent(bin, pull)
			ratio.SetMarkerStyle(20)
			ratio.SetMarkerColor(1)
			ratio.Draw("hist same")

			c.cd()
			c.SaveAs(config.paths["Fits"] + "/figures/" + c.GetName() + ".pdf")
			SetOwnership(top, False)
			SetOwnership(bottom, False)
			SetOwnership(c, False)
			total_background_hist.IsA().Destructor(total_background_hist)
	return 0

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Plot mSD fits")
	parser.add_argument("--signals", type=str, default="Sbb125", help="Signal names to run")
	parser.add_argument("--jet_types", type=str, default="AK8,CA15", help="AK8 or CA15")
	parser.add_argument("--regions", type=str, default="N2CR", help="SR, N2SR, or N2CR")
	parser.add_argument("--pseudodata", action="store_true", help="Run over QCD pseudodata rather than data")
	parser.add_argument("--decidata", action="store_true", help="Run over 1/10 data")
	parser.add_argument("--poly_degrees", type=str, default="", help="nrho1:npt1,nrho2:npt2,...")
	parser.add_argument("--npts", type=str, help="npt")
	args = parser.parse_args()

	poly_degrees_strings = {}
	if args.poly_degrees == "":
		poly_degrees_strings["AK8"] = "2:1,3:1,2:2"
		poly_degrees_strings["CA15"] = "2:1,3:1,4:1,5:1,6:1,4:2,5:2,1:3,2:3,3:3,4:3"
	else:
		poly_degrees_strings["AK8"] = args.poly_degrees
		poly_degrees_strings["CA15"] = args.poly_degrees
	poly_degree_pairs = {"AK8":[], "CA15":[]}
	for jet_type in args.jet_types.split(","):
		for poly_degree_string in poly_degrees_strings[jet_type].split(","):
			nrho = int(poly_degree_string.split(":")[0])
			npt = int(poly_degree_string.split(":")[1])
			poly_degree_pairs[jet_type].append((nrho, npt))

	failed_jobs = []
	for signal in args.signals.split(","):
		for jet_type in args.jet_types.split(","):
			for region in args.regions.split(","):
				for poly_degree in poly_degree_pairs[jet_type]:
					ret = fitb_plots(signal, jet_type, region, poly_degree[0], poly_degree[1], pseudodata=args.pseudodata, decidata=args.decidata, logy=False)
					if ret == -1:
						failed_jobs.append([signal, jet_type, region, poly_degree])
	print "Failed jobs:"
	print failed_jobs
