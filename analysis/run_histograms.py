
if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Produce and plot ieta-iphi histograms to look for buggy events')
	input_group = parser.add_mutually_exclusive_group() 
	input_group.add_argument('--supersamples', type=str, help="List of supersamples to fun, comma separated")
	input_group.add_argument('--files', type=str, help="Input file name(s), comma separated")
	input_group.add_argument('--files_txt', type=str, help="Text files with lists of inputs, comma separated")
	parser.add_argument('--output_folder', type=str, help="Output folder")
	parser.add_argument('--label', type=str, help="Output file label")
	parser.add_argument('--jet_type', type=str, default="AK8", help="AK8 or CA15")
	parser.add_argument('--skim_inputs', action='store_true', help="Run over skim inputs")
	args = parser.parse_args()

	# Make list of input files
	samples = []
	sample_files = {}
	if args.supersamples:
		supersamples = args.supersamples.split(",")
		samples = [] 
		for supersample in supersamples:
			samples.extend(config.samples[supersample])
			for sample in config.samples[supersample]:
				if args.skim_inputs:
					sample_files[sample] = config.skims[sample]
				else:
					sample_files[sample] = config.sklims[sample]
	elif args.files:
		sample_files[args.label] = args.files.split(",")
		samples.append(args.label)
	elif args.files_txt:
		samples.append(args.label)
		sample_files[args.label] = []
		with open(args.files_txt, 'r') as f:
			for line in f:
				sample_files[args.label].append(line.strip())

	for sample in samples:
		print "\n *** Running sample {}".format(sample)
		if "Sbb" in sample or args.skim_inputs or "ZPrime" in sample:
			tree_name = "Events"
		else:
			tree_name = "otree"

		# Sanity check: make sure tree exists in file
		for filename in sample_files[sample]:
			print "[event_selection_histograms] INFO : Checking contents of file {}".format(filename)
			f = ROOT.TFile.Open(filename, "READ")
			t = f.Get(tree_name)
			if not t:
				if tree_name == "otree":
					backup_tree_name = "Events"
				else:
					backup_tree_name = "otree"
				t_backup = f.Get(backup_tree_name)
				if t_backup:
					print "[setup_limits] WARNING : Didn't find tree {} in input file, but did find {}. Changing the tree name, but try to fix this.".format(tree_name, backup_tree_name)
					tree_name = backup_tree_name
				else:
					print "[setup_limits] ERROR : Didn't find tree {} in input file, nor {}. Quitting!".format(tree_name, backup_tree_name)
					sys.exit(1)
			# Check that the "NEvents" histogram is present
			h_NEvents = f.Get("NEvents")
			if not h_NEvents:
				if "data" in sample:
					print "[setup_limits] ERROR : NEvents histogram in not in this file! It is probably corrupt. This is data, so this problem is fatal."
					sys.exit(1)
				else:
					print "[setup_limits] WARNING : NEvents histogram in not in this file! It is probably corrupt. This is MC, so I am skipping the file. But, you probably want to remove from the input list."
					sample_files[sample].remove(filename)
			
		limit_histogrammer = Histograms(sample, tree_name=tree_name, jet_type=args.jet_type)
		output_file_basename ="InputHistograms_{}_{}.root".format(sample, args.jet_type) 
		if args.output_folder:
			limit_histogrammer.set_output_path("{}/{}".format(args.output_folder, output_file_basename))
		else:
			limit_histogrammer.set_output_path("/uscms/home/dryu/DAZSLE/data/LimitSetting/{}".format(output_file_basename))
		for filename in sample_files[sample]:
			print "Input file {}".format(filename)
			limit_histogrammer.add_file(filename)
		#limit_histogrammer.set_jet_type(args.jet_type)
		if "JetHTRun2016" in sample or "SingleMuRun2016" in sample:
			limit_histogrammer.set_data_source("data")
		else:
			limit_histogrammer.set_data_source("simulation")
		if "ps10" in sample:
			limit_histogrammer.set_prescale(10)
		limit_histogrammer.start()
		limit_histogrammer.run()
		limit_histogrammer.finish()

