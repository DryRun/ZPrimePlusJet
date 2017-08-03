import os

analysis_parameters = {}
analysis_parameters["AK8"] = {
	"MASS_BINS":80,
	"MSD":[40, 600],
	"RHO":[-5.5, -2.0],
	"PT":[450., 1000.],
	"DCSV":0.9,
	"N2DDT":0.,
	"PT_BINS":[450., 500.,550.,600.,675.,800.,1000.],
	"BB_SF":0.91,
	"BB_SF_ERR":0.03,
	"V_SF":0.993,
	"V_SF_ERR":0.043,
}
analysis_parameters["CA15"] = {
	"MASS_BINS":80,
	"MSD":[40, 600],
	"RHO":[-4.7, -1.0],
	"PT":[450., 1000.],
	"DCSV":0.85,
	"N2DDT":0.,
	"PT_BINS":[450., 500.,550.,600.,675.,800.,1000.],
	"BB_SF":0.91,
	"BB_SF_ERR":0.03,
	"V_SF":0.993,
	"V_SF_ERR":0.043,
}

# Signal sample bookkeeping
signal_names = [] # All signal names (simulated_signal_names + interpolated_signal_names)
simulated_signal_names = []
signal_models = ["Sbb", "PSbb"]
signal_model_masses = [50,100,125,200,300,350,400,500]  # 400,500
signal_masses = {} # E.g. "Sbb100":100
for model in ["Sbb", "PSbb"]:
	for mass in signal_model_masses:
		this_signal_name = "{}{}".format(model, mass)
		signal_names.append(this_signal_name)
		simulated_signal_names.append(this_signal_name)
		signal_masses[this_signal_name] = mass

interpolated_signal_masses = [x for x in range(50, 525, 25) if not x in signal_model_masses]
interpolated_signal_names = []
for model in ["Sbb", "PSbb"]:
	for mass in interpolated_signal_masses:
		this_signal_name = "{}{}".format(model, mass)
		signal_names.append(this_signal_name)
		interpolated_signal_names.append(this_signal_name)
		signal_masses[this_signal_name] = mass

# Masses to be used in limit setting. AK8 can't go above 400 GeV because of no statistics.
limit_signal_masses = {
	"AK8":range(50,425,25),
	"CA15":range(50,525,25)
}
limit_signal_names = {"AK8":[], "CA15":[]}
for jet_type in ["AK8", "CA15"]:
	for model in ["Sbb", "PSbb"]:
		for mass in limit_signal_masses[jet_type]:
			limit_signal_names[jet_type].append("{}{}".format(model, mass))

# File locations
paths = {
	"data":os.path.expandvars("$HOME/DAZSLE/data/"),
	"LimitSetting":os.path.expandvars("$HOME/DAZSLE/data/LimitSetting")
}

def get_histogram_file(selection, jet_type):
	return paths["LimitSetting"] + "/Xbb_inputs/histograms_{}_{}.root".format(selection, jet_type)

def get_interpolation_file(jet_type):
	return paths["LimitSetting"] + "/Xbb_inputs/interpolations_{}.root".format(jet_type)

def get_datacard_directory(signal_name, jet_type, qcd=False):
	if qcd:
		return paths["LimitSetting"] + "/Xbb_inputs/SR_{}/card_qcd_mcstat/{}".format(jet_type, signal_name)
	else:
		return paths["LimitSetting"] + "/Xbb_inputs/SR_{}/card_mcstat/{}".format(jet_type, signal_name)
