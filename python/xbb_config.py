import os

analysis_parameters = {}
analysis_parameters["AK8"] = {
	"MASS_BINS":80,
	"MSD":[40, 600],
	"RHO":[-6.0, -2.1],
	"PT":[450., 1000.],
	"DCSV":0.9,
	"DCSV_LOOSE":0.8,
	"N2DDT":0.,
	"PT_BINS":[450., 500.,550.,600.,675.,800.,1000.],
	"BB_SF":0.91,
	"BB_SF_ERR":0.03,
	"V_SF":0.993,
	"V_SF_ERR":0.043,
	"MAX_NRHO":3,
	"MAX_NPT":2,
	"DEFAULT_NRHO":2,
	"DEFAULT_NPT":1,
	"VFAIL_SF":1.025,
	"VFAIL_SF_ERR":0.043,
}
analysis_parameters["CA15"] = {
	"MASS_BINS":76,
	"MSD":[68, 600],
	"RHO":[-4.7, -1.0],
	"PT":[450., 1000.],
	"DCSV":0.9,
	"DCSV_LOOSE":0.8,
	"N2DDT":0.,
	"PT_BINS":[450., 500.,550.,600.,675.,800.,1000.],
	"BB_SF":1.0,
	"BB_SF_ERR":0.04,
	"V_SF":0.968,
	"V_SF_ERR":0.058,
	"MAX_NRHO":6,
	"MAX_NPT":2,
	"DEFAULT_NRHO":5,
	"DEFAULT_NPT":1,
	"VFAIL_SF":1.0628742515,
	"VFAIL_SF_ERR":0.058,
}

# Signal sample bookkeeping
signal_names = [] # All signal names (simulated_signal_names + interpolated_signal_names)
simulated_signal_names = []
signal_models = ["Sbb", "PSbb", "ZPrime"]
signal_model_masses = {
	"Sbb":[50,100,125,200,300,350,400,500], 
	"PSbb":[50,100,125,200,300,350,400,500], 
	"ZPrime":[50, 75, 100, 125, 200, 300],
}
signal_masses = {} # E.g. "Sbb100":100
for model in signal_models:
	for mass in signal_model_masses[model]:
		this_signal_name = "{}{}".format(model, mass)
		signal_names.append(this_signal_name)
		simulated_signal_names.append(this_signal_name)
		signal_masses[this_signal_name] = mass

interpolated_signal_masses = {}
interpolated_signal_names = []
for model in signal_models:
	if model == "ZPrime":
		interpolated_signal_masses[model] = [x for x in range(50, 325, 25) if not x in signal_model_masses[model]]
	else:
		interpolated_signal_masses[model] = [x for x in range(50, 525, 25) if not x in signal_model_masses[model]]
	for mass in interpolated_signal_masses[model]:
		this_signal_name = "{}{}".format(model, mass)
		signal_names.append(this_signal_name)
		interpolated_signal_names.append(this_signal_name)
		signal_masses[this_signal_name] = mass
#for mass in range(50, 325, 25):
#	this_signal_name = "{}{}".format("ZPrime", mass)
#	signal_names.append(this_signal_name)
#	interpolated_signal_names.append(this_signal_name)
#	signal_masses[this_signal_name] = mass

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
	"LimitSetting":os.path.expandvars("$HOME/DAZSLE/data/LimitSetting/"),
	"Fits":os.path.expandvars("$HOME/DAZSLE/data/Fits/"),
}

def get_histogram_file(selection, jet_type):
	return paths["LimitSetting"] + "/Xbb_inputs/histograms_{}_{}.root".format(selection, jet_type)

def get_interpolation_file(selection, jet_type):
	return paths["LimitSetting"] + "/Xbb_inputs/interpolations_{}_{}.root".format(selection, jet_type)

def get_datacard_directory(signal_name, jet_type, qcd=False, decidata=False, region="SR"):
	if qcd:
		return paths["LimitSetting"] + "Xbb_inputs/{}_{}/cards_qcd_mcstat/{}".format(region, jet_type, signal_name)
	elif decidata:
		return paths["LimitSetting"] + "Xbb_inputs/{}_{}/cards_ps10_mcstat/{}".format(region, jet_type, signal_name)
	else:
		return paths["LimitSetting"] + "Xbb_inputs/{}_{}/cards_mcstat/{}".format(region, jet_type, signal_name)

def get_ftest_directory(signal_name, jet_type, nrho1, npt1, nrho2, npt2, qcd=False, decidata=False, region="SR", ):
	category = "{}_{}".format(region, jet_type)
	if qcd:
		category += "_pseudodata"
	elif decidata:
		category += "_ps10"
	return paths["Fits"] + "/ftest/{}/r{}p{}_vs_r{}p{}/{}".format(category, nrho1, npt1, nrho2, npt2, signal_name)