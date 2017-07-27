analysis_parameters = {}
analysis_parameters["AK8"] = {
	"MASS_BINS":52,
	"MSD":[40, 404],
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
	"MASS_BINS":52,
	"MSD":[40, 404],
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
	interpolated_signal_masses[model] = 
	for mass in interpolated_signal_masses:
		this_signal_name = "{}{}".format(model, mass)
		signal_names.append(this_signal_name)
		interpolated_signal_names.append(this_signal_name)
		signal_masses[this_signal_name] = mass
