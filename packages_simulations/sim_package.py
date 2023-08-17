import numpy as np
import subprocess
from subprocess import STDOUT


def run_main(
	save_dir = "./files-data/synthetic-data/",
	param_file_name = './packages_simulations/files-input/jet-params.txt',
	write_lor_prof = False,
	track_lor_prof_evo = False,
	write_jet_params = False,
	light_curve = True,
	light_curve_time_int=np.array([-1e10,-1e9]),
	light_curve_energy_int=np.array([-1e10,-1e9]),
	time_log_scale = False,
	time_step = 0.1,
	spectra = True,
	spectra_time_ints=np.array([(-1e10,-1e9)]),
	spectra_energy_int=np.array([-1e10,-1e9]),
	num_energy_bins = 400,
	num_spectra = 0,
	inc_comps = True
	):
	"""
	Method to run main simulation code.

	Attributes:
	(variable)				(type)		(unit)	(description)
	save_dir:				string,				indicates which directory to save everything (default is ../files-data/synthetic-data)
	param_file_name: 		string,				indicates the jet parameter file (default is jet-params.txt)
	write_lor_prof:			bool,				indicates whether to or not write the Lorentz distribution to a file
	track_lor_prof_evo:		bool,				indicates whether or not to write snapshots of the Lorentz distribution at different times
	write_jet_params:		bool,				indicates whether or not to write out the jet parameters calculated during the simulation to a file
	
	light_curve:			bool,				indicates whether to calculate a light curve
	light_curve_time_int:	np.ndarray,	sec		time interval to calculate the light curve over, indicated by a start and stop time 
	light_curve_energy_int:	np.ndarray,	keV		energy interval to calculate the light curve over, indicated by a minimum and maximum energy
	time_log_scale:			bool,				indicates whether the time axis should be generated in linear space (False) or log space (True)
	time_step:				float,				the time step along the time axis

	spectra:				bool,				indicates whether to calculate spectra
	spectra_time_ints:		np.ndarray,	sec		time intervals to calculate the spectra over, indicated by an array of tuples with start and stop times 
	spectra_energy_int:		np.ndarray,	keV		energy interval to calculate the spectra over, indicated by a minimum and maximum energy
	num_energy_bins:		int,				number of energy bins along the energy axis
	num_spectra:			int,				number of spectra to be created
	
	inc_comps:				bool,				indicates whether to record the separate emission components that compose the generated spectra or light curves
	"""


	command_string = "./packages_simulations/main {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17}".format(
		save_dir,
		param_file_name,
		int(write_lor_prof),
		int(track_lor_prof_evo),
		int(write_jet_params),
		int(light_curve),
		light_curve_time_int[0],
		light_curve_time_int[1],
		light_curve_energy_int[0],
		light_curve_energy_int[1],
		int(time_log_scale),
		time_step,
		int(spectra),
		spectra_energy_int[0],
		spectra_energy_int[1],
		num_energy_bins,
		num_spectra,
		int(inc_comps)
		)
		
	for i in range(len(spectra_time_ints)):
		command_string+=" {} {}".format(spectra_time_ints[i][0], spectra_time_ints[i][1])


	subprocess.run([command_string],shell=True,stderr=STDOUT)

	return 0;