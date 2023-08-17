import numpy as np

def write_jet_params( tw = 20,dte = 0.001,E_dot_iso = 1e52,
	theta = 0.1,r_open = 1e6,eps_th = 0.03,sigma = 0.1,
	eps_e_int = 0.3,eps_b_int = 0.3,zeta_int = 0.001,p_int = 2.2,
	eps_e_ext = 0.1,eps_b_ext = 1e-4,zeta_ext = 1,p_ext = 2.2,
	k = 0,rho_not = 1.672e-24 ):
	"""
	Method to write the jet parameter file that is used as an input file for the GRB outflow simulation code.

	Attributes:
	(variable)	(type)			(description)
	tw: 		sec,			duration of the wind
	dte: 		sec,			time between subsequent shell launches
	E_dot_iso: 	erg/s,	 		isotropic equivalent energy injection rate
	theta: 		rad,			jet opening angle
	r_open: 	cm,		 		jet opening radius
	eps_th: 	fraction,		the fraction of jet energy stored as thermal energy
	sigma: 		dimensionless,	the jet magnetization
	eps_e_int: 	fraction,		fraction of jet energy stored in the accelerated electron population (for material within the jet)
	eps_b_int: 	fraction,		fraction of jet energy stored in the magnetic field (for material within the jet)
	zeta_int: 	fraction,		fraction of electrons accelerated by internal shocks
	p_int: 		dimensionless,	power law index of the accelerated electron population distribution (for material within the jet)
	eps_e_ext: 	fraction,		fraction of jet energy stored in the accelerated electron population (for material external to the jet)
	eps_b_ext: 	fraction,		fraction of jet energy stored in the magnetic field (for material external to the jet)
	zeta_ext: 	fraction,		fraction of electrons accelerated by forward shocks
	p_ext: 		dimensionless,	power law index of the accelerated electron population distribution (for material external to the jet)
	k: 			dimensionless,	external medium density profile power law index, i.e., constant = 0, wind = 2 ###
	rho_not: 					density normalization (for k = 0, g cm^-3, i.e., n0*mp = n0*1.672e-24 | for k = 2, g cm^-1, i.e.,A_star * 5.e11 )
	"""

	fn = "./packages_simulations/files-input/jet-params.txt"
	file = open(fn,"w")

	file.write("{}\n".format(tw))
	file.write("{}\n".format(dte))
	file.write("{}\n".format(E_dot_iso))
	file.write("{}\n".format(theta))
	file.write("{}\n".format(r_open))
	file.write("{}\n".format(eps_th))
	file.write("{}\n".format(sigma))
	file.write("{}\n".format(eps_e_int))
	file.write("{}\n".format(eps_b_int))
	file.write("{}\n".format(zeta_int))
	file.write("{}\n".format(p_int))
	file.write("{}\n".format(eps_e_ext))
	file.write("{}\n".format(eps_b_ext))
	file.write("{}\n".format(zeta_ext))
	file.write("{}\n".format(p_ext))
	file.write("{}\n".format(k))
	file.write("{}\n".format(rho_not))

	file.close()


	return 0;	

# LorentzDist methods
def write_lor_dis_linear(gamma_start, gamma_end,fluctuate=True):
	"""
	Method to create a Lorentz distribution from a simple line extending from the starting to ending Lorentz factor.

	Attributes:
	(variable)		(type)	(description)
	gamma_start:	float,	Starting Lorentz factor
	gamma_end:		float,	Ending Lorentz factor
	fluctuate:		bool, 		Indicates whether to include fluctuations in the distribution or not
	"""
	fn_jet_params = "./packages_simulations/files-input/jet-params.txt"
	fn_lor_prof = "./packages_simulations/files-input/lor-prof-params.txt"

	file_jet_params = open(fn_jet_params,"a")
	file_jet_params.write("linear")
	file_jet_params.write(fn_lor_prof)
	file_jet_params.close()

	file_lor_prof = open(fn_lor_prof,"w")
	file_lor_prof.write("{}".format(gamma_start))
	file_lor_prof.write("{}".format(gamma_end))
	file_lor_prof.write("{}".format(fluctuate))
	file_lor_prof.close()

	return 0;

def write_lor_dis_step(gamma_1,gamma_2,mfrac,fluctuate=True):
	"""
	Method to create a Lorentz distribution from a step function.

	Attributes:
	(variable)		(type)	(description)
	gamma_1:		float,	Lorentz factor before transition
	gamma_2:		float,	Lorentz factor after transition
	mfrac:			float,	Fraction of outflow mass at which the Lorentz factor transition begins
	fluctuate:		bool, 	Indicates whether to include fluctuations in the distribution or not 
	"""

	fn_jet_params = "./packages_simulations/files-input/jet-params.txt"
	fn_lor_prof = "./packages_simulations/files-input/lor-prof-params.txt"

	file_jet_params = open(fn_jet_params,"a")
	file_jet_params.write("osci")
	file_jet_params.write(fn_lor_prof)
	file_jet_params.close()

	file_lor_prof = open(fn_lor_prof,"w")
	file_lor_prof.write("{}".format(gamma_1))
	file_lor_prof.write("{}".format(gamma_2))
	file_lor_prof.write("{}".format(mfrac))
	file_lor_prof.write("{}".format(fluctuate))
	file_lor_prof.close()

	return 0;

def write_lor_dis_smoothstep(gamma_1, gamma_2, tfrac, fluctuate = True):
	"""
	Method to create a Lorentz distribution from an smoothed step function.

	Attributes:
	(variable)		(type)	(description)
	gamma_1:		float,	Lorentz factor before transition
	gamma_2:		float,	Lorentz factor after transition
	tfrac:			float,	Fraction of the entire burst duration at which the Lorentz factor transition begins
	fluctuate:		bool, 	Indicates whether to include fluctuations in the distribution or not
	"""
	fn_jet_params = "./packages_simulations/files-input/jet-params.txt"
	fn_lor_prof = "./packages_simulations/files-input/lor-prof-params.txt"

	file_jet_params = open(fn_jet_params,"a")
	file_jet_params.write("osci")
	file_jet_params.write(fn_lor_prof)
	file_jet_params.close()

	file_lor_prof = open(fn_lor_prof,"w")
	file_lor_prof.write("{}".format(gamma_1))
	file_lor_prof.write("{}".format(gamma_2))
	file_lor_prof.write("{}".format(tfrac))
	file_lor_prof.write("{}".format(fluctuate))
	file_lor_prof.close()

	return 0;

def write_lor_dis_square_inject(num_squares, square_starts, square_durations, square_amps,median=0,decay=0, fluctuate=True):
	"""
	Method to create a Lorentz distribution as a combination of square distributions.

	Attributes:
	(variable)			(type)		(description)
	num_squares:		int,		Number of square distributions to include in the Lorentz distribution
	square_starts:		np.ndarray,	Array of start times for each square
	square_durations:	np.ndarray,	Array of duration times for each square
	square_amps:		np.ndarray,	Array of amplitudes for each square
	median:				float,		An optional flat offset from zero
	decay:				float,		An optional exponential decay to the flat offset
	fluctuate:			bool, 		Indicates whether to include fluctuations in the distribution or not 
	"""
	fn_jet_params = "./packages_simulations/files-input/jet-params.txt"
	fn_lor_prof = "./packages_simulations/files-input/lor-prof-params.txt"

	file_jet_params = open(fn_jet_params,"a")
	file_jet_params.write("osci\n")
	file_jet_params.write(fn_lor_prof+"\n")
	file_jet_params.close()

	file_lor_prof = open(fn_lor_prof,"w")
	file_lor_prof.write("{}\n".format(median))
	file_lor_prof.write("{}\n".format(decay))
	file_lor_prof.write("{}\n".format(num_squares))
	
	for i in range(num_squares):
		file_lor_prof.write("{} ".format(square_starts[i]))
	file_lor_prof.write("\n")

	for i in range(num_squares):
		file_lor_prof.write("{} ".format(square_durations[i]))
	file_lor_prof.write("\n")

	for i in range(num_squares):
		file_lor_prof.write("{} ".format(square_amps[i]))
	file_lor_prof.write("\n")

	file_lor_prof.write("{}\n".format(fluctuate))
	file_lor_prof.close()

	return 0;

def write_lor_dis_fred_inject(num_freds,fred_start_arr,tau_1_arr,tau_2_arr,fred_amp_arr,median=0,decay=0,fluctuate=True):
	"""
	Method to create a Lorentz distribution as a combination of Fast Rise Exponential Decay (FRED) distributions.
	FRED functions based on Hakkila and Preece 2014

	Attributes:
	(variable)		(type)		(description)
	num_freds:		int,		Number of FRED distributions to include in the Lorentz distribution
	fred_start_arr: np.ndarray,	Array indicating the time the FRED begins
	tau_1_arr: 		np.ndarray, Array of FRED rise-time parameters
	tau_2_arr: 		np.ndarray,	Array of FRED decay-time parameters
	fred_amp_arr:	np.ndarray,	Array of FRED amplitudes
	median:			float,		An optional flat offset from zero
	decay:			float,		An optional exponential decay to the flat offset
	fluctuate:		bool, 		Indicates whether to include fluctuations in the distribution or not
	"""

	fn_jet_params = "./packages_simulations/files-input/jet-params.txt"
	fn_lor_prof = "./packages_simulations/files-input/lor-prof-params.txt"

	file_jet_params = open(fn_jet_params,"a")
	file_jet_params.write("fred_inject\n")
	file_jet_params.write(fn_lor_prof+"\n")
	file_jet_params.close()

	file_lor_prof = open(fn_lor_prof,"w")
	file_lor_prof.write("{}\n".format(median))
	file_lor_prof.write("{}\n".format(decay))
	file_lor_prof.write("{}\n".format(num_freds))

	for i in range(num_freds):
		file_lor_prof.write("{} ".format(fred_start_arr[i]))
	file_lor_prof.write("\n")

	for i in range(num_freds):
		file_lor_prof.write("{} ".format(tau_1_arr[i]))
	file_lor_prof.write("\n")

	for i in range(num_freds):
		file_lor_prof.write("{} ".format(tau_2_arr[i]))
	file_lor_prof.write("\n")

	for i in range(num_freds):
		file_lor_prof.write("{} ".format(fred_amp_arr[i]))
	file_lor_prof.write("\n")

	file_lor_prof.write("{}\n".format(fluctuate))
	file_lor_prof.close()

	return 0;

def write_lor_dis_gauss_inject(num_gauss,gauss_mean_arr,gauss_amp_arr,gauss_sig_arr,num_sig_arr,median=0,decay=0,fluctuate=True):
	"""
	Method to create a Lorentz distribution as a combination of Gaussian distributions.

	Attributes:
	(variable)		(type)		(description)
	num_gauss:		int,		Number of Gauss distributions to include in the Lorentz distribution
	gauss_mean_arr:	np.ndarray,	Array of Gauss mean time values (i.e., at what time is the Gaussian centered)
	gauss_amp_arr:	np.ndarray,	Array of Gauss amplitudes values
	gauss_sig_arr:	np.ndarray,	Array of Gauss sigma values
	num_sig_arr:	np.ndarray,	Array of integers indicating how many sigma the Gaussian should go out to.
	median:			float,		An optional flat offset from zero
	decay:			float,		An optional exponential decay to the flat offset
	fluctuate:		bool, 		Indicates whether to include fluctuations in the distribution or not
	"""

	fn_jet_params = "./packages_simulations/files-input/jet-params.txt"
	fn_lor_prof = "./packages_simulations/files-input/lor-prof-params.txt"

	file_jet_params = open(fn_jet_params,"a")
	file_jet_params.write("gauss_inject\n")
	file_jet_params.write(fn_lor_prof+"\n")
	file_jet_params.close()

	file_lor_prof = open(fn_lor_prof,"w")
	file_lor_prof.write("{}\n".format(median))
	file_lor_prof.write("{}\n".format(decay))
	file_lor_prof.write("{}\n".format(num_gauss))

	for i in range(num_gauss):
		file_lor_prof.write("{} ".format(gauss_mean_arr[i]))
	file_lor_prof.write("\n")

	for i in range(num_gauss):
		file_lor_prof.write("{} ".format(gauss_amp_arr[i]))
	file_lor_prof.write("\n")

	for i in range(num_gauss):
		file_lor_prof.write("{} ".format(gauss_sig_arr[i]))
	file_lor_prof.write("\n")

	for i in range(num_gauss):
		file_lor_prof.write("{} ".format(num_sig_arr[i]))
	file_lor_prof.write("\n")

	file_lor_prof.write("{}\n".format(fluctuate))
	file_lor_prof.close()

	return 0;

def write_lor_dis_osci(median,amp,freq,decay,fluctuate=True):
	"""
	Method to create a Lorentz distribution from an sinusoidal function with an exponential decay term (see Equation 16 of Hascoet, Daigne, and Mochkovitch 2013).
	
	Parameters from Hascoet et al 2013: median=433.,amp=0.66,freq=5,decay=0.5

	Attributes:
	(variable)		(type)	(description)
	median:			float,	Offset from zero
	amp:			float,	Amplitude of sinusoidal
	freq:			float,	Frequency of sinusoidal function
	decay:			float,	Exponential decay index
	fluctuate:		bool, 	Indicates whether to include fluctuations in the distribution or not
	"""
	fn_jet_params = "./packages_simulations/files-input/jet-params.txt"
	fn_lor_prof = "./packages_simulations/files-input/lor-prof-params.txt"

	file_jet_params = open(fn_jet_params,"a")
	file_jet_params.write("osci\n")
	file_jet_params.write(fn_lor_prof+"\n")
	file_jet_params.close()

	file_lor_prof = open(fn_lor_prof,"w")
	file_lor_prof.write("{}\n".format(median))
	file_lor_prof.write("{}\n".format(amp))
	file_lor_prof.write("{}\n".format(freq))
	file_lor_prof.write("{}\n".format(decay))
	file_lor_prof.write("{}\n".format(fluctuate))
	file_lor_prof.close()

	return 0;
