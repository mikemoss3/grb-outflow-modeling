import numpy as np

def write_jet_params( tw = 20,
	dte = 0.001,
	E_dot_iso = 1e52,
	theta = 0.1,
	r_open = 1e6,
	eps_th = 0.03,
	sigma = 0.1,
	eps_e_int = 0.3,
	eps_b_int = 0.3,
	zeta_int = 0.001,
	p_int = 2.2,
	eps_e_ext = 0.1,
	eps_b_ext = 1e-4,
	zeta_ext = 1,
	p_ext = 2.2,
	k = 0,
	rho_not = 1.672e-24 ):
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

	fn = "jet-params.txt"

	return 0;	

# LorentzDist methods
def write_lor_dis_fred_inject(num_freds=2,fred_start_arr=np.array([-0.6,4]),tau_1_arr=np.array([1,1.5]),tau_2_arr=np.array([1.5,3]),fred_amp_arr=np.array([450,100]),median=0,decay=0,fluctuate=True):
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
	fluctuations:	bool, 		Indicates whether to include fluctuations in the distribution or not
	"""

	fn = "jet-params.txt"


	return 0;

def write_lor_dis_gauss_inject(num_gauss=2,gauss_mean_arr=np.array([1, 3]),gauss_amp_arr=np.array([300,200]),gauss_sig_arr=np.array([2, 1.5]),num_sig_arr=np.array([4,4]),median=0,decay=0,fluctuate=True):
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
	fluctuations:	bool, 		Indicates whether to include fluctuations in the distribution or not
	"""

	return 0;

def write_lor_dis_linear(fn,gamma_start=500., gamma_end=2.,fluctuate=True):
	"""
	Method to create a Lorentz distribution from a simple line extending from the starting to ending Lorentz factor.

	Attributes:
	(variable)		(type)	(description)
	gamma_start:	float,	Starting Lorentz factor
	gamma_end:		float,	Ending Lorentz factor

	"""
	return 0;

def write_lor_dis_osci(median=433.,amp=0.66,freq=5,decay=0.5,fluctuate=True):
	"""
	Method to create a Lorentz distribution from an sinusoidal function with an exponential decay term (see Equation 16 of Hascoet, Daigne, and Mochkovitch 2013).
	
	Attributes:
	(variable)		(type)	(description)
	median:			float,	Offset from zero
	amp:			float,	Amplitude of sinusoidal
	freq:			float,	Frequency of sinusoidal function
	decay:			float,	Exponential decay index
	fluctuations:	bool, 	Indicates whether to include fluctuations in the distribution or not
	"""
	return 0;

def write_lor_dis_smoothstep(gamma_1 = 100, gamma_2 = 400, tfrac = 0.5, fluctuate = True):
	"""
	Method to create a Lorentz distribution from an smoothed step function.

	Attributes:
	(variable)		(type)	(description)
	g1:				float,	Lorentz factor before transition
	g2:				float,	Lorentz factor after transition
	tfrac:			float,	Fraction of the entire burst duration at which the Lorentz factor transition begins
	fluctuations:	bool, 	Indicates whether to include fluctuations in the distribution or not
	"""
	return 0;

def write_lor_dis_square_inject(num_squares = 3, square_starts = np.array([0, 5 , 15]), square_durations = np.array([5, 5, 5]), square_amps = np.array([200, 300,50 ]),median=0,decay=0, fluctuate=True):
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
	fluctuations:		bool, 		Indicates whether to include fluctuations in the distribution or not 
	"""
	return 0;

def write_lor_dis_step(gamma_1=100,gamma_2=400,mfrac=0.5,fluctuate=True):
	"""
	Method to create a Lorentz distribution from a step function.

	Attributes:
	(variable)		(type)	(description)
	g1:				float,	Lorentz factor before transition
	g2:				float,	Lorentz factor after transition
	mfrac:			float,	Fraction of outflow mass at which the Lorentz factor transition begins
	fluctuations:	bool, 	Indicates whether to include fluctuations in the distribution or not 
	"""

	return 0;