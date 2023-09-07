from packages_simulations.file_write_package import *
from packages_simulations.sim_package import run_main

# Set jet parameters
tw = 		10. 	# sec, duration of the wind
dte = 		0.01 	# sec, time between subsequent shell launches
E_dot_iso = 1.e53 	# erg/s, isotropic equivalent energy injection rate
theta = 	0.07 	# rad, jet opening angle
r_open = 	3.e6 	# cm, jet opening radius
eps_th = 	0.02 	# fraction, the fraction of jet energy stored as thermal energy
sigma = 	0.1 	# dimensionless, the jet magnetization
eps_e_int = 0.3 	# fraction, fraction of jet energy stored in the accelerated electron population (for material within the jet)
eps_b_int = 0.07 	# fraction, fraction of jet energy stored in the magnetic field (for material within the jet)
zeta_int = 	1.e-3 	# fraction, fraction of electrons accelerated by internal shocks
p_int = 	2.2 	# dimensionless, power law index of the accelerated electron population distribution (for material within the jet)
eps_e_ext = 0.1 	# fraction, fraction of jet energy stored in the accelerated electron population (for material external to the jet)
eps_b_ext = 1.e-4 	# fraction, fraction of jet energy stored in the magnetic field (for material external to the jet)
zeta_ext = 	1 		# fraction, fraction of electrons accelerated by forward shocks
p_ext = 	2.2 	# dimensionless, power law index of the accelerated electron population distribution (for material external to the jet)
k = 		0 		# dimensionless, external medium density profile power law index, i.e., constant = 0, wind = 2 ###
rho_not = 	1.672e-24 # density normalization (for k = 0, g cm^-3, i.e., n0*mp = n0*1.672e-24 | for k = 2, g cm^-1, i.e.,A_star * 5.e11 )

# Write jet parameters
write_jet_params(tw, dte, E_dot_iso,
	theta, r_open, eps_th, sigma,
	eps_e_int, eps_b_int, zeta_int, p_int,
	eps_e_ext, eps_b_ext, zeta_ext, p_ext,
	k, rho_not)

# Write Lorentz distribution parameters
write_lor_dis_fred_inject(num_freds=2,fred_start_arr=np.array([-0.6,5]),tau_1_arr=np.array([1,2]),tau_2_arr=np.array([3,4]),fred_amp_arr=np.array([450,350]),median=0,decay=0,fluctuate=True)


# Run simulation code
run_main(
	save_dir = "./files-data/synthetic-data/2023-09-07/",
	param_file_name = './packages_simulations/files-input/jet-params.txt',
	write_lor_prof = True,
	track_lor_prof_evo = False,
	write_jet_params = False,
	light_curve = True	,
	light_curve_time_int=np.array([0,10]),
	light_curve_energy_int=np.array([8,4e4]),
	time_log_scale = False,
	time_step = 0.1,
	spectra = True,
	num_spectra = 2,
	spectra_time_ints=np.array([(1,2),(3,4)]),
	spectra_energy_int=np.array([8,4e4]),
	num_energy_bins = 400,
	inc_comps = True
	)