import numpy as np
import matplotlib.pyplot as plt 
import scipy.integrate as integrate 


c = 3*10**10
a = 7.566 * np.power(10.,-15.); 
kb_kev = 8.617*np.power(10.,-8. )

# Define constants that are used between methods
omega_m = 0.3 # matter density of the universe
omega_lam = 0.7 # dark energy density of the universe
H0 = 67.4*np.power(10,5) # Hubbles Constant cm/s/Mpc

def beta(gamma):
	return 1 - 1./2./np.power(gamma,2.)

def make_lor_prof(t_w,dt):

	num_shells = int(t_w / dt)
	lor_prof = np.zeros(shape = num_shells,dtype=[("TI",float),("GAMMA",float),("RADIUS",float)])


	lor_prof['TI'] = dt*np.linspace(start=0,stop = num_shells,num=num_shells)
	lor_prof['GAMMA'] = 333*(1+(2/3)*np.cos(5*np.pi*(1-(lor_prof['TI']/t_w))))*np.exp(- lor_prof['TI']/2/t_w)
	lor_prof['RADIUS'] =  - lor_prof['TI'] * beta(lor_prof['GAMMA'])

	return lor_prof

def calc_therm_params(lor_prof,E_iso,theta,ell,t_w,eps_therm,sigma):

	therm_lum = np.zeros(shape = len(lor_prof),dtype=[("TE",float),("TA",float),("DT",float),("RADIUS",float),("TEMP",float),("LUM",float)])

	E_iso_dot = E_iso / t_w
	E_dot = E_iso_dot * np.power(theta,2.)/4.
	E_dot_therm = E_dot * eps_therm

	therm_lum['RADIUS'] = 2.9 * 1.e13 * (E_iso_dot/1.e53) / (1.+sigma) /  np.power(lor_prof['GAMMA']/100.,3.) / c
	# therm_lum['RADIUS'] = 0.2 * (E_dot/lor_prof['GAMMA']/c**2) / (8*np.pi * c * lor_prof['GAMMA']**2) / c
	
	therm_lum['TE'] = (therm_lum['RADIUS'] - lor_prof['RADIUS']) / beta(lor_prof['GAMMA'])
	therm_lum['TA'] = therm_lum['TE'] - therm_lum['RADIUS']
	therm_lum['DT'] = therm_lum['RADIUS'] / 2 / np.power(lor_prof['GAMMA'],2)

	Phi = np.power(theta,-2/3) * np.power(therm_lum['RADIUS']*c,-2/3) * np.power(ell,2/3) * np.power(lor_prof['GAMMA'],2/3)

	T0 = np.power(E_dot_therm / c / a / (np.pi * ell**2), 0.25)
	therm_lum['TEMP'] = T0 * Phi

	# therm_lum['LUM'] = E_dot_therm * Phi
	# therm_lum['LUM'] = np.power(lor_prof['GAMMA'],8/3) * eps_therm * c**2 * np.power(E_dot,1/3) * np.power(8*np.pi*ell/0.2/theta,2/3) # Using M_dot = E_dot / Gamma c^2
	therm_lum['LUM'] = np.power(lor_prof['GAMMA'],8/3) * eps_therm * c**2 * E_dot * np.power(E_iso_dot,-2/3) * np.power(8*np.pi*ell/0.2/theta,2/3) # Using M_dot = E_iso_dot / Gamma c^2

	return therm_lum

def calc_therm_spec(therm_spectrum,therm_lum,tmin,tmax):
	
	spectrum_sum = 0
	therm_spectrum['RATE'] *= 0

	# Calculate the normalization
	norm = 0.; 
	e_hi = 5e3;
	e_lo = 1e-6;
	norm_num_bin= int(20.*np.log10(e_hi/e_lo));
	norm_energy_axis = np.logspace(start=np.log10(e_lo),stop=np.log10(e_hi),num=norm_num_bin);
	norm_de_axis = np.zeros(shape=len(norm_energy_axis))
	norm_de_axis[0:-1] = [(norm_energy_axis[i+1] - norm_energy_axis[i]) for i in range(len(norm_energy_axis)-1)]
	norm_de_axis[-1] = therm_spectrum['DE'][-2]

	for i in range(len(therm_lum)):
		if (therm_lum['TA'][i] <= tmax) and ((therm_lum['TA'][i]+therm_lum['DT'][i] >= tmin) ):
			
			en_curr=0.
			norm=0.
			for j in range(norm_num_bin-1):
				en_curr = norm_energy_axis[j];
				if(en_curr < 5.*np.power(10.,3.)):
					norm += norm_de_axis[j] * ThermalSpec(en_curr,therm_lum['TEMP'][i]);

			tmp_val=0.;
			for j in range(len(therm_spectrum)):
				en_curr = therm_spectrum['ENERGY'][j];
					
				if(en_curr < 5.*pow(10.,3.)):
					tmp_val = therm_lum['LUM'][i] * ThermalSpec(en_curr,therm_lum['TEMP'][i]) / norm;
				else:
					tmp_val=0;

				# Add the contribution to the total spectrum according to a Left-Riemann-Sum
				
				spectrum_sum += therm_spectrum['DE'][j] * tmp_val;
				# Check if the spectrum rate per energy bin is requested and store it if so. 
				therm_spectrum['RATE'][j] += tmp_val;

	return spectrum_sum

def ThermalSpec(energy,temp, alpha=0.4,norm=1):

	return norm*np.power(energy/temp/kb_kev,1.+alpha)/(np.exp(energy/temp/kb_kev) - 1.);


def calc_therm_light_curve(therm_lum,emin,emax):

	light_curve = np.zeros(shape=len(therm_lum),dtype=[("TIME",float),("RATE",float)])
	light_curve['TIME'] = therm_lum['TA']

	spec_e_min = emin
	spec_e_max = emax
	num_e_bins = 100

	therm_spectrum = np.zeros(shape = num_e_bins, dtype=[("ENERGY",float),("RATE",float),("DE",float)])
	therm_spectrum['ENERGY'] = np.logspace(start=np.log10(spec_e_min), stop=np.log10(spec_e_max),num=num_e_bins)
	therm_spectrum['DE'][0:-1] = [(therm_spectrum['ENERGY'][i+1] - therm_spectrum['ENERGY'][i]) for i in range(len(therm_spectrum['ENERGY'])-1)]
	therm_spectrum['DE'][-1] = therm_spectrum['DE'][-2]

	e_hi = 5e3;
	e_lo = 1e-6;
	norm_num_bin= int(20.*np.log10(e_hi/e_lo));
	norm_energy_axis = np.logspace(start=np.log10(e_lo),stop=np.log10(e_hi),num=norm_num_bin);
	norm_de_axis = np.zeros(shape=len(norm_energy_axis))
	norm_de_axis[0:-1] = [(norm_energy_axis[i+1] - norm_energy_axis[i]) for i in range(len(norm_energy_axis)-1)]
	norm_de_axis[-1] = therm_spectrum['DE'][-2]

	for i in range(len(light_curve)-1):
		light_curve['RATE'][i] = calc_therm_spec(therm_spectrum=therm_spectrum, therm_lum=therm_lum, tmin=light_curve['TIME'][i], tmax=light_curve['TIME'][i+1])/(light_curve['TIME'][i+1]-light_curve['TIME'][i])

	return light_curve

def lum_dis(z: float):
	""" 
	Caclulate luminosity distance for redshift z
	"""
	lum_dis_Mpc = ((1+z)*c/H0 ) * integrate.quad(lambda zi: 1/np.sqrt( ((omega_m*np.power(1+zi,3) )+omega_lam) ),0,z)[0]
	lum_dis_cm = lum_dis_Mpc * 3.086e24 # Mpc -> cm
	return lum_dis_cm

if __name__ == "__main__":

	t_w = 10
	dt = 0.01
	E_iso = 1e54
	theta = 0.1
	ell = 3e6
	eps_therm = 0.03
	sigma = 0.1

	z = 1

	lor_prof = make_lor_prof(t_w = t_w, dt = dt)

	# ax = plt.figure().gca()
	# ax.scatter(-lor_prof['RADIUS'],lor_prof['GAMMA']/100)

	therm_lum = calc_therm_params(lor_prof=lor_prof, E_iso=E_iso, theta=theta, ell=ell, t_w=t_w, eps_therm=eps_therm, sigma=sigma)

	# ax = plt.figure().gca()
	# ax.scatter(therm_lum['TA']*(1+z),therm_lum['LUM'])
	# ax = plt.figure().gca()
	# ax.scatter(therm_lum['TA']*(1+z),therm_lum['TEMP']*kb_kev/(1+z))
	# ax.set_yscale("log")
	# ax = plt.figure().gca()
	# ax.scatter(therm_lum['TA'],lor_prof['GAMMA']**(8./3.))

	spec_e_min = 0.1
	spec_e_max = 1e4
	num_e_bins = 100
	therm_spectrum = np.zeros(shape = num_e_bins, dtype=[("ENERGY",float),("RATE",float),("DE",float)])
	therm_spectrum['ENERGY'] = np.logspace(start=np.log10(spec_e_min), stop=np.log10(spec_e_max),num=num_e_bins)
	therm_spectrum['DE'][0:-1] = [(therm_spectrum['ENERGY'][i+1] - therm_spectrum['ENERGY'][i]) for i in range(len(therm_spectrum['ENERGY'])-1)]
	therm_spectrum['DE'][-1] = therm_spectrum['DE'][-2]
	tmin = 0
	tmax = 30

	ax = plt.figure().gca()
	calc_therm_spec(therm_spectrum=therm_spectrum, therm_lum=therm_lum,tmin= tmin,tmax= tmax)
	therm_spectrum['RATE']/=(tmax-tmin)
	ax.scatter(therm_spectrum['ENERGY']/(1+z),therm_spectrum['RATE']*(therm_spectrum['ENERGY']**2),label="Time Int")

	# t_arr = np.array([(3.9,4.1),(4.6,4.8),(11.8,12.)])
	# for i in range(len(t_arr)):
	# 	therm_spectrum['RATE'] *= 0.
	# 	tmin = t_arr[i][0] / (1+z)
	# 	tmax = t_arr[i][1] / (1+z)
	# 	calc_therm_spec(therm_spectrum=therm_spectrum, therm_lum=therm_lum,tmin= tmin,tmax= tmax)
	# 	therm_spectrum['RATE']/=(tmax-tmin)
	# 	ax.scatter(therm_spectrum['ENERGY']/(1+z),therm_spectrum['RATE']*(therm_spectrum['ENERGY']**2),label="{:.2f} - {:.2f} s".format(t_arr[i][0],t_arr[i][1]))

	ax.set_xscale("log")
	ax.set_yscale("log")

	ax.set_ylim(1e48,1e52)
	ax.set_xlim(0.1,1e5)
	ax.legend()

	# plt.savefig("figs/2023-08-23-spectra-time-res.png")
	
	"""
	lc_t_min = 0
	lc_t_max = 10
	dt = 0.01
	lc_e_min = 8
	lc_e_max = 4e4

	light_curve = calc_therm_light_curve(therm_lum=therm_lum,emin=lc_e_min,emax=lc_e_max)

	light_curve['TIME'] *= (1+z)
	light_curve['RATE'] /= (4*np.pi*lum_dis(z)**2)

	ax = plt.figure().gca()
	# colors = ["C1","C2","C3"]
	# for i in range(len(t_arr)):
	# 	ax.axvspan(xmin=t_arr[i][0],xmax=t_arr[i][1],alpha=0.5,color=colors[i])

	ax.scatter(light_curve['TIME'],light_curve['RATE'],marker=".")

	# plt.savefig("figs/2023-08-23-light-curves-v2.png")
	"""

	plt.show()

"""
Some notes:
1. In the equation for R_ph (Eq 9 in H13), is M_dot the isotropic equivalent? 
(because I get better agreement when I use the isotropic equivalent power
than the beaming corrected power)


"""