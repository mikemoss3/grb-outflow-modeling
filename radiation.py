"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-03-16


This code defines an Emission class which contains all the emission produced during the prompt jet simulation

This script contains the various spectral functions used to produced the spectra for the burst.
Additionally, this script has methods for creating plots of the spectra and light curve (when supplied with a spectrum)


"""

import numpy as np
import cosmologicalconstants as cc
import utils as ut
from scipy.integrate import quad

# from scipy.integrate import quad


class Emission(object):
	"""
	Spectrum class.
	"""

	def __init__(self,components=None,Tmin=None,Tmax=None,Emin=0.1,Emax=5e5,dE=70):
		"""
		Defines the default parameters of a spectrum.
		"""
		if components is not None:
			self.components = components
		else:
			self.components = []
		# if len(self.components) > 0:
			# self.make_spec()[0]

	def load_spec(self,file_name):
		"""
		Method to load a synchrotron spectrum from a file, where the file path is specified by file_name. These files are created from the promptjet.py simulation code.
		A redshift is given by the user so the spectral values can be calculated in the observer frame (the simulation code calculates in the source frame). 
		"""
		
		# Grab the length of the first line in the file
		first_line_len=len(open(file_name).readline().rstrip().split())

		# Check which type of component this is. We do this by checking the length of the first line
		if first_line_len == 11:
			# Load spectrum component and add it into the components list
			spec_synch = np.genfromtxt(file_name,dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float),('gammar',float),('e',float),('delt',float),('tau',float),('relvel',float)])
			self.components.append(spec_synch)

			# Return the spectral component 
			return spec_synch

		elif first_line_len == 6:
			# Load spectrum component and add it into the components list
			spec_therm = np.genfromtxt(file_name,dtype=[('te',float),('ta',float),('delt',float),('T',float),('Flux',float),('Rphot',float)])
			self.components.append(spec_therm)

			# Return the spectral component 
			return spec_therm

	def make_spec(self,comp_num=None,Tmin=None,Tmax=None,Emin=0.1,Emax=5e5,num_bins=None,z=0,sum_only=False,nuFnu=False):
		"""
		Method to create spectra from the specified component(s). The spectrum is calculated in the observer frame.

		Attributes:
		comp_num=None, User given integer to specify the component to use, if None is given all components will be summed
		
		Tmin=None, Beginning of the time interval over which the spectra will be computed
		Tmax=None, End of the time interval over which the spectra will be computed
		Emin=0.1, Beginning of the energy range over which the spectra will be computed
		Emax=5e5, End of the energy range over which the spectra will be computed
		num_bins=None, Number of bins to use in the summation. By default the number of energy bins will be close to 10 bins per order of magnitude.
		
		z=0, Redshift of the observed source
		sum_only=False, Set this parameter to True if the user only wants the summation of the spectrum across the given time interval and energy range 
			(For instance, this is used to reduce the number of calculation when creating a light curve)
		nuFnu=False, reports the nuFnu spectrum instead of the photon spectrum
		"""

		# Find start and stop time if a time interval is not specified
		# Both start and stop times should be corrected for time dilation
		# The start time will be the first emission from any of the components
		if Tmin is None:
			Tmin=1e99
			for comp in self.components:
				if Tmin > comp['ta'][0]*(1+z):
					Tmin = comp['ta'][0]*(1+z)
		# The stop time is the final (emission start + emission duration)
		if Tmax is None:
			Tmax=-1e99
			for comp in self.components:
				if Tmax < np.max(comp['ta']+comp['delt'])*(1+z):
					Tmax = np.max(comp['ta']+comp['delt'])*(1+z)

		# If num_bins is None, calculate the number of energy bins to ensure approximately 10 bins per order of magnitude 
		if num_bins is None:
			num_bins = int( 100*np.log10(Emax/Emin) )

		# Make array to store spectrum and initialize the summation of the spectrum as zero 
		spectrum_dE = np.zeros(shape=num_bins,dtype=[("ENERGY",float),("RATE",float)])
		spectrum_sum = 0.

		# Create energy axis for the spectrum
		energy_axis = np.logspace(np.log10(Emin),np.log10(Emax),num_bins)
		spectrum_dE['ENERGY'] = energy_axis

		# If comp_num = None, then include contributions from all stored components
		if comp_num == None:
			# For each component stored in our Emission class, sum the spectral contribution
			for comp in self.components:
				# Grab only the emission which occurred between Tmin and Tmax (i.e., emission with has (ta + delt_a)*(1+z) overlapping with (Tmin,Tmax) )
				partial_comp = comp[(comp['ta']*(1+z)<=Tmax) & ( (comp['ta']+comp['delt'])*(1+z)>=Tmin)] 
				
				add_spec_dE, add_spec_sum = self.__calc_spectrum_contr(partial_comp,energy_axis,z,sum_only=sum_only,nuFnu=nuFnu)
				spectrum_dE['RATE'] += add_spec_dE
				spectrum_sum += add_spec_sum

			return spectrum_dE, spectrum_sum

		# Else, record the spectrum for the particular component 
		elif comp_num is not None:
			# For each component stored in our Emission class, sum the spectral contribution
			comp = self.components[comp_num]
			# Grab only the emission which occurred between Tmin and Tmax (i.e., emission with has (ta + delt_a)*(1+z) overlapping with (Tmin,Tmax) )
			partial_comp = comp[(comp['ta']*(1+z)<=Tmax) & ( (comp['ta']+comp['delt'])*(1+z)>=Tmin)] 

			add_spec_dE, add_spec_sum = self.__calc_spectrum_contr(partial_comp,energy_axis,z,sum_only=sum_only,nuFnu=nuFnu)
			spectrum_dE['RATE'] = add_spec_dE
			spectrum_sum = add_spec_sum

			return spectrum_dE, spectrum_sum

	def __calc_spectrum_contr(self,partial_comp,energy_axis,z,sum_only,nuFnu):
		"""
		A private method used by the make_spec() method to calculate the contribution of a particular component to the spectrum

		Attributes:
		partial_comp = The subset of the component that occurs within a desired time interval (this is found in make_spec())
		energy_axis = The energy axis over which to calculate the spectrum
		z = Redshift of the source, used to calculate the spectrum in the observer frame
		sum_only=False, Set this parameter to True if the user only wants the summation of the spectrum across the given time interval and energy range 
			(For instance, this is used to reduce the number of calculation when creating a light curve)
		nuFnu=False, reports the nuFnu spectrum instead of the photon spectrum

		"""

		# Make array to store spectrum and initialize the summation of the spectrum as zero 
		spectrum_dE = np.zeros(shape=len(energy_axis))
		spectrum_sum = 0.

		# Check which type of component this is. 
		# We do this by checking if the component has specific fields unique to component types

		# If the component contains a temperature field "T", then it is a thermal component 
		if 'T' in partial_comp.dtype.names:
			# For each elementary emission which occurred in this partial component,
			# find the contribution to the spectrum
			for i in range(len(partial_comp)):
				add_spec_dE, add_spec_sum = calc_therm_contr(partial_comp[i],energy_axis,z,sum_only=sum_only,nuFnu=nuFnu)
				spectrum_dE += add_spec_dE
				spectrum_sum += add_spec_sum

		# If the component contains an alpha_synchrotron field "asyn", then it is a synchrotron component
		if 'asyn' in partial_comp.dtype.names:
			partial_comp = partial_comp[(partial_comp['tau']<1) & (partial_comp['relvel']>0.1)]
			# For each elementary emission which occurred in this partial component,
			# find the contribution to the spectrum
			for i in range(len(partial_comp)):
				add_spec_dE, add_spec_sum = calc_synch_contr(partial_comp[i],energy_axis,z,sum_only=sum_only,nuFnu=nuFnu)
				spectrum_dE += add_spec_dE
				spectrum_sum += add_spec_sum

		return spectrum_dE, spectrum_sum

	def make_light_curve(self,comp_num=None,Tmin=None,Tmax=None,dt=0.05,Emin=0.1,Emax=5e5,z=0):
		"""
		Method to make a light curve from the emission of the specified component(s). The light curve is calculated in the observer frame.

		Attributes:
		comp_num=None, User given integer to specify the component to use, if None is given all components will be summed
		
		Tmin=None, Beginning of the time interval over which the spectra will be computed
		Tmax=None, End of the time interval over which the spectra will be computed
		dt=0.05, Time resolution of the light curve
		Emin=0.1, Beginning of the energy range over which the spectra will be computed
		Emax=5e5, End of the energy range over which the spectra will be computed
		num_bins=None, Number of bins to use in the summation. By default the number of energy bins will be close to 10 bins per order of magnitude.
		
		z=0, Redshift of the observed source
		"""
		
		# Find start and stop time if a time interval is not specified
		# Both start and stop times should be corrected for time dilation
		# The start time will be the first emission from any of the components		
		if Tmin is None:
			Tmin=1e99
			for comp in self.components:
				if Tmin > comp['ta'][0]*(1+z):
					Tmin = comp['ta'][0]*(1+z)
		# The stop time is the final (emission start + emission duration)
		if Tmax is None:
			Tmax=-1e99
			for comp in self.components:
				if Tmax < np.max(comp['ta']+comp['delt'])*(1+z):
					Tmax = np.max(comp['ta']+comp['delt'])*(1+z)
	
		# Make time array is in the observer frame.
		time_arr = np.arange(start=Tmin,stop=Tmax+5,step=dt) # Add 5 seconds after the light curve, for better visualization.
		# Make light curve array, both a time axis and count rate axis
		light_curve = np.zeros(shape=len(time_arr),dtype=[('TIME',float),('RATE',float)])
		light_curve['TIME'] = time_arr
		
		# Move through each time bin
		for idt in range(len(light_curve)-1):
			# Sum the spectral contribution within this time bin
			light_curve['RATE'][idt] = self.make_spec(comp_num=comp_num, Tmin=light_curve['TIME'][idt], Tmax=light_curve['TIME'][idt+1], Emin=Emin, Emax=Emax,z=z,sum_only=True)[1]

		return light_curve

	def plot_spectrum(self, ax, comp_num=None, plot_comps=False, nuFnu=True, Tmin=None,Tmax=None,Emin=0.1,Emax=5e5,num_bins=None,z=0,fontsize=14,fontweight='bold'):
		"""
		Method to plot the spectrum from the specified component(s). The spectrum is calculated in the observer frame.

		Attributes:
		comp_num=None, User given integer to specify the component to use, if None is given all components will be summed
		plot_comps=False, Whether or not to plot each spectral component in addition to the total spectrum 

		Tmin=None, Beginning of the time interval over which the spectra will be computed
		Tmax=None, End of the time interval over which the spectra will be computed
		Emin=0.1, Beginning of the energy range over which the spectra will be computed
		Emax=5e5, End of the energy range over which the spectra will be computed
		num_bins=None, Number of bins to use in the summation. By default the number of energy bins will be close to 10 bins per order of magnitude.
		
		z=0, Redshift of the observed source
		
		fontsize=14, Font size of the labels and ticks of the plot 
		fontweight='bold', Font weight of the labels and ticks of the plot 
		"""
		
		if comp_num is not None:
			spectrum = self.make_spec(comp_num=comp_num,Tmin=Tmin,Tmax=Tmax,num_bins=num_bins,Emin=Emin,Emax=Emax,z=z,nuFnu=nuFnu)[0]			
			ax.plot(spectrum['ENERGY'],spectrum['RATE'],label='Comp. {}'.format(comp_num))
		
		else:
			spectrum = self.make_spec(Tmin=Tmin,Tmax=Tmax,num_bins=num_bins,Emin=Emin,Emax=Emax,z=z,nuFnu=nuFnu)[0]
			ax.plot(spectrum['ENERGY'],spectrum['RATE'],label='Total')

			if plot_comps is True:
				for i in range(len(self.components)):
					comp_spec = self.make_spec(comp_num=i,Tmin=Tmin,Tmax=Tmax,num_bins=num_bins,Emin=Emin,Emax=Emax,z=z,nuFnu=nuFnu)[0]
					ax.plot(comp_spec['ENERGY'],comp_spec['RATE'],label='Comp. {}'.format(i))

			if nuFnu is True:
				ax.set_ylabel(r'E$^2$N(E)',fontsize=fontsize,fontweight=fontweight)
			else:
				ax.set_ylabel(r'N(E)',fontsize=fontsize,fontweight=fontweight)

	def plot_light_curve(self, ax,comp_num=None, plot_comps=True, Tmin=None,Tmax=None,dt=0.05, Emin=0.1,Emax=5e5,z=0):
		"""
		Method to plot the light curve from the specified component(s). The light curve is calculated in the observer frame.

		Attributes:
		comp_num=None, User given integer to specify the component to use, if None is given all components will be summed
		plot_comps=False, Whether or not to plot each spectral component in addition to the total spectrum 

		Tmin=None, Beginning of the time interval over which the spectra will be computed
		Tmax=None, End of the time interval over which the spectra will be computed
		dt=0.05, Time resolution of the light curve
		Emin=0.1, Beginning of the energy range over which the spectra will be computed
		Emax=5e5, End of the energy range over which the spectra will be computed
		
		z=0, Redshift of the observed source
		"""

		# Plot complete light curvce
		if comp_num is not None:
			light_curve = self.make_light_curve(comp_num=comp_num,Tmin=Tmin, Tmax=Tmax, Emin=Emin, Emax=Emax, dt=dt,z=z)
			ax.scatter(light_curve['TIME'],light_curve['RATE'],label='Total')

		else:
			light_curve = self.make_light_curve(Tmin=Tmin, Tmax=Tmax, Emin=Emin, Emax=Emax, dt=dt,z=z)
			ax.scatter(light_curve['TIME'],light_curve['RATE'],label='Total')

			# Plot component light curves
			if plot_comps is True:
				for comp_num in range(len(self.components)):
					light_curve = self.make_light_curve(comp_num=comp_num, Tmin=Tmin, Tmax=Tmax, Emin=Emin, Emax=Emax, dt=dt,z=z)
					ax.scatter(light_curve['TIME'],light_curve['RATE'],label='Comp. {}'.format(comp_num))

def calc_therm_contr(partial_comp,energy_axis,z,sum_only=True,nuFnu=False):
	"""
	Method to calculate a thermal spectrum for a given flux and temperature.
	Currently, the thermal spectrum is assumed to be a cut-off power law with low energy index alpha=-0.4.

	Attributes:
	partial_comp,
	energy_axis,
	z = Redshift of the source, used to calculate the spectrum in the observer frame
	sum_only=False, Set this parameter to True if the user only wants the summation of the spectrum across the given time interval and energy range 
		(For instance, this is used to reduce the number of calculation when creating a light curve)
	nuFnu=False, reports the nuFnu spectrum instead of the photon spectrum
	"""
	
	# Put the temperature in the observer frame
	T_obs = partial_comp['T']/(1+z)

	# Calculate the normalization of the spectrum 
	norm = 1 / quad(lambda x, temp: x*thermal(x,temp), a=0, b=np.Infinity, args=(T_obs) )[0]

	# Make array to store spectrum and initialize the summation of the spectrum as zero 
	spectrum_dE = np.zeros(shape=len(energy_axis))
	spectrum_sum = 0.

	# If the nuFnu spectrum is desired:
	if nuFnu is True:
		# Total spectrum across the entire energy range:
		spectrum_sum = np.sum((energy_axis[1:]-energy_axis[0:-1]) * energy_axis[0:-1]**2 * thermal(energy_axis[0:-1],T_obs))
		# Apply the necessary normalizations and conversions
		spectrum_sum *=  norm * partial_comp['Flux']

		# To save time we can omit calculating the spectrum at each energy bin
		# else we enter the if statement below
		if sum_only is False:
			# For each energy bin along the energy axis, calculate the contribution to the spectrum
			spectrum_dE = energy_axis**2 * thermal(energy_axis,T_obs)			
			# Apply the necessary normalizations and conversions
			spectrum_dE *= norm * partial_comp['Flux']
		
		return spectrum_dE, spectrum_sum

	else:
		# Total spectrum across the entire energy range:
		spectrum_sum = np.sum((energy_axis[1:]-energy_axis[0:-1])*thermal(energy_axis[0:-1],T_obs))
		# Apply the necessary normalizations and conversions
		spectrum_sum *=  norm * partial_comp['Flux']

		# To save time we can omit calculating the spectrum at each energy bin
		# else we enter the if statement below
		if sum_only is False:
			# For each energy bin along the energy axis, calculate the contribution to the spectrum
			spectrum_dE = thermal(energy_axis,T_obs)
			# Apply the necessary normalizations and conversions
			spectrum_dE *= norm * partial_comp['Flux']
		
		return spectrum_dE, spectrum_sum

def calc_synch_contr(partial_comp,energy_axis,z,sum_only=True,nuFnu=False,alpha=-1,beta=-2.5):
	"""
	Method to calculate a synchrotron spectrum for a given dissipated energy and Synchrotron Energy.
	Currently, the synchrotron spectrum is assumed to be a smoothly broken power law with low energy index alpha and high energy index beta.

	Attributes:
	partial_comp,
	energy_axis,
	z = Redshift of the source, used to calculate the spectrum in the observer frame
	sum_only=False, Set this parameter to True if the user only wants the summation of the spectrum across the given time interval and energy range 
		(For instance, this is used to reduce the number of calculation when creating a light curve)
	nuFnu=False, reports the nuFnu spectrum instead of the photon spectrum
	alpha=-2/3, The low energy power law index
	beta=-2.5, The high energy power law index
	"""
	
	# Put the Synchrotron energy into the observer frame 
	E_sync_obs = partial_comp['Esyn']/(1+z)

	# Calculate the normalization of the spectrum 
	norm = 1 / quad(lambda x, esyn,gamma_1,gamma_2: x*synchrotron(x, esyn,gamma_1,gamma_2),a=0,b=np.Infinity,args=(E_sync_obs,alpha,beta) )[0]

	# Make array to store spectrum and initialize the summation of the spectrum as zero 
	spectrum_dE = np.zeros(shape=len(energy_axis))
	spectrum_sum = 0.
	
	# If the nuFnu spectrum is desired:
	if nuFnu is True:
		# Total spectrum across the entire energy range:
		spectrum_sum = np.sum((energy_axis[1:]-energy_axis[0:-1]) * energy_axis[0:-1]**2 * synchrotron(energy_axis[0:-1], E_sync_obs, alpha, beta))
		# Apply the necessary normalizations and conversions
		spectrum_sum *= cc.kev_to_erg * norm * partial_comp['e'] / (partial_comp['delt']*(1+z))

		# To save time we can omit calculating the spectrum at each energy bin
		# else we enter the if statement below
		if sum_only is False:
			# For each energy bin along the energy axis, calculate the contribution to the spectrum
			spectrum_dE = energy_axis**2 * synchrotron(energy_axis, E_sync_obs, alpha=alpha, beta=beta)
			# Apply the necessary normalizations and conversions
			spectrum_dE *= cc.kev_to_erg * norm * partial_comp['e'] / (partial_comp['delt']*(1+z))

		return spectrum_dE, spectrum_sum
	
	else:
		# Total spectrum across the entire energy range:
		spectrum_sum = np.sum( (energy_axis[1:]-energy_axis[0:-1]) * synchrotron(energy_axis[0:-1], E_sync_obs, alpha=alpha, beta=beta))
		# Apply the necessary normalizations and conversions
		spectrum_sum *= norm * partial_comp['e'] / (partial_comp['delt']*(1+z))

		# To save time we can omit calculating the spectrum at each energy bin
		# else we enter the if statement below
		if sum_only is False:
			# For each energy bin along the energy axis, calculate the contribution to the spectrum
			spectrum_dE = synchrotron(energy_axis, E_sync_obs, alpha=alpha, beta=beta) 
			# Apply the necessary normalizations and conversions
			spectrum_dE *= norm * partial_comp['e'] / (partial_comp['delt']*(1+z))

		return spectrum_dE, spectrum_sum

def thermal(x,temp,alpha=0.4):
	"""
	Method to produce a Planck spectrum over a given energy range given and a specified temperature
	
	Attributes:
	x: The energy value(s) to compute the spectrum for
	temp: the temperature of the thermal source
	alpha: photon spectral index, default = 0.4 (a value which differs from a typical Planck spectrum to account for relativistic corrections)
	"""

	kb = cc.kb_kev # Boltzmann constant, in units keV/K	
	
	return np.power(x/(kb*temp),1+alpha)/(np.exp(x/(kb*temp))-1)

	# E_peak = 3.9*kb*temp # peak energy, in units of keV 
	# return (1/E_peak) * cutpl(x=x,ep=E_peak,alpha=0.4)

def synchrotron(x,Esyn,alpha=-2/3,beta=-2.5):
	"""
	Method to produce a synchrotron spectrum over a given energy range and specified synchrotron energy and energy dissipated
	Currently, the synchrotron spectrum is assumed to be a broken power law with low energy index alpha and high energy index beta.
	
	Attributes:
	x: The energy value(s) to compute the spectrum for
	Esyn: The typical synchrotron energy of the emission
	alpha = -2/3, The low energy power law index
	beta = -2.5, The high energy power law index
	"""

	# return (1/Esyn) * bpl(x,Esyn,alpha=alpha,beta=beta)
	return (1/Esyn) * band(x,Esyn,alpha=alpha,beta=beta)

def bpl(x,ep,alpha=-2/3,beta=-2.5):
	"""
	Method to produce a sharply broken power law with a specified peak energy and over a given energy range (or particular energy value)
	
	Attributes:
	x: The energy value(s) to compute the spectrum for
	ep: The peak energy of the spectrum 
	alpha=-2/3, The low energy power law index
	beta=-2.5, The high energy power law index
	"""

	# If an array of energies was supplied for computation, it will have a length attribute.
	# Calculate the spectrum over the specified energy range
	if hasattr(x,'__len__'):
		
		# Create an array to store the spectrum
		f = np.zeros(shape=len(x))
	
		# If the peak energy is greater than the energy range maximum, calculate the entire spectrum in the low energy power law
		if ep>x[-1]:
			f = np.power(x/ep,alpha)
		else:
			# Calculate the index where the energy axis surpasses the peak energy
			ind_xp = np.argmax(x>ep)

			# Calculate the power law:
			# Below the peak energy use the low energy index
			f[:ind_xp] = np.power(x[:ind_xp]/ep,alpha)
			
			# Above the peak energy use the high energy index
			f[ind_xp:] = np.power(x[ind_xp:]/ep,beta)

		return f

	# Compute the power law spectrum at a particular energy
	else:
		# If the energy is below the peak energy
		if x < ep:
			return np.power(x/ep,alpha)
		# If the energy is above the peak energy
		elif x>= ep:
			return np.power(x/ep,beta)

def cutpl(x,ep,alpha=-2/3):
	"""
	Method to produce a cut off power law over a given energy range
	
	Attributes:
	x: The energy value(s) to compute the spectrum for
	ep: The peak energy of the spectrum 
	alpha: Low energy power law index
	"""
	return np.power(x/ep,alpha) * np.exp(-x/ep) 

def band(x,e0,alpha=-2/3,beta=-2.5):
	"""
	Method to produce a Band function with a specified peak energy and over a given energy range (or particular energy value)
	
	Attributes:
	x: The energy value(s) to compute the spectrum for
	ep: The peak energy of the spectrum 
	alpha=-2/3, The low energy power law index
	beta=-2.5, The high energy power law index
	"""

	if hasattr(x,'__len__'):
		
		# Create an array to store the spectrum
		f = np.zeros(shape=len(x))

		# If the peak energy is greater than the energy range maximum, calculate the entire spectrum in the low energy power law
		if e0>x[-1]:
			f = (x/100)**alpha * np.exp(-x/e0)
		else:
			# Calculate the index where the energy axis surpasses the peak energy
			ind_xp = np.argmax(x>e0*(alpha-beta))

			# Calculate the power law:
			# Below the peak energy use the low energy index
			f[:ind_xp] = (x[:ind_xp]/100)**alpha * np.exp(-x[:ind_xp]/e0)
			
			# Above the peak energy use the high energy index
			f[ind_xp:] = np.power((alpha-beta)*e0/100,alpha-beta)*np.exp(beta-alpha)*np.power(x[ind_xp:]/100,beta)

		return f
	
	# Compute the power law spectrum at a particular energy
	else:
		# If the energy is below the peak energy
		if x <= (alpha-beta)*e0:
			return (x/100)**alpha * np.exp(-x/e0)
		# If the energy is above the peak energy
		else:
			return np.power((alpha-beta)*e0/100,alpha-beta)*np.exp(beta-alpha)*np.power(x/100,beta)
