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

# from scipy.integrate import quad


class Emission(object):
	"""
	Spectrum class.
	"""

	def __init__(self,components=None,Tmin=None,Tmax=None,Emin=0.1,Emax=5e5,dE=1000):
		"""
		Defines the default parameters of a spectrum.
		"""
		self.components = []
		if len(self.components) > 0:
			self.make_spec()

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
			spec_synch = np.genfromtxt('sim_results/ordlor_spectrum_synch.txt',dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float),('gammar',float),('e',float),('delt',float),('tau',float),('relvel',float)])
			self.components.append(spec_synch)

			# Return the spectral component 
			return spec_synch

		elif first_line_len == 5:
			# Load spectrum component and add it into the components list
			spec_therm = np.genfromtxt('sim_results/ordlor_spectrum_therm.txt',dtype=[('te',float),('ta',float),('delt',float),('T',float),('Flux',float)])
			self.components.append(spec_therm)

			# Return the spectral component 
			return spec_therm


	def make_spec(self,comp_num=None,Tmin=None,Tmax=None,Emin=0.1,Emax=5e5,num_bins=1000,z=0):
		"""
		Method to create a spectrum from all active components as a function of time

		Attributes:
		dt = time resolution, seconds 
		Emin, Emax = specifies the min and max of the desired energy range, keV
		"""

		# Find start and stop time if a time interval is not specified
		if Tmin is None:
			Tmin=1e99
			for comp in self.components:
				if Tmin > comp['ta'][0]*(1+z):
					Tmin = comp['ta'][0]*(1+z)

		if Tmax is None:
			Tmax=-1e99
			for comp in self.components:
				if Tmax < np.max(comp['ta']+comp['delt'])*(1+z):
					Tmax = np.max(comp['ta']+comp['delt'])*(1+z)

		# Make array to store spectrum
		spectrum = np.zeros(shape=num_bins,dtype=[("ENERGY",float),("RATE",float)])

		# Create energy axis of spectrum
		spectrum['ENERGY'] = np.logspace(np.log10(Emin),np.log10(Emax),num_bins) 


		# If comp_num = None, then include contributions from all stored components
		if comp_num == None:
			# For each component stored in our Emission class, sum the spectral contribution
			for comp in self.components:
			# Grab only the emission which occurred between Tmin and Tmax (i.e., emission with has (ta + delt_a)*(1+z) overlapping with (Tmin,Tmax) )
				partial_comp = comp[(comp['ta']*(1+z)<=Tmax) & ( (comp['ta']+comp['delt'])*(1+z)>=Tmin)] 


				# Check which type of component this is. We do this by checking if the component has specific fields unique to component types
				if 'T' in partial_comp.dtype.names:
					# For each elementary in this partial component, find the contribution to the spectrum (and normalize by the duration of the emission)
					for i in range(len(partial_comp)):
						spectrum['RATE']+= 1e20*( partial_comp['Flux'][i]/4/np.pi/ut.lum_dis(z)**2 )*thermal(spectrum['ENERGY'],partial_comp['T'][i]/(1+z))/(partial_comp['delt'][i]*(1+z))

				if 'asyn' in partial_comp.dtype.names:
					# For each elementary in this partial component, find the contribution to the spectrum (and normalize by the duration of the emission)
					for i in range(len(partial_comp)):
						spectrum['RATE']+= partial_comp['e'][i]*synchrotron(spectrum['ENERGY'],partial_comp['Esyn'][i]/(1+z))/(partial_comp['delt'][i]*(1+z))

			return spectrum
		# Else, record the spectrum for the particular component 
		elif comp_num is not None:
			# For each component stored in our Emission class, sum the spectral contribution
			comp = self.components[comp_num]
			# Grab only the emission which occurred between Tmin and Tmax (i.e., emission with has (ta + delt_a)*(1+z) overlapping with (Tmin,Tmax) )
			partial_comp = comp[(comp['ta']*(1+z)<Tmax) & ((comp['ta']+comp['delt'])*(1+z)>Tmin)]


			# Check which type of component this is. We do this by checking if the component has specific fields unique to component types
			if 'T' in partial_comp.dtype.names:
				# For each elementary in this partial component, find the contribution to the spectrum (and normalize by the duration of the emission)
				for i in range(len(partial_comp)):
					spectrum['RATE']+= 1e20*(partial_comp['Flux'][i]/4/np.pi/ut.lum_dis(z)**2 )*thermal(spectrum['ENERGY'],partial_comp['T'][i]/(1+z))/(partial_comp['delt'][i]*(1+z))

			if 'asyn' in partial_comp.dtype.names:
				# For each elementary in this partial component, find the contribution to the spectrum (and normalize by the duration of the emission)
				for i in range(len(partial_comp)):
					spectrum['RATE']+= partial_comp['e'][i]*synchrotron(spectrum['ENERGY'],partial_comp['Esyn'][i]/(1+z))/(partial_comp['delt'][i]*(1+z))

			return spectrum

	def make_light_curve(self,comp_num=None,Tmin=None,Tmax=None,dt=0.05,Emin=0.1,Emax=5e5,z=0):
		# Find start and stop time if a time interval is not specified
		
		if Tmin is None:
			Tmin=1e99
			for comp in self.components:
				if Tmin > comp['ta'][0]*(1+z):
					Tmin = comp['ta'][0]*(1+z)

		if Tmax is None:
			Tmax=-1e99
			for comp in self.components:
				if Tmax < np.max(comp['ta']+comp['delt'])*(1+z):
					Tmax = np.max(comp['ta']+comp['delt'])*(1+z)
	
		# Make light curve array
		time_arr = np.arange(start=Tmin,stop=Tmax+5,step=dt) # Add 5 seconds after the light curve, This time array is in the observer frame.
		light_curve = np.zeros(shape=len(time_arr),dtype=[('TIME',float),('RATE',float)])
		light_curve['TIME'] = time_arr
		
		# If comp_num = None, then include contributions from all stored components
		if comp_num == None:
			# Move through each time bin
			for idt in range(len(light_curve)-1):
				# Sum the spectral contribution within this time bin
				light_curve['RATE'][idt] = np.sum(self.make_spec(comp_num=None, Tmin=light_curve['TIME'][idt], Tmax=light_curve['TIME'][idt+1], Emin=Emin, Emax=Emax,z=z)['RATE'])

		# Else, record the spectrum for the particular component 
		elif comp_num is not None:
			# Move through each time bin
			for idt in range(len(light_curve)-1):
				light_curve['RATE'][idt] = np.sum(self.make_spec(comp_num=comp_num, Tmin=light_curve['TIME'][idt], Tmax=light_curve['TIME'][idt+1], Emin=Emin, Emax=Emax,z=z)['RATE'])

		return light_curve

	
	def plot_spectrum(self, ax, comp_num=None, plot_comps=False, nuFnu=True, Tmin=None,Tmax=None,Emin=0.1,Emax=5e5,num_bins=1000,z=0):
		"""
		Method to plot the stored spectrum over a given energy range
		"""
	
		# For axis labels
		fontsize=14
		fontweight='bold'
	
		if comp_num is not None:
			spectrum = self.make_spec(comp_num=comp_num,Tmin=Tmin,Tmax=Tmax,num_bins=num_bins,Emin=Emin,Emax=Emax,z=z)
			
			if nuFnu is True:
				ax.plot(spectrum['ENERGY'],spectrum['RATE']*spectrum['ENERGY']**2,label='Comp. {}'.format(comp_num))
			if nuFnu is False:
				ax.plot(spectrum['ENERGY'],spectrum['RATE'],label='Comp. {}'.format(comp_num))
		
		else:
			spectrum = self.make_spec(Tmin=Tmin,Tmax=Tmax,num_bins=num_bins,Emin=Emin,Emax=Emax,z=z)

			if nuFnu is True:
				ax.plot(spectrum['ENERGY'],spectrum['RATE']*spectrum['ENERGY']**2,label='Total')
				if plot_comps is True:
					for i in range(len(self.components)):
						comp_spec = self.make_spec(comp_num=i,Tmin=Tmin,Tmax=Tmax,num_bins=num_bins,Emin=Emin,Emax=Emax,z=z)
						ax.plot(comp_spec['ENERGY'],comp_spec['RATE']*comp_spec['ENERGY']**2,label='Comp. {}'.format(i))
				ax.set_ylabel(r'E$^2$N(E)',fontsize=fontsize,fontweight=fontweight)
			if nuFnu is False:
				ax.plot(spectrum['ENERGY'],spectrum['RATE'],label='Total')
				if plot_comps is True:
					for i in range(len(self.components)):
						comp_spec = self.make_spec(comp_num=i,Tmin=Tmin,Tmax=Tmax,num_bins=num_bins,Emin=Emin,Emax=Emax,z=z)
						ax.plot(comp_spec['ENERGY'],comp_spec['RATE'],label='Comp. {}'.format(i))
				ax.set_ylabel(r'N(E)',fontsize=fontsize,fontweight=fontweight)
	
	
	def plot_light_curve(self, ax,comp_num=None, plot_comps=True, Tmin=None,Tmax=None,dt=0.05, Emin=0.1,Emax=5e5,z=0):
		"""
		Method to plot the simulated light curve
		"""

		# Plot complete light curvce
		if comp_num is not None:
			light_curve = self.make_light_curve(comp_num=comp_num,Tmin=Tmin, Tmax=Tmax, Emin=Emax, Emax=Emax, dt=dt,z=z)
			ax.step(light_curve['TIME'],light_curve['RATE'],label='Total')

		else:
			light_curve = self.make_light_curve(Tmin=Tmin, Tmax=Tmax, Emin=Emax, Emax=Emax, dt=dt,z=z)
			ax.step(light_curve['TIME'],light_curve['RATE'],label='Total')

			# Plot component light curves
			if plot_comps is True:
				for comp_num in range(len(self.components)):
					light_curve = self.make_light_curve(comp_num=comp_num, Tmin=Tmin, Tmax=Tmax, Emin=Emax, Emax=Emax, dt=dt,z=z)
					ax.step(light_curve['TIME'],light_curve['RATE'],label='Comp. {}'.format(comp_num))
	

def thermal(x,temp):
	"""
	Method to produce a thermal spectrum over a given energy range given and a specified temperature and
	"""

	kb = cc.kb/1e3 # Boltzmann constant, in units keV/K

	if hasattr(x,'__len__'):
		f = np.zeros(shape=len(x))

		# If below the break energy:
		# low_ind_i=0
		low_ind_f=np.argmax(x>4*kb*temp)
		f[:low_ind_f] = 2* x[:low_ind_f]**1.4 / (cc.h**2 * cc.c**2) / (np.exp(x[:low_ind_f]/kb/temp)-1)
		# If below the Pair-Production energy but above the break energy, around 2.5 MeV
		high_ind_f=np.argmax(x>2.5*1e3)
		f[low_ind_f:high_ind_f] = 2* x[low_ind_f:high_ind_f]**1.4 * np.exp(-x[low_ind_f:high_ind_f]/kb/temp) / (cc.h**2 * cc.c**2)
		# If above the PP energy, then the emission is just zero, which is the default value of the array

		return f

	else:

		# If below the break energy:
		if x<3.9*kb*temp:
			return 2* x**1.4 / (cc.h**2 * cc.c**2) / (np.exp(x/kb/temp)-1)
		# If below the PP energy but above the break energy, around 2.5 MeV
		elif (x>=4*kb*temp) and (x<2.5*1e3):
			return 2* x**1.4 * np.exp(-x/kb/temp) / (cc.h**2 * cc.c**2)
		else:
			return 0

def synchrotron(x,Esyn):
	"""
	Method to produce a synchrotron spectrum over a given energy range and specified synchrotron energy and energy dissipated
	"""
	dNE_sync = (1/Esyn) * bpl(x,Esyn)

	return dNE_sync

def bpl(x,ep):
	if hasattr(x,'__len__'):
		f = np.zeros(shape=len(x))
		
		alpha = -2/3
		f[:np.argmax(x>ep)-1] = np.power(x[:np.argmax(x>ep)-1]/ep,alpha)
		beta= -2.5
		f[np.argmax(x>ep)-1:] = np.power(x[np.argmax(x>ep)-1:]/ep,beta)

		return f

	else:
		if x < ep:
			alpha = -2/3
			return np.power(x/ep,alpha)
		elif x>= ep:
			beta= -2.5
			return np.power(x/ep,beta)