"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-03-16


This code defines an Emission class which contains all the emission produced during the prompt jet simulation

This script contains the various spectral functions used to produced the spectra for the burst.
Additionally, this script has methods for creating plots of the spectra and light curve (when supplied with a spectrum)


"""

import numpy as np
import matplotlib.pyplot as plt
import cosmologicalconstants as cc



class Emission(object):
	"""
	Spectrum class.
	"""

	def __init__(self):
		"""
		Defines the default parameters of a spectrum.
		"""
		self.spec_therm = None 
		self.spec_synch = None 

	def add_synch_contribution(self,te,ta,asyn,Beq,gammae,Esyn,gammar,e,delt):
		"""
		Add a contribution to the synchrotron spectrum
		"""
		# If this is the first contribution to the spectrum, the spectrum must be initialized. 
		if self.spec_synch is None:
			self.spec_synch = np.array((te,ta,asyn,Beq,gammae,Esyn,gammar,e,delt), dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float),('gammar',float),('e',float),('delt',float)])
		# Otherwise, just append to the already created spectrum
		else:
			self.spec_synch = np.append(self.spec_synch, np.array((te,ta,asyn,Beq,gammae,Esyn,gammar,e,delt),dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float),('gammar',float),('e',float),('delt',float)])) 

	def add_therm_contribution(self,te,ta,delt,T,L):
		"""
		Add a contribution to the thermal spectrum
		"""

		if self.spec_therm is None:
			self.spec_therm = np.array((te,ta,delt,T,L), dtype=[('te',float),('ta',float),('delt',float),('T',float),('L',float)])
		# Otherwise, just append to the already created spectrum
		else:
			self.spec_therm = np.append(self.spec_therm, np.array((te,ta,delt,T,L),dtype=[('te',float),('ta',float),('delt',float),('T',float),('L',float)] ) )


	def make_tot_spec(self,spec_therm=None,spec_synch=None,dt=0.01,Emin=4e4,Emax=5e6,dE=100):
		"""
		Method to create a spectrum from all active components as a function of time

		Attributes:
		dt = time resolution, seconds 
		Emin, Emax = specifies the min and max of the desired energy range, eV
		"""


		# Record which components are active
		if (spec_synch is None) & (self.spec_synch is not None):
			spec_synch = self.spec_synch
		elif (spec_therm is None) & (self.spec_therm is not None):
			spec_therm = self.spec_therm

		# When does the spectrum start and stop?
		tstart = 0
		tend = 0 # will be found from the active spectral components
		if spec_synch is not None:
			if np.max(spec_synch['ta']+spec_synch['delt']) > tend:
				tend  = np.max(spec_synch['ta']+spec_synch['delt'])
		elif spec_therm is not None:
			if np.max(spec_therm['ta']+spec_therm['delt']) > tend:
				tend  = np.max(spec_therm['ta']+spec_therm['delt'])
		else: 
			# all components are empty
			print("At least one spectral component must be provided.")
			return 0

		# Make time axis
		time_ax = np.arange(start=tstart,stop=tend+dt,step=dt)
		# Make spec axis
		num_bins = int((Emax-Emin)/dE) # Number of bins along the energy axis
		en_ax = np.logspace(np.log10(Emin),np.log10(Emax),num_bins) # Create energy axis of spectra
		spec_axis = np.zeros(shape=(len(time_ax),len(en_ax)) )


		if spec_synch is not None:
			for i in range(len(spec_synch)):
				print(i)
				time_ind_start = np.argmax(time_ax>spec_synch['ta'][i])-1
				time_ind_end = np.argmax(time_ax>spec_synch['ta'][i]+spec_synch['delt'][i])


				spec_axis[time_ind_start:time_ind_end,:]+=synchrotron(en_ax,spec_synch['Esyn'][i]*10,spec_synch['e'][i]/spec_synch['delt'][i]/dt)

		if spec_therm is not None:
			for i in range(len(spec_therm)):
				print(i)
				time_ind_start = np.argmax(time_ax>spec_therm['ta'][i])-1
				time_ind_end = np.argmax(time_ax>spec_therm['ta'][i]+spec_therm['delt'][i])

				spec_axis[time_ind_start:time_ind_end,:]+=thermal(en_ax,spec_therm['T'][i])


		return time_ax,en_ax, spec_axis


def thermal(energy_bins,temp):
	"""
	Method to produce a thermal spectrum over a given energy range given and a specified temperature and
	"""

	kb = cc.kb # Boltzmann constant

	# Initialize array for thermal spectrum 
	dNE_therm = np.zeros(shape=len(energy_bins))

	# Index at maximum energy for thermal, this is at 2 MeV (pair production)
	pp_ind = np.argmax(energy_bins>2*1e6)
	# Index at critical point between power law and exp cut off
	break_ind = np.argmax(energy_bins>4*kb*temp)

	dNE_therm[0:break_ind] +=2* energy_bins[0:break_ind]**1.4 / (cc.h**2 * cc.c**2) / (np.exp(energy_bins[0:break_ind]/kb/temp)-1)
	if pp_ind == 0:
		dNE_therm[break_ind:] += 2* energy_bins[break_ind:]**1.4 * np.exp(-energy_bins[break_ind:]/kb/temp) / (cc.h**2 * cc.c**2)
	else: 
		dNE_therm[break_ind:pp_ind] += 2* energy_bins[break_ind:pp_ind]**1.4 * np.exp(-energy_bins[break_ind:pp_ind]/kb/temp) / (cc.h**2 * cc.c**2)

	return dNE_therm

def synchrotron(energy_bins,Esyn,endiss):
	"""
	Method to produce a synchrotron spectrum over a given energy range and specified synchrotron energy and energy dissipated
	"""

	dNE_sync = np.zeros(shape=len(energy_bins))

	# Which index is the break at?
	break_ind = np.argmax(energy_bins>Esyn)

	# Before break
	x= -2/3
	dNE_sync[0:break_ind] += (endiss/Esyn)*np.power(energy_bins[0:break_ind]/Esyn, x)

	# After break
	x=-2.5
	dNE_sync[break_ind:] += (endiss/Esyn)*np.power(energy_bins[break_ind:]/Esyn, x)


	return dNE_sync


def plot_spectrum( ax, spec_therm=None, spec_synch=None,nuFnu=True, num_bins=1000,emin=100,emax=1e9):
	"""
	Method to plot the stored spectrum over a given energy range
	"""

	# Specify energy range
	emin = emin # eV 
	emax = emax # eV 
	# Make energy bins (x-axis)
	enlogbins = np.logspace(np.log10(emin),np.log10(emax),num_bins)

	# Initialize array for total spectrum 
	dNE = np.zeros(shape=len(enlogbins))

	# If all spectrum types are empty, then there is nothing to plot
	if spec_therm is None and spec_synch is None:
		print("You must supply a spectrum to be plotted.")
		return 0 

	# If a thermal spectrum has been supplied
	if spec_therm is not None:

		# Initialize array for thermal spectrum 
		dNE_therm = np.zeros(shape=len(enlogbins))

		for j in range(len(spec_therm)):
			therm_contr = thermal(enlogbins,spec_therm['T'][j])
			dNE_therm += therm_contr
			dNE += therm_contr

	# If a synchrotron spectrum has been supplied 
	if spec_synch is not None:

		# Initialize array for synchrotron spectrum 
		dNE_sync = np.zeros(shape=len(enlogbins))

		for j in range(len(spec_synch['Esyn'])):		
			synch_contr = synchrotron(enlogbins,spec_synch['Esyn'][j],spec_synch['e'][j])
			dNE_sync += synch_contr
			dNE += synch_contr


	# For axis labels
	fontsize=14
	fontweight='bold'

	if nuFnu is True:
		# Plot each spectral component and the total spectrum
		ax.plot(enlogbins,dNE*enlogbins**2,label='Total')
		if spec_therm is not None:
			ax.plot(enlogbins,dNE_therm*enlogbins**2,label='Therm')
		if spec_synch is not None:
			ax.plot(enlogbins,dNE_sync*enlogbins**2,label='Synch')
		ax.set_ylabel(r'E$^2$N(E)',fontsize=fontsize,fontweight=fontweight)
	if nuFnu is False:
		# Plot each spectral component and the total spectrum
		ax.plot(enlogbins,dNE,label='Total')
		if spec_therm is not None:
			ax.plot(enlogbins,dNE_therm,label='Therm')
		if spec_synch is not None:
			ax.plot(enlogbins,dNE_sync,label='Synch')
		ax.set_ylabel(r'N(E)',fontsize=fontsize,fontweight=fontweight)


	# Plot aesthetics
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('E (eV)',fontsize=fontsize,fontweight=fontweight)

	for tick in ax.xaxis.get_major_ticks():
	    tick.label1.set_fontsize(fontsize=fontsize)
	    tick.label1.set_fontweight(fontweight)
	for tick in ax.yaxis.get_major_ticks():
	    tick.label1.set_fontsize(fontsize=fontsize)
	    tick.label1.set_fontweight(fontweight)


def plot_light_curve( ax, em_therm=None, em_synch=None, dt=0.03, emin=8000,emax=4e6):
	"""
	Method to plot the simulated light curve
	"""

	# If all spectrum types are empty, then there is nothing to plot
	if em_therm is None and em_synch is None:
		print("You must supply an emission file to be plotted.")
		return 0 

	# Specify energy range
	emin = emin # eV 
	emax = emax # eV

	num_bins=1000
	enlogbins = np.logspace(np.log10(emin),np.log10(emax),num_bins)

	# Time range in observer frame (set by the first and last emission arrival time)	
	T_min = 0
	
	T_max_therm = 0
	T_max_synch = 0
	if em_therm is not None:
		T_max_therm = np.max(em_therm['ta'])
	if em_synch is not None:
		T_max_synch = np.max(em_synch['ta'])

	T_max = np.max([T_max_therm, T_max_synch])

	# Observer time axis
	time_arr = np.arange(start=T_min,stop=T_max,step=dt)

	# Initialize count array
	count_arr = np.zeros(shape=len(time_arr))
	if em_therm is not None:
		count_arr_therm = np.zeros(shape=len(time_arr))
	if em_synch is not None:
		count_arr_synch = np.zeros(shape=len(time_arr))

	
	# For each observer time bin, find the emission contribution
	for i in range(len(time_arr)-1):

		# Initialize array for total spectrum 
		dNE = np.zeros(shape=len(enlogbins))
	
		# If a thermal spectrum has been supplied
		if em_therm is not None:
			dt_em_therm = em_therm[(em_therm['ta']>time_arr[i]) & (em_therm['ta']<time_arr[i+1]) ]

			# Initialize array for thermal spectrum 
			dNE_therm = np.zeros(shape=len(enlogbins))

			for j in range(len(dt_em_therm)):
				therm_contr = thermal(enlogbins,dt_em_therm['T'][j])
				dNE_therm += therm_contr
				dNE += therm_contr

			# Take the sum of the spectrum and add it to the count array in the correct time bin
			count_arr_therm[i] += np.sum(dNE_therm)
			count_arr[i] += np.sum(dNE)

		# If a synchrotron spectrum has been supplied 
		if em_synch is not None:

			dt_em_synch = em_synch[(em_synch['ta']>time_arr[i]) & (em_synch['ta']<time_arr[i+1]) ]

			# Initialize array for synchrotron spectrum 
			dNE_synch = np.zeros(shape=len(enlogbins))

			for j in range(len(dt_em_synch['Esyn'])):
				synch_contr = synchrotron(enlogbins,dt_em_synch['Esyn'][j],dt_em_synch['e'][j])
				dNE_synch += synch_contr
				dNE += synch_contr

			# Take the sum of the spectrum and add it to the count array in the correct time bin
			count_arr_synch[i] += np.sum(dNE_synch)
			count_arr[i] += np.sum(dNE)


	# For axis labels
	fontsize=14
	fontweight='bold'

	if em_therm is not None:
		ax.step(time_arr,count_arr_therm,label='Thermal',where='mid')
	if em_synch is not None:
		ax.step(time_arr,count_arr_synch,label='Synchrotron',where='mid')
	
	ax.step(time_arr,count_arr,label='Total',where='mid')

	# Plot aesthetics
	ax.set_xlabel('Time (sec)',fontsize=fontsize,fontweight=fontweight)
	ax.set_ylabel('Counts/Sec',fontsize=fontsize,fontweight=fontweight)

	for tick in ax.xaxis.get_major_ticks():
	    tick.label1.set_fontsize(fontsize=fontsize)
	    tick.label1.set_fontweight(fontweight)
	for tick in ax.yaxis.get_major_ticks():
	    tick.label1.set_fontsize(fontsize=fontsize)
	    tick.label1.set_fontweight(fontweight)