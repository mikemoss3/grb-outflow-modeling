"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-03-02


This code defines a LightCurve object class. This class contains all relevant information and definitions of a light curve created by
GRB prompt emission.

"""

import numpy as np
import matplotlib.pyplot as plt
import cosmologicalconstants as cc

import spectrum as sp




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
				therm_contr = sp.thermal(enlogbins,dt_em_therm['T'][j])
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
				synch_contr = sp.synchrotron(enlogbins,dt_em_synch['Esyn'][j],dt_em_synch['e'][j])
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