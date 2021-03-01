"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2020-09-01


This code defines a Spectrum object class. This class contains all relevant information and definitions of a spectrum created by
GRB prompt emission.

"""

import numpy as np
import matplotlib.pyplot as plt
import cosmologicalconstants as cc


class Spectrum(object):
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

	def add_therm_contribution(self,te,ta,T,L):
		"""
		Add a contribution to the thermal spectrum
		"""

		if self.spec_therm is None:
			self.spec_therm = np.array((te,ta,T,L), dtype=[('te',float),('ta',float),('T',float),('L',float)])
		# Otherwise, just append to the already created spectrum
		else:
			self.spec_therm = np.append(self.spec_therm, np.array((te,ta,T,L),dtype=[('te',float),('ta',float),('T',float),('L',float)] ) ) 


def plot_spectrum( ax, spec_therm=None, spec_synch=None,nuFnu=True, num_bins=1000,emin=100,emax=1e9):
	"""
	Method to plot the stored spectrum
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
		
		kb = cc.kb # Boltzmann constant

		for j in range(len(spec_therm)):
			# a thermal spectrum should cutoff before 2 MeV because of pair opacity.
			for i in range(len( enlogbins[0:np.argmax(enlogbins>2*1e6) ] )):
			# for i in range(len(enlogbins)):
				if enlogbins[i] < 4*kb*spec_therm['T'][j]:
					val = 2* enlogbins[i]**1.4 / (cc.h**2 * cc.c**2) / (np.exp(enlogbins[i]/kb/spec_therm['T'][j])-1)
				else: 
					val = 2* enlogbins[i]**1.4 * np.exp(-enlogbins[i]/kb/spec_therm['T'][j]) / (cc.h**2 * cc.c**2)
				dNE_therm[i] += val
				dNE[i] += val

	# If a synchrotron spectrum has been supplied 
	if spec_synch is not None:

		# Initialize array for synchrotron spectrum 
		dNE_sync = np.zeros(shape=len(enlogbins))

		for j in range(len(spec_synch['Esyn'])):
			for i in range(len(enlogbins)):
				# Low energy power law
				if enlogbins[i] < spec_synch['Esyn'][j]:
					x = -2/3
				# High energy power law
				elif enlogbins[i] > spec_synch['Esyn'][j]:
					x = -2.5
				# val = (0.01*0.33*3e-3)*(spec_synch['e'][j]/spec_synch['Esyn'][j])*np.power(enlogbins[i]/spec_synch['Esyn'][j], x)
				val = (spec_synch['e'][j]/spec_synch['Esyn'][j])*np.power(enlogbins[i]/spec_synch['Esyn'][j], x)
				dNE_sync[i] += val
				dNE[i] += val


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

