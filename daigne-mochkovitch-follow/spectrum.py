"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2020-09-01


This code defines a Spectrum object class. This class contains all relevant information and definitions of a spectrum created by
GRB prompt emission.

"""

import numpy as np
import matplotlib.pyplot as plt


class Spectrum(object):
	"""
	Spectrum class.
	"""

	def __init__(self):
		"""
		Defines the default parameters of a spectrum.
		"""
		self.spectrum = None 

	def add_contribution(self,te,ta,asyn,Beq,gammae,Esyn,gammar,e,delt):
		"""
		Add a contribution to the spectrum
		"""
		# If this is the first contribution to the spectrum, the spectrum must be initialized. 
		if self.spectrum is None:
			self.spectrum = np.array((te,ta,asyn,Beq,gammae,Esyn,gammar,e,delt), dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float),('gammar',float),('e',float),('delt',float)])
		# Otherwise, just append to the already created spectrum
		else:
			self.spectrum = np.append(self.spectrum, np.array((te,ta,asyn,Beq,gammae,Esyn,gammar,e,delt),dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float),('gammar',float),('e',float),('delt',float)])) 

	def load_spectrum(self,input_filename):
		"""
		Method to load a spectrum from a text file
		"""

		self.spectrum = np.genfromtxt(input_filename,dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float),('gammar',float),('e',float),('delt',float)])
		return self.spectrum

def plot_spectrum(ax,spectrum, nuFnu=True, num_bins=1000,emin=1,emax=1e7):
	"""
	Method to plot the stored spectrum
	"""

	emin = emin # eV 
	emax = emax # eV 
	enlogbins = np.logspace(np.log10(emin),np.log10(emax),num_bins)
	dNE = np.zeros(shape=len(enlogbins))
	for i in range(len(enlogbins)):
		for j in range(len(spectrum['Esyn'])):
			if enlogbins[i] < spectrum['Esyn'][j]:
				x = 2/3
			elif enlogbins[i] > spectrum['Esyn'][j]:
				x = 2.5
			dNE[i] += (spectrum['e'][j]/spectrum['Esyn'][j])*np.power(enlogbins[i]/spectrum['Esyn'][j], -x)


	# For axis labels
	fontsize=14
	fontweight='bold'

	if nuFnu is True:
		ax.scatter(enlogbins,dNE*enlogbins**2)
		ax.set_ylabel(r'E$^2$N(E)',fontsize=fontsize,fontweight=fontweight)
	if nuFnu is False:
		ax.scatter(enlogbins,dNE)
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

