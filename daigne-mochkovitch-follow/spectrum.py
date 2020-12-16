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

	def add_contribution(self,te,ta,asyn,Beq,gammae,Esyn):
		"""
		Add a contribution to the spectrum
		"""
		# If this is the first contribution to the spectrum, the spectrum must be initialized. 
		if self.spectrum is None:
			self.spectrum = np.array((te,ta,asyn,Beq,gammae,Esyn), dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float)])
		# Otherwise, just append to the already created spectrum
		else:
			self.spectrum = np.append(self.spectrum, np.array((te,ta,asyn,Beq,gammae,Esyn),dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float)])) 

	def load_spectrum(self,input_filename):
		"""
		Method to load a spectrum from a text file
		"""

		self.spectrum = np.genfromtxt(input_filename,dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float)])
		return self.spectrum

	def plot_spectrum(self,num_bins=1000):
		"""
		Method to plot the stored spectrum
		"""

		# Make spectrum bins
		emin = np.min(self.spectrum['Esyn'])
		emax = np.max(self.spectrum['Esyn'])
		enlogbins = np.logspace(np.log10(emin),np.log10(emax),num_bins)

		plt.hist(self.spectrum['Esyn'],bins=enlogbins)
		plt.xscale('log')

