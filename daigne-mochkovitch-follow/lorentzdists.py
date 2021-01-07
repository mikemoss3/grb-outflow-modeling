"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2020-08-17


This code defines a number of Lorentz distribution shapes distribution shapes to be used in the simulation of a GRB prompt jet 
made of n consecutive shells 

"""

import numpy as np
import matplotlib.pyplot as plt 

def step(g1=100, g2=400, numshells=5000, mfrac=0.5):
	"""
	Distribute the Lorentz factors of the shells into a step function. 
	The left (right) side of the step will have a step function of g1 (g2), the step between the two occurs 
	at a mass fraction equal to mfrac.

	The mass fraction is defined as the sum of the shell masses up to step divided by the total mass.
	So, the fraction of mass for the region with g1. 

	The first layer is the farthest along in the jet.
	"""

	# number of shells with Lorentz factor g1: 
	n1 =  int(numshells / ( ((1-mfrac)*g2/(mfrac*g1)) + 1 ))
	# number of shells with Lorentz factor g2: 
	n2 = int(numshells - n1)

	# Make array of shells
	shell_arr = np.ndarray(shape=numshells,dtype=[('RADIUS',float),('GAMMA',float),('MASS',float)])

	# Start all shell radii at the initial radii
	shell_arr['RADIUS'] = 10e10

	# Set the Lorentz factors and masses for each section of the step distribution
	shell_arr[0:n1]['GAMMA'] = np.ones(shape=n1)*g1
	shell_arr[0:n1]['MASS'] = np.ones(shape=n1)*mfrac/n1
	shell_arr[n1::]['GAMMA'] = np.ones(shape=n2)*g2
	shell_arr[n1::]['MASS'] = np.ones(shape=n2)*(1-mfrac)/n2
	
	return shell_arr

def load_lorentz_dist(input_filename):
	"""
	Method to upload a Lorentz distribution from a file
	"""

	ld = np.genfromtxt(input_filename,dtype=([('RADIUS',float),('GAMMA',float),('MASS',float)]) )
	return ld

def plot_lorentz_dist(lorentz_arr,title=None):
	"""
	Method to plot the given Lorentz factor distribution
	"""

	# To match paper graphics
	flipped_mass_arr = np.flip(lorentz_arr['MASS'])
	flipped_gamma_arr = np.flip(lorentz_arr['GAMMA'])

	masscum = np.cumsum(flipped_mass_arr)
	massfraccum = masscum/masscum[-1]

	plt.plot(massfraccum,flipped_gamma_arr)
	if title is not None:
		plt.title(title)
	plt.xlabel(r'M/M$_{tot}$')
	plt.ylabel(r'$\Gamma$')


