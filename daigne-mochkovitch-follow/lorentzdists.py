"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2020-08-17


This code defines a number of Lorentz distribution shapes distribution shapes to be used in the simulation of a GRB prompt jet 
made of n consecutive shells 

"""

import numpy as np
import matplotlib.pyplot as plt 

from shell import Shell

def step(g1=400, g2=100, numshells=5000, mfrac=0.5):
	"""
	Distribute the Lorentz factors of the shells into a step function. 
	The left (right) side of the step will have a step function of g1 (g2), the step between the two occurs 
	at a mass fraction equal to mfrac.

	The mass fraction is defined as the sum of the shell masses up to step divided by the total mass.
	"""

	# number of shells with Lorentz factor g1: 
	n1 =  int(numshells / ( ((1-mfrac)*g2/(mfrac*g1)) + 1 ))
	# number of shells with Lorentz factor g2: 
	n2 = int(numshells - n1)

	# Make array of shells
	shell_arr = np.ndarray(shape=numshells,dtype=object)

	# Set the Lorentz factors and masses for each section of the step distribution
	# I vectorize the Shell initialization method, this is used instead of looping over the shell array and creating a new Shell instance 
	# for each element. 
	vShell = np.vectorize(Shell)
	shell_arr[0:n1] = vShell(gamma=np.ones(shape=n1)*g1,mass=mfrac/n1)
	shell_arr[n1::] = vShell(gamma=np.ones(shape=n2)*g2,mass=(1-mfrac)/n2)

	
	return shell_arr

def plot_lorentz_dist(lorentz_arr,title=None):
	"""
	Method to plot the given Lorentz factor distribution
	"""

	lorlist = np.zeros(shape=len(lorentz_arr))
	masslist = np.zeros(shape=len(lorentz_arr))
	for i in range(len(lorentz_arr)):
		lorlist[i] = lorentz_arr[i].gamma
		masslist[i] = lorentz_arr[i].mass

	masscum = np.cumsum(masslist)
	massfraccum = masscum/masscum[-1]

	plt.plot(massfraccum,lorlist)
	if title is not None:
		plt.title(title)


