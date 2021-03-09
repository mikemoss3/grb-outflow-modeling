"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-02-10

Some useful equations that different parts of the code will need.
"""

import numpy as np 
import cosmologicalconstants as cc
import scipy.integrate as integrate 


def vel(gamma):
	"""
	Method to calculate velocity (m/s) from Lorentz factor (gamma)
	"""

	return cc.c * np.sqrt(1-(1/gamma**2) )

def beta(gamma):
	"""
	Method to calculate the beta factor from the Lorentz factor (beta == v/c)
	"""

	# Approximation
	return 1- 1/(2*gamma**2)
	# Exact
	# return np.sqrt(1- 1/gamma**2)

def rel_vel(gamma_1, gamma_2):
	"""
	Method to calculate the relative velocity of two shells with lorentz factor gamma_1 and gamma_2, where gamma_1 > gamma_2.
	"""

	return (gamma_1**2 - gamma_2**2)/(gamma_1**2 + gamma_2**2)

def lum_dis(z: float):
	""" 
	Caclulate luminosity distance for redshift z
	"""
	lum_dis_Mpc = ((1+z)*cc.c/(cc.H0) ) * integrate.quad(lambda zi: 1/np.sqrt( ((cc.omega_m*np.power(1+zi,3) )+cc.omega_lam) ),0,z)[0]
	lum_dis_cm = lum_dis_Mpc * 3.086e24 # Mpc -> cm
	return lum_dis_cm