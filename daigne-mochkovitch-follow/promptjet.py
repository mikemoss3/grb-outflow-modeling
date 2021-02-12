"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2020-08-17


This code defines a prompt-jet object class. The jet is comprised of many shells (or layers), each shell has a different mass, 
Lorentz factor, and distance from the central engine. Depending on respective parameters of different shell, shells may collide, 
producing internal shocks. 

This code aims to recreate the results of Daigne and Mochkovitch 1998.  

"""

import numpy as np
import matplotlib.pyplot as plt 

from utils import *
import shelldists as sd
import spectrum as sp
import cosmologicalconstants as cc


class PromptJet(object):
	"""
	Prompt Jet class.
	"""

	def __init__(self,numshells=5000,dte=0.002,shelldist=sd.step):
		"""
		Defines the default parameters of the jet.

		Parameters:
		numshells = number of shells to be launched by the central engine
		dte = time between shell launches, this can be specific by a single float to apply a constant time step through out the jet evolution or can be a array of the shell emission times
		shelldist = the shape of the lorentz distribution of the shells, this should be on of the functions described in lorentzdists.py 

		"""
		# Number of shells
		self.numshells = numshells
		
		# Check if a single time step was given or a list of launch times
		# If a list of launch times was given was given
		if hasattr(dte,"__len__"):
			# Check if the list is the same size as the number of shells
			if len(dte) != numshells:
				print("The list of shell launch times must be the same size as the number of shells.")

			# Total duration of central engine activity (total time to produce numshells), in seconds
			# This will just be the final jet launch time
			self.tw = dte[-1] 
		# Else if a single constant difference between launch time
		else:
			self.dte = dte
			
			# Total duration of central engine activity (total time to produce numshells), in seconds 
			self.tw = self.dte*(self.numshells+1)


		# Make the list of shells with a Lorentz distributions (this is the distribution at t = tw)
		self.jet_shells = shelldist(dte=dte,numshells=self.numshells)
		np.savetxt('./sim_results/t0_shells.txt',self.jet_shells)


		# Initialize a spectrum object 
		self.spectrum = sp.Spectrum()

	def jet_evolution(self,tb=1,jet_shells=None):
		"""
		Method to emulate the evolution of a jet made of consecutive shells (including the collisions of shells)
		
		Params:
		tb = break out time, time that the first shell breaks out of the photosphere

		"""

		if jet_shells is None:
			jet_shells = self.jet_shells

		# Time until next shell is launched
		t_shell_launch = self.dte
		# Next shell to be launched
		activate_ind = 1

		# Keep track of global time
		true_t = tb

		# Initialize test flags:
		t3e4_flag=False
		t2e5_flag=False
		t5e5_flag=False

		# Flag to mark if all the shells in the jet are ordered (no more collisions occur at this point )
		ord_lorentz = False
		while ord_lorentz == False:
			
			# Calculate the time of collision between all adjacent (and active) shells
			t_coll_lowest = 1e99 # Place holder time until next collision
			for i in range(len(jet_shells)-1):
				# We only need to check the active shells 
				if jet_shells['STATUS'][i] == 1:

					# Find the closest up stream shell that is active
					found_nxt = False
					ind_nxt_shell=0
					tmp_ind = i
					while found_nxt == False:
						tmp_ind+=1
						if jet_shells['STATUS'][tmp_ind]==1:
							found_nxt = True
							ind_nxt_shell = tmp_ind

					# If the down stream shell must have a Lorentz factor smaller than the upstream shell or they will never collide
					if jet_shells['GAMMA'][i] < jet_shells['GAMMA'][ind_nxt_shell]: 
						# The down stream shell is farther out, but slower
						# The up stream shell is closer to the central engine out, but faster
						tmp_ds_r = jet_shells['RADIUS'][i]
						tmp_ds_b = beta(jet_shells['GAMMA'][i])
						tmp_us_r = jet_shells['RADIUS'][ind_nxt_shell]
						tmp_us_b = beta(jet_shells['GAMMA'][ind_nxt_shell])

						# Calculate the time until collision between these two shells
						t_coll_i = (tmp_ds_r - tmp_us_r) / (cc.c * (tmp_us_b - tmp_ds_b) )
						
						# Check if this time is shorter than our previous time until next collision
						if t_coll_i < t_coll_lowest:
							t_coll_lowest = t_coll_i

							# Record which shells these were
							ind_s_ds = i
							ind_s_us = ind_nxt_shell

			true_t += t_coll_lowest

			# These are the Lorentz factors, masses, and radii of the two colliding shells
			# Down stream shell 
			shell_ds_g = jet_shells[ind_s_ds]['GAMMA']
			shell_ds_m = jet_shells[ind_s_ds]['MASS']
			shell_ds_r = jet_shells[ind_s_ds]['RADIUS']
			# Up stream shell
			shell_us_g = jet_shells[ind_s_us]['GAMMA']
			shell_us_m = jet_shells[ind_s_us]['MASS']
			shell_us_r = jet_shells[ind_s_us]['RADIUS']

			# Radius at which collision occurred
			rad_coll = jet_shells['RADIUS'][ind_s_ds] + cc.c*beta(shell_ds_g)*t_coll_lowest 
			# rad_coll = jet_shells['RADIUS'][ind_s_us] + cc.c*beta(shell_us_g)*t_coll_lowest # Alternative way to calculate
			# rad_coll = (shell_ds_r - shell_us_r) * ( beta(shell_us_g) / (beta(shell_us_g) - beta(shell_ds_g) )) # Alternative way to calculate, but suffers from numerical instabilities

			# print(t_coll_lowest)
			# print(ind_s_ds, ind_s_us)
			# print(jet_shells['RADIUS'][ind_s_ds], jet_shells['RADIUS'][ind_s_us])
			# print(jet_shells['GAMMA'][ind_s_ds], jet_shells['GAMMA'][ind_s_us])
			# print(rad_coll)
			# print('\n')
			

			# Now that we have found the correct next collision, find contribution to emission:
			### Calculate Contribution to Spectrum ### 

			gamma_int = 0.5*(np.sqrt(shell_ds_g/shell_us_g)+np.sqrt(shell_us_g/shell_ds_g)) # Lorentz factor for internal motion in shocked material
			eps = (gamma_int-1)*cc.mp*cc.c**2 # erg, Average proton factor from the collision of two shells

			# Calculate the characteristic electron Lorentz factor for different assumption: 
			# A) Synchrotron emission from accelerated electrons 
			# B) the electron Lorentz factor is large enough to directly produce Synchroton radiation

			# Define constants assumed by model and used by both assumptions
			alpha_e = 1/3 # Fraction of dissipated energy that goes into the electrons 
			alpha_m = 1 # Between 0.1 - 1, fraction of the dissipated energy which goes into magnetic fluctuation
			alpha_b = 1/3 # Fraction of dissipated energy that goes into the magnetic field 
			ksi = 1e-3 # Fraction of electrons which are accelerated 
			mu = 1.75 # Between 1.5 - 2, index of the fluctuation spectrum 
			E_dot = 1e52 # erg/s, Injected energy rate

			# Characteristic Lorentz factor of accelerated electrons					
			# gamma_e = alpha_e *eps / (cc.me* cc.c**2) # Average electron lorentz factor, for electrons in equipartition with protons
			# gamma_e = np.power((alpha_m/ksi)*(eps/cc.me/cc.c**2),1/(3-mu)) # Average electron lorentz factor, electrons directly produce synchrotron radiation (via turbulent magnetic fields)
			gamma_e = 1e4 # Average electron lorentz factor, constant
			gamma_r = np.sqrt(shell_ds_g * shell_us_g) # Approximate resulting Lorentz factor from the shell collision
			gamma_bar = (shell_ds_g + shell_us_g)/2 # Average Lorentz factor of jet shells 


			n = E_dot / (4*np.pi * rad_coll**2 * gamma_bar**2 * cc.mp * cc.c**3) # 1/cm^3, Comoving proton number density
			# n = E_dot / (4*np.pi * rad_coll**2 * (gamma_bar*gamma_r) * cc.mp * cc.c**3) # 1/cm^3, Comoving proton number density, new implementation of Frederic
			Beq = np.sqrt(8*np.pi*alpha_b*n*eps) # (erg/cm^3)^1/2, Magnetic field density, assuming B = Beq  

			# Synchrotron energy emitted by accelerated electron
			E_syn_ev = 50*(gamma_r/300)*(Beq/1e3)*(gamma_e/100)**2 # eV, Synchrotron energy in the rest frame
			E_syn = E_syn_ev * 1.60218e-12 # convert to erg
			E_0syn = E_syn/gamma_r # erg, Synchrotron energy in the comoving frame

			# Check if Klein-Nishina limit or Thompson (i.e., w>>1 or w<<1)
			w = gamma_e*E_0syn/(cc.me*cc.c**2) # Critical value between Thomson and Klein-Nishina
			# w_alt = 33*(Beq/1e3)*(gamma_e/1e4)**3 # Alternative relation fro w
			# print('w={} , wp={}'.format(w, w_alt) )
			
			t_syn = 6*np.power(gamma_e/100,-1)*np.power(Beq/1e3,-2) # sec, synchrotron time-scale

			# Q_IC == tau_star * gamma_e**2 == Y == the Compton Parameter
			# Write equation 22 in a slightly simplified way
			c = -(3/2)*(0.2/np.pi/rad_coll**2)*(E_dot/gamma_bar**2 / cc.c**2)*(gamma_e)*np.power(Beq/1000,-2)*100
			if w >= 1:
				# Klein-Nishina
				# Calculate the left hand side of equation 29
				# Notice the negative sign to move it to the right hand side (to use the quadratic formula)
				# c = -8*1e-4*(1+np.log(2*w))*(E_dot/1e52)*np.power(t_var/1,-2)*np.power(gamma_bar/300,-6)*np.power(Beq/1e3,-4)*np.power(gamma_e/1e4,-5)
				# Q_ICw = (-1 + np.sqrt(1-4*c))/2
				# alpha_ic = Q_ICw/(1+Q_ICw) # Fraction of energy that goes into Inverse Compton electrons
				
				Q_IC = (-1 + np.sqrt(1-4*c))/2
				alpha_ic = (Q_IC/w)/(1+(Q_IC/w)) # Fraction of energy that goes into Inverse Compton electrons
			else: 
				# w < 1 
				# Thomson
				# c = -0.2 *(E_dot/cc.c**2/gamma_bar**2)*t_syn*gamma_e**2 / (4*np.pi*(cc.c*t_var*gamma_bar**2)**2)
				Q_IC = (-1 + np.sqrt(1-4*c))/2
				alpha_ic = Q_IC/(1+Q_IC) # Fraction of energy that goes into Inverse Compton electrons

			alpha_syn = 1 - alpha_ic # Fraction of energy that goes into Synchrotron electrons 
			

			# Other efficiency quantities to record:
			# f_dyn = (np.sqrt(shell_us_g/shell_ds_g)-1)**2 / ((shell_us_g/shell_ds_g)+1), # dynamic efficiency
			# f_dyn = 1 - (shell_ds_m+shell_us_m)*gamma_comb/( shell_ds_m*shell_ds_g + shell_us_m*shell_us_g) # Equation 3 from Kobayashi, Piran, and Sari 1997 
			# f_dyn = 1

			# Energy dissipated in the shock, assumes that the energy dissipated occurs when upstream shell sweeps up a mass equal to its own (or sweeps up the entire downstream shell mass)
			en_diss = np.min(np.array([shell_ds_m,shell_us_m])) * cc.c**2 * (shell_ds_g + shell_us_g - 2*gamma_r)

			# Variability timescale 
			delt = rad_coll/(2*cc.c*gamma_r**2) # sec, Dynamical time scale of the shell

			# If the emission is efficient, add the contribution
			# E.g., If the emission time is less than the shell expansion (i.e., the dynamical scale of the shell)
			t_syn = 6*np.power(gamma_e/100,-1)*np.power(Beq/1000,-2) # sec, synchrotron time-scale
			if t_syn < ((1+Q_IC)*rad_coll/cc.c/gamma_bar):
				# The relative velocity between the shells must be greater than the local sound speed
				if rel_vel(shell_us_g,shell_ds_g) > 0.1:
					ta = true_t - (rad_coll/cc.c)
					
					# The wind most also be transparent to the produced
					if ind_s_ds == 0:
						tau = 0
					else:
						M_sol = 2*1e33 # g, mass of the sun
						M_typical = 0.001 * M_sol # g, typical mass in prompt jet
						M_typical = 1
						tau = 0.2 * np.sum( M_typical*jet_shells['MASS'][jet_shells['STATUS']==1][:ind_s_ds]/(4 * np.pi * jet_shells['RADIUS'][jet_shells['STATUS']==1][:ind_s_ds]**2) )
					# print(tau)
					if tau < 1:	
						# Add contribution to the spectrum 
						self.spectrum.add_contribution(te=true_t,ta=ta,asyn=alpha_syn,Beq=Beq,gammae=gamma_e,Esyn=E_syn_ev, gammar=gamma_r, e=en_diss, delt = delt)
					

			# De-active the up-stream jet shell
			jet_shells[ind_s_us]['STATUS'] = 0 

			# Move shells forward 
			for i in range(len(jet_shells)):
				# But only if the shells are active and launched
				if jet_shells[i]['STATUS'] == 1:
					jet_shells[i]['RADIUS'] += cc.c*beta(jet_shells[i]['GAMMA'])*t_coll_lowest


			# Final Lorentz factor of the combined shell (after complete redistribution of momenta)
			gamma_comb = shell_coll_gamma(shell_ds_g, shell_us_g, shell_ds_m, shell_us_m)

			# Set the Lorentz factor of the merged shells.
			jet_shells[ind_s_ds]['GAMMA'] = gamma_comb
			# Mass of the combined shell (add the mass of the other shell)
			jet_shells[ind_s_ds]['MASS'] += jet_shells[ind_s_us]['MASS']




			# Record shell positions at specific times (to match D&M 1998)
			if true_t > 3e4 and t3e4_flag==False:
				np.savetxt('./sim_results/t3e4_shells.txt',jet_shells)
				t3e4_flag=True 
			if true_t > 2e5 and t2e5_flag==False:
				np.savetxt('./sim_results/t2e5_shells.txt',jet_shells)
				t2e5_flag=True 
			if true_t > 5e5 and t5e5_flag==False:
				np.savetxt('./sim_results/t5e5_shells.txt',jet_shells)
				t5e5_flag=True 

			# Check if all the shells have been launched
			# and 
			# Check if all the active shells have an ordered Lorentz factors are ordered
			# No more collisions will be possible at this point
			if is_sorted(jet_shells[jet_shells['STATUS']==1]['GAMMA'] ) == True:
				ord_lorentz = True
				np.savetxt('./sim_results/ordlor_shells.txt',jet_shells)
				np.savetxt('./sim_results/ordlor_spectrum.txt',self.spectrum.spectrum)
				print("At time t={} s, all shells have been launched and Lorentz factors are ordered.".format(true_t))


def calc_t_coll(shell_1,shell_2):
	"""
	Function to calculate the time of collision between two shells.
	Params:
	shell_1 = (downsteam shell) contains the mass, radius, Lorentz factor, and emission time of shell 1 
	shell_2 = (upstream shell) contains the mass, radius, Lorentz factor, and emission time of shell 2 
	# We assume (Radius_shell_1 > Radius_shell_2) and (Gamma_shell_1 < Gamma_shell_2)
	"""

	# Current radius of shell 1 and shell 2 
	r10 = shell_1['RADIUS']
	r20 = shell_2['RADIUS']
	# Lorentz factors of shell 1 and shell 2
	g1 = shell_1['GAMMA']
	g2 = shell_2['GAMMA']

	# Calculate the collision time 
	t_coll = (2*g2**2)*(r10-r20)/(cc.c*((g2/g1)**2 - 1))
	# t_coll = (r10-r20)/cc.c/(beta(g2)-beta(g1)) # Alternative relation from Kobayashi, Piran, Sari 1997


	# Check for nan value
	if t_coll != t_coll:
		t_coll = 1e99

	return t_coll

def shell_coll_gamma(s1g,s2g,s1m,s2m):
	"""
	Find the resulting Lorentz factor of two collided shells from the Lorentz factor and Mass of shell 1 and 2
	"""

	# The approximate Lorentz factor of the combined layer after complete redistribution of momentum and energy
	gamma_comb = np.sqrt( s1g*s2g * (s1m*s1g + s2m*s2g) / (s1m*s2g + s2m*s1g)  )

	return gamma_comb



def is_sorted(x):
	"""
	Method to test whether the array, x, is sorted
	Sorted in this context means x[i] > x[i+1]
	"""
	return (np.diff(x)<=0).all()





