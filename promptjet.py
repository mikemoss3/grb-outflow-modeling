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

	def __init__(self,numshells=5000,dte=0.002,shelldist=sd.step,alpha_e=1/3,alpha_m=1,alpha_b=1/3,ksi=1e-3,mu=1.74,E_dot=1e52,tb=1,theta=0.1,r_open=1e7,r_sat=1e8,eps_th=0.03,sigma=0.1):
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
		self.jet_shells = shelldist(dte=dte,numshells=self.numshells,E_dot=E_dot)

		# Define constants of the jet environment
		self.alpha_e = alpha_e # Fraction of dissipated energy that goes into the electrons 
		self.alpha_m = alpha_m # Between 0.1 - 1, fraction of the dissipated energy which goes into magnetic fluctuation
		self.alpha_b = alpha_b # Fraction of dissipated energy that goes into the magnetic field 
		self.ksi = ksi # Fraction of electrons which are accelerated 
		self.mu = mu # Between 1.5 - 2, index of the fluctuation spectrum 
		self.E_dot = E_dot # erg/s, Injected energy rate
		self.tb = tb # Break out time of the jet (when the first shell breaks out of the ejecta material surrounding the compact object)
		self.theta = theta # radians, Half-opening angle of the jet
		self.r_open = r_open # cm, Opening radius of the jet
		self.r_sat = r_sat # cm, Saturation radius, i.e., where the acceleration of the shells is complete
		self.eps_th = eps_th # Fraction of energy in the outflow in the form of thermal energy 
		self.sigma = sigma # Magnetization of the outflow 

		# Save initial state of shells
		np.savetxt('./sim_results/t0_shells.txt',self.jet_shells)

		# Initialize a spectrum object 
		self.spectrum = sp.Spectrum()

	def jet_evolution(self,jet_shells=None,alpha_e=None,alpha_m=None,alpha_b=None,ksi=None,mu=None,E_dot=None,tb=None,theta=0.1,r_open=1e7,r_sat=1e8,eps_th=0.03,sigma=0.1):
		"""
		Method to emulate the evolution of a jet made of consecutive shells (including the collisions of shells)
		

		"""

		# The following my have been assigned when the class instance was created, but each can be changed for each jet_evolution run (without having to make a new jet)
		if jet_shells is None:
			jet_shells = self.jet_shells
		if alpha_e is None:
			alpha_e = self.alpha_e
		if alpha_m is None:
			alpha_m = self.alpha_m
		if alpha_b is None:
			alpha_b = self.alpha_b
		if ksi is None:
			ksi = self.ksi
		if mu is None:
			mu = self.mu
		if E_dot is None:
			E_dot = self.E_dot
		if tb is None:
			tb = self.tb	
		if theta is None:
			theta = self.theta
		if r_open is None:
			r_open = self.r_open
		if r_sat is None:
			r_sat = self.r_sat
		if eps_th is None:
			eps_th = self.eps_th
		if sigma is None:
			sigma = self.sigma

		true_t = tb # Initialize the time in the reference frame of the jet. 


		# Open file to record collision dynamics
		file_shell_dyn = open('./sim_results/shell_dyn.txt','w')
		file_shell_dyn.write("TimeObs (s) \tTimeJet (s) \tRcoll (cm) \tDSGamma \tDSMass (g) \tUPGamma \tUPMass (g) \tGammaf \ttau\n")
		file_shell_dyn.close()
		file_shell_dyn = open('./sim_results/shell_dyn.txt','a')
		# Initialize flags, used for creating plots:
		t3e4_flag=False
		t2e5_flag=False
		t5e5_flag=False


		### Calculate Thermal considerations
		# Thermal considerations do not depend on shell collisions, but only on the energy in a particular shell as it passes the photosphere

		# Calculate the photospheric radius for each jet shell, Equation 9 of Hascoet 2013
		r_phot = (0.2*E_dot / ( 8*np.pi*(cc.c**3) *(jet_shells['GAMMA']**3)) ) 
		# r_phot = 2.9*10**13 * (0.2*(E_dot/1e53) / ( (1+sigma)*(jet_shells['GAMMA']/100)) )

		# Times when each shell will cross the photoshpere
		t_cross_phot = (r_phot - jet_shells['RADIUS']) / ( cc.c*beta(jet_shells['GAMMA']) )
		# Observer time 
		t_obs = t_cross_phot-r_phot/cc.c
		# Calculate useful constant for next calculations
		Phi = np.power(theta,-2/3)*np.power(r_phot,-2/3)*np.power(r_open,2/3)*np.power(jet_shells['GAMMA'],2/3)

		# Temperature at photosphere
		a = 7.566 * 1e-15 # erg cm^-3 K^-4, Radiation constant
		T0 = np.power(E_dot*eps_th*theta**2 / (16*np.pi**2 * a * cc.c * r_open**2),1/4) # K 
		# T0 = (2/3)*np.power(eps_th,1/4)*np.power(theta/0.1,1/2)*np.power(E_dot/1e53,1/4)*np.power(r_open/1e7,-1/2) # MeV
		T_phot = T0 * Phi

		# Luminosity at photosphere 
		L_phot = (eps_th*theta**2*E_dot/(16*np.pi)) * Phi

		for i in range(len(jet_shells)):
			self.spectrum.add_therm_contribution(t_cross_phot[i],t_obs[i],T_phot[i],L_phot[i])
		


		### Simulate Internal Shocks
		# Flag to mark if all the shells in the jet are ordered (no more collisions occur at this point )
		ord_lorentz = False

		# Check if all the active shells have ordered Lorentz factors
		# If true, no more collisions will be possible at this point
		if is_sorted(jet_shells[jet_shells['STATUS']==1]['GAMMA'] ) == True:
			ord_lorentz = True
			print("Jet shells are already ordered.")


		while ord_lorentz == False:
			
			# Calculate the time of collision between all adjacent (and active) shells
			t_coll_lowest = 1e99 # Place holder time until next collision
			
			# Grab all the indices for the active shells
			active_shell_inds = np.where(jet_shells["STATUS"]==1)[0]
			# For each active shell, find the next active shell
			# I don't want to look at the very final active shell, because there will be no shells after it. 
			for i in range(len(active_shell_inds)-1):

				# The current active shell is the downstream shell
				ind_ds_shell = active_shell_inds[i]
				# The next active shell will be the up stream shell
				ind_us_shell = active_shell_inds[i+1]

				# If the down stream shell must have a Lorentz factor smaller than the upstream shell or they will never collide
				if jet_shells['GAMMA'][ind_ds_shell] < jet_shells['GAMMA'][ind_us_shell]: 

					# The down stream shell is farther out, but slower
					# The up stream shell is closer to the central engine out, but faster
					tmp_ds_r = jet_shells['RADIUS'][ind_ds_shell]
					tmp_ds_b = beta(jet_shells['GAMMA'][ind_ds_shell])
					tmp_us_r = jet_shells['RADIUS'][ind_us_shell]
					tmp_us_b = beta(jet_shells['GAMMA'][ind_us_shell])

					# Calculate the time until collision between these two shells
					t_coll_tmp = (tmp_ds_r - tmp_us_r) / (cc.c * (tmp_us_b - tmp_ds_b) )


					# Check if this time is shorter than our previous time until next collision
					if t_coll_tmp < t_coll_lowest:
						t_coll_lowest = t_coll_tmp

						# Record which shells these were
						ind_s_ds = ind_ds_shell
						ind_s_us = ind_us_shell

			true_t += t_coll_lowest

			# These are the Lorentz factors, masses, and radii of the two colliding shells
			# Down stream shell 
			shell_ds_g = jet_shells[ind_s_ds]['GAMMA']
			shell_ds_m = jet_shells[ind_s_ds]['MASS']
			# shell_ds_r = jet_shells[ind_s_ds]['RADIUS']
			# Up stream shell
			shell_us_g = jet_shells[ind_s_us]['GAMMA']
			shell_us_m = jet_shells[ind_s_us]['MASS']
			# shell_us_r = jet_shells[ind_s_us]['RADIUS']

			# Radius at which collision occurred
			rad_coll = jet_shells['RADIUS'][ind_s_ds] + cc.c*beta(shell_ds_g)*t_coll_lowest 
			# rad_coll = jet_shells['RADIUS'][ind_s_us] + cc.c*beta(shell_us_g)*t_coll_lowest # Alternative way to calculate
			# rad_coll = (shell_ds_r - shell_us_r) * ( beta(shell_us_g) / (beta(shell_us_g) - beta(shell_ds_g) )) # Alternative way to calculate, but suffers from numerical instabilities

			### Calculate Contribution to Spectrum ### 

			## Synchrotron considerations:
			gamma_int = 0.5*(np.sqrt(shell_ds_g/shell_us_g)+np.sqrt(shell_us_g/shell_ds_g)) # Lorentz factor for internal motion in shocked material
			eps = (gamma_int-1)*cc.mp*cc.c**2 # erg, Average proton factor from the collision of two shells

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
				
				# Calculate the Compton parameter
				Q_IC = (-1 + np.sqrt(1-4*c))/2
				alpha_ic = (Q_IC/w)/(1+(Q_IC/w)) # Fraction of energy that goes into Inverse Compton electrons
			else: 
				# w < 1 
				# Thomson
				# c = -0.2 *(E_dot/cc.c**2/gamma_bar**2)*t_syn*gamma_e**2 / (4*np.pi*(cc.c*t_var*gamma_bar**2)**2)
				
				# Calculate the Compton parameter
				Q_IC = (-1 + np.sqrt(1-4*c))/2
				alpha_ic = Q_IC/(1+Q_IC) # Fraction of energy that goes into Inverse Compton electrons

			alpha_syn = 1 - alpha_ic # Fraction of energy that goes into Synchrotron electrons 
			

			# Other efficiency quantities to record:
			f_dyn = (np.sqrt(shell_us_g/shell_ds_g)-1)**2 / ((shell_us_g/shell_ds_g)+1) # Dynamic efficiency, Equation 2 from Bosnjak, Daigne, Dubus 2009 
			# f_dyn = 1 - (shell_ds_m+shell_us_m)*gamma_comb/( shell_ds_m*shell_ds_g + shell_us_m*shell_us_g) # Dynamic efficiency, Equation 3 from Kobayashi, Piran, and Sari 1997 
			# f_dyn = 1

			# Energy dissipated in the shock, assumes that the energy dissipated occurs when upstream shell sweeps up a mass equal to its own (or sweeps up the entire downstream shell mass)
			en_diss =  np.min(np.array([shell_ds_m,shell_us_m])) * cc.c**2 * (shell_ds_g + shell_us_g - 2*gamma_r)
			# en_diss *= (alpha_syn * alpha_e * f_dyn * 3e-3)

			# Variability timescale 
			delt = rad_coll/(2*cc.c*gamma_r**2) # sec, Dynamical time scale of the shell

			# Observer arrival time
			ta = true_t - (rad_coll/cc.c)

			# Calculate optical depth 
			if ind_s_ds == 0:
				tau = 0
			else:
				tau = 0.2 * np.sum( jet_shells['MASS'][jet_shells['STATUS']==1][:ind_s_ds]/(4 * np.pi * jet_shells['RADIUS'][jet_shells['STATUS']==1][:ind_s_ds]**2) )

			# If the emission is efficient, add the contribution
			# E.g., If the emission time is less than the shell expansion (i.e., the dynamical scale of the shell)
			t_syn = 6*np.power(gamma_e/100,-1)*np.power(Beq/1000,-2) # sec, synchrotron time-scale
			if t_syn < ((1+Q_IC)*rad_coll/cc.c/gamma_bar):
				# The relative velocity between the shells must be greater than the local sound speed
				if rel_vel(shell_us_g,shell_ds_g) > 0.1:
					# The wind must be transparent in order to produce observable emission 
					if tau < 1:	
						# Add contribution to the spectrum 
						self.spectrum.add_synch_contribution(te=true_t,ta=ta,asyn=alpha_syn,Beq=Beq,gammae=gamma_e,Esyn=E_syn_ev, gammar=gamma_r, e=en_diss, delt = delt)	


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
			# Set the radius of the merged shells.
			jet_shells[ind_s_ds]['RADIUS'] = rad_coll
			# Mass of the combined shell (add the mass of the other shell)
			jet_shells[ind_s_ds]['MASS'] += jet_shells[ind_s_us]['MASS']




			# Record collision dynamics for each collision
			save_coll_arr = np.array([ta,true_t,rad_coll,shell_ds_g,shell_ds_m,shell_us_g,shell_us_m,gamma_comb,tau])
			[file_shell_dyn.write("{:1.5e}\t".format(save_coll_arr[i])) for i in range(len(save_coll_arr))]
			file_shell_dyn.write("\n")

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


			# Check if all the active shells have ordered Lorentz factors
			# If true, no more collisions will be possible at this point
			if is_sorted(jet_shells[jet_shells['STATUS']==1]['GAMMA'] ) == True:
				ord_lorentz = True
				np.savetxt('./sim_results/ordlor_shells.txt',jet_shells)
				np.savetxt('./sim_results/ordlor_spectrum_synch.txt',self.spectrum.spec_synch)
				np.savetxt('./sim_results/ordlor_spectrum_therm.txt',self.spectrum.spec_therm)
				print("At time t={} s, all shells have been launched and Lorentz factors are ordered.".format(true_t))

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





