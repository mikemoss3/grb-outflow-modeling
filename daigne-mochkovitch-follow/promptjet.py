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
	
			# Move through each shell starting at the second farthest out (2nd most downstream) shell and moving upstream
			## We can skip the farthest out shell because it doesn't have any shells to collide with
			# Calculate time of collisions between shell_i and all down stream
			# Find the shortest time of collision 
			# Analytically calculate collision 
			# Calculate contribution to spectra 
			# Move all shells forward by the time step found (i.e., the shortest time of collision) 

			# First check if the shells have ordered Lorentz factors, no collisions will be possible
			if is_sorted(jet_shells[jet_shells['STATUS']==1]['GAMMA'] ) == True:
				# if so, but there are still more shells to be launched, then launch the next shell.
				if jet_shells[-1]['STATUS'] == 2:
					# Move shells forward by the amount of time it took for the new shell to be launched
					for i in range(len(jet_shells)):
						# But only if the shells are active
						if jet_shells[i]['STATUS'] == 1:
							jet_shells[i]['RADIUS'] += vel(jet_shells[i]['GAMMA'])*t_shell_launch


					# Launch (activate) next shell
					jet_shells[activate_ind]['STATUS'] = 1

					# Keep track of the next shell to be launched.
					activate_ind +=1

				# If so, but there are no more shells, then end the simulation
				else:
					ord_lorentz = True
					print("All shells have been launched and Lorentz factors are ordered.")
			else: 
				# Find the next collision time
				t_coll_lowest = 1e99
				for i in range(1,len(jet_shells)):					
					# Only evaluate the shells that are active and launched 
					if jet_shells[i]['STATUS'] == 1:
						# Calculate the time of collision between this shell and the next shell down stream of it
							
						# Find the next active down stream shell
						# argmax quickly finds the first element that is equal to 'True', which is why we flip the Active array and then subtract the index from the length
						flipped_ind = np.argmax(np.flip(jet_shells[0:i]['STATUS'])==1)
						# By the way, if argmax does not find any match, it returns 0 (i.e., False), so I have to check if argmax is trying to return index=0 or False.
						# Honestly, I cannot believe it was coded like that. 
						if np.flip(jet_shells[0:i]['STATUS'])[flipped_ind] == 1:

							# The index of the next active shell 
							# The -1 is needed because of array indexing starting at 0 when flipped. 
							nxt_act_shell_ind = len(jet_shells[0:i]) - 1 - flipped_ind
							# If the next shell has a Lorentz factor greater than or equal to, the two shells cannot collide
							if jet_shells[nxt_act_shell_ind]['GAMMA'] < jet_shells[i]['GAMMA']:
								# The total time of collision equals the time it will the take the two shells to collide plus the breakout time
								t_coll = calc_t_coll(jet_shells[nxt_act_shell_ind],jet_shells[i])


								# check if this is the lowest collision time so far
								if t_coll < t_coll_lowest:
									# Record the current lowest collision time
									t_coll_lowest = t_coll
									# Keep track of which shells these were
									ind_shell_1 = nxt_act_shell_ind
									ind_shell_2 = i

				# We have found when the next collision will be, but we should now check if the next collision will be before or after a new shells is launched
				
				# If the collision happens after a new shell is launched and there are still more shells to be launched
				# Then launch (activate) that a new shell has been launched (it will be positioned at the photosphere)
				# and move things forward by the time it took for the new shell to launch.
				# Check if this new shell will collide with the shell downstream of it before the other collision time we found. 
				if t_coll_lowest > t_shell_launch and jet_shells[-1]['STATUS'] == 2:
					# Then the change of time will be the time between shell launched:
					delta_t = t_shell_launch
					true_t += delta_t 


					# Move shells forward by the amount of time it took for the new shell to be launched
					for i in range(len(jet_shells)):
						# But only if the shells are active
						if jet_shells[i]['STATUS'] == 1:
							jet_shells[i]['RADIUS'] += vel(jet_shells[i]['GAMMA'])*delta_t

					# Launch (activate) next shell
					jet_shells[activate_ind]['STATUS'] = 1

					# Keep track of the next shell to be launched.
					activate_ind +=1 
				else:
					if t_coll_lowest < t_shell_launch: 
						# The next time a shell is launched is now reduced by the 
						t_shell_launch -= t_coll_lowest			

					# The time increment is the lowest collision time
					delta_t = t_coll_lowest
					true_t += delta_t 


					# Now that we have found the correct next collision:
					# Proceed with evaluating the collision and moving all launched shells forward 			
					# To 'combine' the two shells, we delete the old shell by turning it off, and combine their Lorentz factors and masses
					jet_shells[ind_shell_1]['STATUS'] = 0

					# Move shells forward 
					for i in range(len(jet_shells)):
						# But only if the shells are active and launched
						if jet_shells[i]['STATUS'] == 1:
							jet_shells[i]['RADIUS'] += vel(jet_shells[i]['GAMMA'])*t_coll_lowest


					# Calculate the Lorentz factor of collided shell
					s1g = jet_shells[ind_shell_1]['GAMMA']
					s2g = jet_shells[ind_shell_2]['GAMMA']
					s1m = jet_shells[ind_shell_1]['MASS']
					s2m = jet_shells[ind_shell_2]['MASS']
					gamma_comb = shell_coll_gamma(s1g, s2g, s1m, s2m)
					jet_shells[ind_shell_2]['GAMMA'] = gamma_comb

					# Mass of the new shell (add the mass of the other shell)
					jet_shells[ind_shell_2]['MASS'] += jet_shells[ind_shell_1]['MASS']

					### Calculate Contribution to Spectrum ### 
					gamma_int = 0.5*(np.sqrt(jet_shells['GAMMA'][ind_shell_1]/jet_shells['GAMMA'][ind_shell_2])+np.sqrt(jet_shells['GAMMA'][ind_shell_2]/jet_shells['GAMMA'][ind_shell_1])) # Lorentz factor for internal motion in shocked material
					eps = (gamma_int-1)*cc.mp*cc.c**2 # # Average proton factor from the collision of two shells

					# Calculate the characteristic electron Lorentz factor for different assumption: 
					# A) protons and electrons are in equipartition
					# B) the electron Lorentz factor is large enough to directly produce Synchroton radiation

					# Start assume A) 
					alpha_e = 1/3 # Fraction of dissipated energy that goes into the electrons 
					# gamma_e = alpha_e *eps / (cc.me* cc.c**2) # Average electron lorentz factor

					# End assume A)
					
					# Start assume B) 
					# Define assumed constants of the model
					alpha_m = 0.5 # Between 0.1 - 1, fraction of the dissipated energy which goes into magnetic fluctuation
					ksi = 1e-3 # Fraction of electrons which are accelerated 
					mu = 1.75 # Between 1.5 - 2, index of the fluctuation spectrum 
					E_dot = 1e52/(4*np.pi) # erg/s, Injected energy rate

					# Calculate necessary parameters to find calculate observable parameters 
					gamma_e = np.power((alpha_m/ksi)*(eps/cc.me/cc.c**2),1/(3-mu)) # Average electron lorentz factor

					n = E_dot / (4*np.pi * jet_shells[ind_shell_2]['RADIUS']**2 * gamma_comb**2 * cc.mp * cc.c**3) # Comoving proton number density
					alpha_b = 1/3 # Fraction of dissipated energy that goes into the magnetic field 
					B = np.sqrt(8*np.pi*alpha_b*n*eps) # Magnetic field strength, assuming B = Beq  
					gamma_r = np.sqrt(jet_shells['GAMMA'][ind_shell_1]*jet_shells['GAMMA'][ind_shell_2]) # Approximate resulting Lorentz factor from the shell collision
					E_syn = 50*(gamma_r/300)*(B/1000)*(gamma_e/100)**2 # eV, Synchrotron energy in the rest frame

					E_0syn = E_syn/gamma_r # Synchrotron energy in the comoving frame
					w = gamma_e*E_0syn/(cc.me*cc.c**2) # Critical value between Thomson and Klein-Nishina
					t_var = jet_shells[ind_shell_2]['RADIUS']/(cc.c*gamma_r**2) # Dynamical time scale of the shell

					# Calculate the optical depth, check if Klein-Nishina limit or Thompson (i.e., w>>1 or w<<1)
					# Q_IC == tau_star * gamma_e**2 == Y == the Compton Parameter
					if w > 10:
						# Klein-Nishina
						# Calculate the left hand side of equation 29
						# Notice the negative sign to move it to the right hand side (to use the quadratic formula)
						c = -8*1e-4*(1+np.log(2*w))*(E_dot/1e52)*np.power(t_var/1,-2)*np.power(gamma_comb/300,-6)*np.power(B/1000,-4)*np.power(gamma_e/1e4,-5)
						Q_IC = w * (-1 + np.sqrt(1-4*c))/2
						alpha_ic = (Q_IC/w)/(1+(Q_IC/w)) # Fraction of energy that goes into Inverse Compton electrons 
					else: 
						# Thomson
						c = - 0.3 *alpha_e/alpha_b
						Q_IC = (-1 + np.sqrt(1-4*c))/2
						alpha_ic = Q_IC/(1+Q_IC) 


					alpha_syn = 1 - alpha_ic # Fraction of energy that goes into Synchrotron electrons 
					# End assume B)
					
					# If the emission is efficient, add the contribution
					# If the emission time is less than the shell expansion (i.e., the dynamical scale of the shell)
					t_syn = 6*np.power(gamma_e/100,-1)*np.power(B/1000,-2) # sec, synchrotron time-scale
					if t_syn < (gamma_comb*t_var*(1+Q_IC)):
						ta = true_t - (gamma_comb/cc.c)
						self.spectrum.add_contribution(te=true_t,ta=ta,asyn=alpha_syn,Beq=B,gammae=gamma_e,Esyn=E_syn)



				# Testing
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
				if jet_shells[-1]['STATUS'] == 1 and is_sorted(jet_shells[jet_shells['STATUS']==1]['GAMMA'] ) == True:
					ord_lorentz = True
					np.savetxt('./sim_results/ordlor_shells.txt',jet_shells)
					np.savetxt('./sim_results/ordlor_spectrum.txt',self.spectrum.spectrum)
					print("At time t={} s, all shells have been launched and Lorentz factors are ordered.".format(true_t))


def calc_t_coll(shell_1,shell_2):
	"""
	Function to calculate the time of collision between two shells.
	Params:
	shell_1 = contains the mass, radius, Lorentz factor, and emission time of shell 1 
	shell_2 = contains the mass, radius, Lorentz factor, and emission time of shell 2 
	"""

	# Current radius of shell 1 and shell_2['RADIUS']
	r10 = shell_1['RADIUS']
	r20 = shell_2['RADIUS']
	# Lorentz factors of shell 1 and shell 2
	g1 = shell_1['GAMMA']
	g2 = shell_2['GAMMA']

	# Calculate the collision time 
	t_coll = (2*g1**2)*(r20-r10)/(cc.c*((g1/g2)**2 - 1))

	# Check for nan value
	if t_coll != t_coll:
		t_coll = 1e99

	return t_coll

def shell_coll_gamma(s1g,s2g,s1m,s2m):
	"""
	Find the resulting Lorentz factor of two collided shells from the Lorentz factor and Mass of shell 1 and 2
	"""

	# The approximate Lorentz factor of the combined layer after complete redistribution of momentum and energy
	gamma_comb = np.sqrt(s1g*s2g*((s1m*s1g + s2m*s2g)/(s1m*s2g + s2m*s1g)))

	return gamma_comb

def vel(gamma):
	"""
	Method to calculate velocity (m/s) from Lorentz factor (gamma)
	"""

	return cc.c * np.sqrt(1-(1/gamma**2) )


def is_sorted(x):
	"""
	Method to test whether the array, x, is sorted
	Sorted in this context means x[i] > x[i+1]
	"""
	return (np.diff(x)<=0).all()