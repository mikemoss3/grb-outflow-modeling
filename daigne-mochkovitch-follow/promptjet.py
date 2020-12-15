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

import lorentzdists as ld
import spectrum as sp
import cosmologicalconstants as cc


class PromptJet(object):
	"""
	Prompt Jet class.
	"""

	def __init__(self,numshells=5000,dt=0.002,shelldist=ld.step):
		"""
		Defines the default parameters of the jet.
		"""
		# Number of shells
		self.numshells = numshells
		# Time in between successive shell emissions from the central engine.
		self.dt = dt
		# Total duration of central engine activity (total time to produce numshells), in seconds 
		self.tw = self.dt*(self.numshells+1) # Plus a little lee-way 


		# Make the list of shells with a Lorentz distributions (this is the distribution at t = tw)
		self.shell_lorentz_arr = shelldist(numshells=self.numshells)

		# Initialize a spectrum object 
		self.spectrum = sp.Spectrum()

	def jet_evolution(self):
		"""
		Method to emulate the evolution of a jet made of consecutive shells (including the collisions of shells)
		"""

		# Define time step, in seconds
		dts = self.dt/2. # use some dts < dt 

		# Initialize start time to 0 sec
		t = 0

		# Initialize the jet with a first shell
		jet_shells = np.ndarray(shape=1,dtype=[('RADIUS',float),('GAMMA',float),('MASS',float)])
		i = 0 # Keeps track of the shell we a grabbing from the shell distribution 
		jet_shells['RADIUS'] = self.shell_lorentz_arr[i]['RADIUS'] 
		jet_shells['GAMMA'] = self.shell_lorentz_arr[i]['GAMMA'] 
		jet_shells['MASS'] = self.shell_lorentz_arr[i]['MASS'] 
		i+=1

		# For each time step:
		# Move all launched shells forward (the distance they move depends on their respective Lorentz factors)
		# Check for a collision between shells, if there is a collision, make a new combined shell
		# Check if a new layer is produced (occurs every 2 ms.), this will only occur if the current number of shells is less than the total
		# Check if the shells are in ordered Lorentz factor, then no more collisions are possible

		# Flag to mark if all the shells in the jet are ordered (no more collisions occur at this point )
		ord_lorentz = False
		while ord_lorentz is False:
			# Increase time step
			t += dts

			# To add another layer, first check if the number of shells currently launched is lower than the total 
			# number of shells 
			if (len(jet_shells) < self.numshells):
				# Then check if its time to launch another shell
				if (np.mod(t,dts)-self.dt)<1e-7:
					r = self.shell_lorentz_arr[i]['RADIUS'] 
					g = self.shell_lorentz_arr[i]['GAMMA'] 
					m = self.shell_lorentz_arr[i]['MASS'] 
					i+=1
					jet_shells = np.append(jet_shells,np.array((r,g,m),dtype=[('RADIUS',float),('GAMMA',float),('MASS',float)]))
			else:
				# If all the shells have been launched, its time to check if the Lorentz factors shells are ordered 
				# There will be no more collisions once this occurs.
				if all((1 - (jet_shells['GAMMA'][k]/jet_shells['GAMMA'][k+1])) < 1e-4 for k in range(len(jet_shells['GAMMA'])-1)):
					ord_lorentz = True
					print('Ordered Lorentz factors in the jet.')


			# Move shells forward
			jet_shells['RADIUS'] = jet_shells['RADIUS']+ (vel(jet_shells['GAMMA'])*dts)

			# Check for collisions, for each shell
			for j in range(len(jet_shells)-1):
				# Check if the shells current position is now <= to the shell launched after it
				if jet_shells['RADIUS'][j] < jet_shells['RADIUS'][j+1]:
					# Calculate the combined Lorentz factor of the two shells
					gamma_comb = shell_coll_gamma(s1g=jet_shells[j]['GAMMA'],s2g=jet_shells[j+1]['GAMMA'],s1m=jet_shells[j]['MASS'],s2m=jet_shells[j+1]['MASS'])
					

					### Calculate spectral addition ### 
					gamma_int = 0.5*(np.sqrt(jet_shells['GAMMA'][j]/jet_shells['GAMMA'][j+1])+np.sqrt(jet_shells['GAMMA'][j+1]/jet_shells['GAMMA'][j])) # Lorentz factor for internal motion in shocked material
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

					n = E_dot / (4*np.pi * jet_shells[j]['RADIUS']**2 * gamma_comb**2 * cc.mp * cc.c**3) # Comoving proton number density
					alpha_b = 1/3 # Fraction of dissipated energy that goes into the magnetic field 
					B = np.sqrt(8*np.pi*alpha_b*n*eps) # Magnetic field strength, assuming B = Beq  
					gamma_r = np.sqrt(jet_shells['GAMMA'][j]*jet_shells['GAMMA'][j+1]) # Approximate resulting Lorentz factor from the shell collision
					E_syn = 50*(gamma_r/300)*(B/1000)*(gamma_e/100)**2 # eV, Synchrotron energy in the rest frame

					E_0syn = E_syn/gamma_r # Synchrotron energy in the comoving frame
					w = gamma_e*E_0syn/(cc.me*cc.c**2) # Critical value between Thomson and Klein-Nishina
					t_var = jet_shells[j]['RADIUS']/(cc.c*gamma_r**2) # Dynamical time scale of the shell

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
						ta = t - (gamma_comb/cc.c)
						self.spectrum.add_contribution(te=t,ta=ta,asyn=alpha_syn,Beq=B,gammae=gamma_e,Esyn=E_syn)


					### Calculate shell dynamics ### 
					# Calculate the radius of the two shells
					jet_shells['RADIUS'][j+1] = jet_shells['RADIUS'][j+1] - (vel(jet_shells['GAMMA'][j+1])*dts) + (vel(gamma_comb)*dts)
					jet_shells['RADIUS'][j] = jet_shells['RADIUS'][j+1]
					# Set the Lorentz factor
					jet_shells[j]['GAMMA']=gamma_comb
					jet_shells[j+1]['GAMMA']=gamma_comb


					


			
			# Testing Purposes:
			# if t > 1e4:
				# ord_lorentz = True
				# plt.figure()
				# ld.plot_lorentz_dist(jet_shells,title='Time = {} sec'.format(t))
			

		# plt.figure()
		# ld.plot_lorentz_dist(jet_shells,title='Time = {} sec'.format(t))


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

	return 3e10 * np.sqrt(1-(1/gamma**2) )