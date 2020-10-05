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

from shell import Shell
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
		self.tw = self.dt*self.numshells


		# Make the list of shells with a Lorentz distributions (this is the distribution at t = tw)
		self.shell_lorentz_arr = shelldist(numshells=self.numshells)

		plt.figure()
		ld.plot_lorentz_dist(self.shell_lorentz_arr,title='initial')

		# Initialize a spectrum object 
		self.spectrum = sp.Spectrum()

	def jet_evolution(self):
		"""
		Method to emulate the evolution of a jet made of consecutive shells (including the collisions of shells)
		"""

		# Define time step, in seconds
		# dts = 0.001
		dts = self.dt/2
		# dts = self.dt

		# Start time of the simulation, in seconds
		t = self.dt # Start at dt, this starts us with one initially launched shell. 

		# The new shell number keeps track of which shell to grab from the initial Lorentz Distribution 
		new_shell_num = 1

		# Array to store the currently launched shells in the jet 
		jet_shells = np.ndarray(shape=1,dtype=object)

		# Set up the first shell. 
		jet_shells[0] = self.shell_lorentz_arr[-new_shell_num]
		new_shell_num+=1
		jet_shells[0].radius += dts*vel(jet_shells[0].gamma)

		# For each time step after the start:
		# Move all launched shells forward (the distance they move depends on their respective Lorentz factors)
		# Check for a collision between shells, if there is a collision, make a new combined shell
		# Check if a new layer is produced (occurs every 2 ms.), this will only occur if the current number of shells is less than the total
		# Check if the shells are in order of decreasing Lorentz factor, then no more collisions are possible

		# The Lorentz factors are not in order at the beginning of the simulation
		lorentz_ord = False

		while lorentz_ord is False:
		# while t<90:
			# print('### Time = {} ###'.format(t))
			# for i in range(len(jet_shells)):
			# 	print(jet_shells[i].radius)
			# 	print(jet_shells[i].gamma)

			# Increase the time by one time step 
			t +=dts

			# Move shells forward
			for i in range(len(jet_shells)):
				jet_shells[i].radius += dts*vel(jet_shells[i].gamma)


			# Check for and handle collisions
			for i in range(len(jet_shells)-1):
				# If a shell has a radius greater than the next shell, they have collided
				if jet_shells[i].radius > jet_shells[i+1].radius:
					# print("Collision!")
					# Because there was a collision, the shock will produce radiation.
					gammaint = 0.5 * (np.sqrt(jet_shells[i].gamma/jet_shells[i+1].gamma)+np.sqrt(jet_shells[i+1].gamma/jet_shells[i].gamma)) # Lorentz factor for internal motion in shocked material
					if gammaint != 1:
						# Add the radiation to the spectrum
						ae = 0.3 # fraction of dissipated energy that goes into electrons
						ab = 0.3 # fraction of dissipated energy that goes into magnetic field
						Edot = 1e52/(4*np.pi) # Energy injection rate, erg/s/sr

						te = t # emission time, source frame
						ta = te - (jet_shells[i].radius/cc.c) # arrival time
						Qic = (-1 + np.sqrt(1+ (4*0.3*ae/ab)))/2 # Compton parameter 
						
						mu = 1.75 # index of the electron spectrum caused by magnetic fluctuations, range = 1.5 - 2.0
						am = 0.5 # fraction of dissipated energy which goes into magnetic fluctuations, range = 0.1 - 1.0
						ksi = 1e-3 # fraction of electrons which are accelerated
						eps = (gammaint - 1) * cc.mp * cc.c**2 # average dissipated energy per photon in a shock
						gammae = np.power(am*eps/(ksi*cc.me* cc.c**2),1/(3-mu)) # average Lorentz factor of the electron population accelerated from turbulent magnetic fields

						shared_gamma = shell_coll_gamma(jet_shells[i],jet_shells[i+1]) # combined Lorentz factor of the two colliding shells
						shared_radius = (jet_shells[i].radius + jet_shells[i+1].radius)/2 # radius of the collision
						n = Edot/(4*np.pi* shared_radius**2 * shared_gamma**2 * cc.mp * cc.c**3) # comoving number density 
						Beq = np.sqrt(8*np.pi*ab*n*eps) # Magnetic field 

						w =  33 * (Beq/1000) * (gammae/1e4)**3 # Compton scattering limit 
						# taus = Qic / gammae**2 # optical depth of shell 
						aic = Qic/w / (1+(Qic/w)) # fraction of energy dissipated due to Inverse Compton 
						asyn = 1 - aic # fraction of energy dissipated due to synchrotron radiation

						Esyn = 50*(shared_gamma/300)*(Beq/1000)*(gammae/100)**2 # synchrotron energy, in eV

						self.spectrum.add_contribution(te, ta, asyn, Beq, gammae, Esyn)

					# As for the mechanics of the two shells 
					# Set the Lorentz factor of the two shells to be the same resulting value
					jet_shells[i].gamma = shared_gamma
					jet_shells[i+1].gamma = shared_gamma

					# Set the two shells to the same (averaged) radius
					jet_shells[i].radius = shared_radius
					jet_shells[i+1].radius = shared_radius


			# Check if there are still shells to grab from the Lorentz distribution and it is time to launch a new jet
			# Checking if it is the correct time to launch a jet is affected by floating point errors, so the if statement is messy.
			if 	new_shell_num<len(self.shell_lorentz_arr)+1 and (np.mod(t,self.dt) == 0 or np.mod(t,self.dt) < 1e-10 or self.dt - np.mod(t,self.dt) < 1e-5):
				# Add new shell to the list
				jet_shells = np.insert(jet_shells,0,self.shell_lorentz_arr[-new_shell_num])
				# Increase the new shell count
				new_shell_num+=1 

			# Check Lorentz factor order (this must only be checked after all shells are launched!)
			if len(jet_shells) == self.numshells:
				if(all(np.abs((jet_shells[i].gamma/jet_shells[i + 1].gamma)-1)<1e-6 for i in range(len(jet_shells)-1))): 
				    print('Jet is ordered.')
				    lorentz_ord = True

			if t > 1e5:
				print('Reached t>1e5.')
				lorentz_ord=True

		# Three conditions must be met in order to contribute to the jet emission 
		# 1. The energy dissipation must occur on a timescale smaller than the shell-expansion time
		# 2. The relative velocity of the two layers must collide at a velocity larger than the sound speed in the medium
		# 3. The wind must be transparent to the emitted photons

			if (t > 3.4e4 and t < 3.4e4 + dts):
				plt.figure()
				ld.plot_lorentz_dist(jet_shells,title='Distribution at 3.4e4')

			if (t >= 1.7e5 and t <= 1.7e5 + dts):
				plt.figure()
				ld.plot_lorentz_dist(jet_shells,title='Distribution at 1.7e5')

		plt.figure()
		ld.plot_lorentz_dist(jet_shells,title='Final Distribution')


	def make_new_shell(self,shell_1: Shell, shell_2: Shell):
		"""
		This method assumes that shell_1 has a Lorentz factor > than the Lorentz factor of shell_2.
		(which also means the masses of the shells m_1 < m_2)
		"""

		# Make new shell
		shell_comb = Shell()

		# The approximate Lorentz factor of the combined layer after complete redistribution of momentum and energy
		shell_comb.gamma = np.sqrt(shell_1.gamma*shell_2.gamma*((shell_1.mass*shell_1.gamma + shell_2.mass*shell_2.gamma)/(shell_1.mass*shell_2.gamma + shell_2.mass*shell_1.gamma)))
		# The location of the new shell is where shell_1 and shell_2 collided
		# This is a bit of an approximation
		shell_comb.radius = shell_2.radius
		# Total combined shell mass 
		shell_comb.mass = shell_1.mass + shell_2.mass

		return shell_comb

def shell_coll_gamma(shell_1: Shell, shell_2: Shell):
	"""
	Find the resulting Lorentz factor of two collided shells
	"""

	# The approximate Lorentz factor of the combined layer after complete redistribution of momentum and energy
	gamma_comb = np.sqrt(shell_1.gamma*shell_2.gamma*((shell_1.mass*shell_1.gamma + shell_2.mass*shell_2.gamma)/(shell_1.mass*shell_2.gamma + shell_2.mass*shell_1.gamma)))

	return gamma_comb

def vel(gamma):
	"""
	Method to calculate velocity (m/s) from Lorentz factor (gamma)
	"""

	return 3e8 * np.sqrt(1-(1/gamma**2) )