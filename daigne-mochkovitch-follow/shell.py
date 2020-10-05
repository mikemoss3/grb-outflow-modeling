"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2020-08-17


This code defines a Shell object class. This class contains all relevant information and definitions of each shell/layer in a 
GRB prompt jet simulation. Many of these Shell objects are used to create the total jet.

"""

class Shell(object):
	"""
	Shell class.
	"""

	def __init__(self,radius=10e12,gamma=1,mass=1):
		"""
		Defines the default parameters of a shell.
		"""
		# Radius in the jet (distance from the central engine), in cm
		self.radius = radius
		# Lorentz factor
		self.gamma = gamma
		# Mass, in kg
		self.mass = mass

		# 
		self.temperature = 10.