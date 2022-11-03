import numpy as np
import ModelComponent


# Definition of Band component 

class BAND(ModelComponent):
	def __init__(self,e0=1e4,alpha=-1,beta=-2.5,norm=1):
		ModelComponent.__init__(self)

		# Model Parameters: 
		# 	Peak energy
		#	alpha
		#	beta
		#	normalization
		self.param_num = 4
		self.param_vals = np.array( (e0,alpha,beta,norm), dtype = [("e0",float),("alpha",float),("beta",float),("norm",float)] )

	def spectral_function(self,energy):
		"""
	    Compute the broken power law spectrum at a particular energy
		"""
		if energy < (self.param_vals["alpha"] - self.param_vals['beta'])*self.param_vals["e0"]:
			return self.param_vals["norm"] * np.power(energy/100., self.param_vals["alpha"]) * np.exp(- energy / self.param_vals["e0"])
		else:
			return self.param_vals["norm"] * np.power((self.param_vals['alpha'] - self.param_vals['beta']) * self.param_vals['e0']/100., self.param_vals["alpha"] - self.param_vals['beta']) * np.exp(self.param_vals['beta'] - self.param_vals['alpha']) * np.pow(energy /100,self.param_vals['beta'])
