import numpy as np
import ModelComponent


# Definition of BlackBody component

class BB(ModelComponent):
	def __init__(self,Tp=50,alpha=0.4,norm=1):
		ModelComponent.__init__(self)

		# Model Parameters: 
		# 	peak temperature (keV)
		#	alpha
		#	normalization
		self.param_num = 3
		self.param_vals = np.array( (Tp,alpha,norm), dtype = [("Tp",float),("alpha",float),("norm",float)] )

	def spectral_function(self,energy):
		"""
	    Compute the black body spectrum at a particular energy
		"""
		return self.param_vals["norm"] * np.power(energy/self.param_vals["Tp"],1.+self.param_vals["alpha"])/(np.exp(energy/self.param_vals["Tp"]) - 1.)


