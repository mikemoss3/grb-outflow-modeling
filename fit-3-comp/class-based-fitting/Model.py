import numpy as np
import ModelComponent

# File used to define Model class

class Model(object):
	def __init__(self):
		self.num_comps = 0
		self.model_comps = np.zeros(shape=0,dtype=ModelComponent)
		self.model_params = None
		self.spectral_function = None

	def add_component(self, comp : ModelComponent):
		self.num_comps = +1 
		self.model_params.append(comp.params)

		self._make_spectral_function()

	def _make_spectral_function(self):
		# Make composite spectral function.
		# Add the spectral function of each component to the total spectral function

		self.spectral_function = lambda energy: [self.model_comps[i].spectral_function(energy) for i in range(self.num_comps)]






