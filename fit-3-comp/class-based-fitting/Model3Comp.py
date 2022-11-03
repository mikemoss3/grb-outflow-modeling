import numpy as np
import Model

import ModelComponent
import BAND
import BB

# File to define the three component model

class Model3Comp(Model):
	def __init__(self):
		Model.__init__(self)

		self.add_component(BB())
		self.add_component(BAND(alpha=-0.6))
		self.add_component(BAND(e0=1e8,alpha=-1.1))