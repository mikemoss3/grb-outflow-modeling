import numpy as np
import scipy.optimize as opt

import Model
import Data

class FitPackage(object):
	def __init__(self):
		# Initialize Fit type
		# Initialize Fit statistic value
		self.fit_stat = 0.


	def fit(self,model : Model, data : Data):
		"""
		Perform a fit between the input model and data
		"""

		opt.curve_fit(model, data['ENERGY'], data['RATE'])