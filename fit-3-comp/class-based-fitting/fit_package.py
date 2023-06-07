import numpy as np
from astropy.modeling import fitting

from model_package import *
from data_package import *

class FitPackage(object):
	def __init__(self):
		# Initialize Fit type
		self.fit_type = 0 # chi-squared
		# Initialize Fit statistic value
		self.fit_stat = 0.

	def set_fit_type(self,fit_type):
		"""
		Method to set the fit type to use. 
		"""
		self.fit_type = fit_type # chi-squared 

	def fit(self,model, data : Data):
		"""
		Perform a fit between the input model and data
		"""

		# Fit data using specified model
		fitter = fitting.LevMarLSQFitter()
		best_fit_model = fitter(model, data['ENERGY'],data['RATE'],weights=1.0/data['ERR'])

		# Calculate the fit statistic for this fit 
		fit_stat = self.calc_fit_stat(best_fit_model,data)

		return best_fit_model, fit_stat

	def calc_fit_stat(self,model, data : Data):
		"""
		Calculate the fit statistic from the current model and data
		"""

		if self.fit_type == 0:
			self.fit_stat = np.sum( np.power( best_fit_model(data['ENERGY']) - data['RATE'], 2) / np.power( data['ERROR'], 2) ) 

		return self.fit_stat
