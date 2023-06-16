import numpy as np
from astropy.modeling import fitting

from model_package import *
from data_package import *

class FittedModel(object):
	def __init__(self,fitter = fitting.LMLSQFitter(), fit_stat = 0):
		# Initialize fitter type
		self.fitter = fitter
		# Initialize fit statistic type and value
		self.fit_stat_type = 0
		self.fit_stat = fit_stat

		# Initialize best fit model and param uncertanties
		self.best_fit_model = None
		self.param_uncs = None

	def set_fitter(self,fitter):
		"""
		Method to set the fit type to use. 
		"""
		self.fitter = fitter

	def fit(self,model, data : Data,verbose=False):
		"""
		Perform a fit between the input model and data
		"""

		# Fit data using specified model
		self.best_fit_model = self.fitter(model, data['ENERGY'],data['RATE'],maxiter=5000,weights=1.0/data['ERR'])

		# Calculate the new parameter uncertainties
		self.set_param_uncs()

		# Calculate the fit statistic for this fit 
		self.fit_stat = self.calc_fit_stat(self.best_fit_model,data)

		if verbose == True:
			print(self.best_fit_model.model_name())
			print("Best fit parameters:")
			for i in range(len(self.best_fit_model.parameters)):
				print("\t{} = {:.3f} +/- {:.3f}".format(self.best_fit_model.param_names[i],self.best_fit_model.parameters[i],getattr(self.best_fit_model,self.best_fit_model.param_names[i]).unc))
			print("Fit statistic = {:.3f}".format(self.fit_stat))


		return self.best_fit_model, self.fit_stat

	def set_param_uncs(self):
		"""
		Calculate the uncertainty on the parameters of the fitted model
		"""
		self.param_uncs = np.zeros(shape=len(self.best_fit_model.parameters))
		
		cov = np.diag(self.fitter.fit_info['param_cov']) # Grab covariance matrix
		j=0

		for i in range(len(self.best_fit_model.parameters)):
			if getattr(self.best_fit_model,self.best_fit_model.param_names[i]).fixed == False:
				getattr(self.best_fit_model,self.best_fit_model.param_names[i]).unc = cov[j]**0.5
				j+=1


	def calc_fit_stat(self,model, data : Data):
		"""
		Calculate the fit statistic from the current model and data
		"""

		if self.fit_stat_type == 0:
			self.fit_stat = np.sum( np.power( model(data['ENERGY']) - data['RATE'], 2) / np.power( data['ERR'], 2) )/(len(data)-len(model.parameters))

		return self.fit_stat
