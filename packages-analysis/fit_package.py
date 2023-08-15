import numpy as np
from astropy.modeling import fitting
from collections import defaultdict
from asymmetric_uncertainty import UncertaintyArray, neg_errors, pos_errors
import corner

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
		Perform a fit to the data using the supplied model
		"""

		# Fit data using specified model
		self.best_fit_model = self.fitter(model, data['ENERGY'],data['RATE'],maxiter=5000)

		# Calculate the new parameter uncertainties
		self.set_param_uncs()

		# Calculate the fit statistic for this fit 
		self.fit_stat = self.calc_fit_stat(self.best_fit_model,data)

		if verbose == True:
			if hasattr(self.best_fit_model, "model_name"):
				print(self.best_fit_model.model_name())
			else:
				print(self.best_fit_model.name)
			print("Best fit parameters:")
			for i in range(len(self.best_fit_model.parameters)):
				print("\t{} = {:.3f} +/- {:.3f}".format(self.best_fit_model.param_names[i],self.best_fit_model.parameters[i],getattr(self.best_fit_model,self.best_fit_model.param_names[i]).unc))
			print("Fit statistic = {:.3f}".format(self.fit_stat))


		return self.best_fit_model, self.fit_stat

	def fit_mc(self,model,data: Data, verbose = False,iterations=1000,num_sigma=5,discretization=100,ret_corner=False):
		"""
		Perform a Monte-Carlo exploration fit of the data using the supplied model 
		"""

		# Initializations
		dists = defaultdict(list) # list to hold the temporary best-fit parameters 
		likelihood_max = -1e99 # max likelihood
		best_fit_params_ind = None # best-fit parameters

		for i in range(iterations):
			# Create  data set that is the base data, but with a random offset scaled to match the uncertainties we used to create the data in the first place
			temp_data = np.zeros(shape=len(data))

			# Randomly fluctuate data
			for j in range(len(data)):
				locs_neg = np.linspace(data['RATE'][j].value-(num_sigma*data['RATE'][j].minus),data['RATE'][j].value,discretization)
				locs_pos = np.linspace(data['RATE'][j].value,data['RATE'][j].value+(num_sigma*data['RATE'][j].plus),discretization)
				locs = np.array(list(locs_neg)+list(locs_pos))

				# Generate random number to sample CDF 
				r = np.random.uniform(low=0, high=1)
				# Find the index of the location along the CDF that has the same value as the randomly generated number
				ind = np.argmax( (data['RATE'][j].cdf(locs) >= r) & (data['RATE'][j].cdf(locs) > r ) )
				# Store the randomly generated location 
				temp_data[j] = locs[ind]

			# Randomly select new initial condition for model parameters in the next fit 
			# I guess here I could implement a prior, but lets stick with uniform
			# for j in range(len(model.parameters)):
			# 	# Array of values within the bound
			# 	poss_vals = np.linspace(getattr(model,model.param_names[j]).min,
			# 	 						getattr(model,model.param_names[j]).max,
			# 							num=100)
			# 	# Random integer
			# 	rand_int = np.random.randint(low=0,high=100)
			# 	# Set random value 
			# 	setattr(model,model.param_names[j],poss_vals[rand_int])

			# Make a temporary fit that we just discard at the end of this loop after extracting the parameter values
			modeli = self.fitter(model, data['ENERGY'], temp_data,maxiter=5000)

			# Calculate likelihood of this fitted model instance.
			likelihood_tmp = self.calc_fit_stat(model=modeli,data=data,fit_stat_type=1)

			# Update max likelihood and best-fit parameters
			if likelihood_tmp > likelihood_max:
				likelihood_max = likelihood_tmp
				best_fit_params_ind = i

			for nm in modeli.param_names:
				dists[nm].append(getattr(modeli, nm).value)

		# Store the best-fit parameters in a list.
		dists = {k: np.array(v) for k,v in dists.items()}
		distarr = np.array([dists[nm] for nm in model.param_names]).T

		# Update best fit model 
		self.best_fit_model = model.deepcopy()
		self.best_fit_model.parameters = distarr[best_fit_params_ind]

		# Calculate the new parameter uncertainties
		self.set_param_uncs()

		# Make a corner plot if desired
		if ret_corner is True:
			corn_plot = corner.corner(distarr, labels=model.param_names, show_titles=True, title_fmt='.3f', truths=model.parameters)

		# Print out likelihood regions and best fit parameters if desired
		if verbose == True:
			if hasattr(self.best_fit_model, "model_name"):
				print(self.best_fit_model.model_name())
			else:
				print(self.best_fit_model.name)
			print("Best fit parameters:")
			for i in range(len(self.best_fit_model.parameters)):
				print("\t{} = {:.3f} +/- {:.3f}".format(self.best_fit_model.param_names[i],self.best_fit_model.parameters[i],getattr(self.best_fit_model,self.best_fit_model.param_names[i]).unc))
			print("Fit statistic = {:.3f}".format(self.fit_stat))

		if ret_corner is True:
			return self.best_fit_model, self.fit_stat, corn_plot
		else:
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


	def calc_fit_stat(self,model, data : Data,fit_stat_type=0):
		"""
		Calculate the fit statistic from the current model and data
		"""

		if fit_stat_type is None:
			fit_stat_type = self.fit_stat_type

		# Chi-squared
		if fit_stat_type == 0:
			# Determine which uncertainty to use by finding if the model above or below the respective data point
			sigma_arr = np.asarray(neg_errors(data['RATE']))
			sigma_arr[ np.flatnonzero(data['RATE'] > model(data['ENERGY'])) ] = np.asarray(pos_errors(data['RATE'][ data['RATE'] > model(data['ENERGY']) ]))

			# Calculate chi squared
			self.fit_stat = np.sum( np.power( model(data['ENERGY']) - UncertaintyArray(data['RATE']).values, 2) / np.power( sigma_arr, 2) )/(len(data)-len(model.parameters))

		# Likelihood
		if fit_stat_type == 1:
			# Determine which uncertainty to use by finding if the model above or below the respective data point
			sigma_arr = np.asarray(neg_errors(data['RATE']))
			sigma_arr[ np.flatnonzero(data['RATE'] > model(data['ENERGY'])) ] = np.asarray(pos_errors(data['RATE'][ data['RATE'] > model(data['ENERGY']) ]))

			# Calculate chi squared terms
			chi2_terms =  np.power( (model(data['ENERGY']) - UncertaintyArray(data['RATE']).values) /sigma_arr,2)
			# Calculate likelihood
			self.fit_stat = np.sum((-chi2_terms/2)-np.log(np.sqrt(2*np.pi*np.power(sigma_arr,2))))

		return self.fit_stat
