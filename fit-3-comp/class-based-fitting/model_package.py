import numpy as np
from astropy.modeling import Fittable1DModel, CompoundModel, Parameter

class ModelComponent():
	""" 
	Definition of ModelComponent class
	"""
	color = "C0"

	def print_params(self):
		"""
		Print current parameters
		"""

		for param in self.param_names:
			print("Name: {}\n\tdescription: {}\n\tvalue = {}\n".format(param,getattr(type(self), param)._description,getattr(type(self), param).value))

		return 0; 

class Band(Fittable1DModel,ModelComponent):
	"""
	One dimensional Band function.

	Parameters
	----------
	e0 : float
		Break energy
	alpha : float
		Low energy power law index
	beta : float
		High energy power law index
	norm : float
		Model normalization
	"""
	# name = "Band2"
	# _n_inputs = 1
	# _inputs = ('x',)
	# _n_outputs = 1
	# _outputs = ('y',)
	e0 = Parameter(default=1e3,name="Break Energy", description="Break energy",min=1e-99,max=1e10)
	alpha = Parameter(default=-1,name="alpha",description="Low energy power law index",min=-10,max=5)
	beta = Parameter(default=-2.5,name="beta",description="High energy power law index",min=-10,max=5)
	norm = Parameter(default=1,name="Normalization",description="Normalization",min=1e-99)
	
	@staticmethod
	def evaluate(energy, e0, alpha, beta, norm):
		"""
		Compute the broken power law spectrum at a particular energy
		"""

		# Initialize the return value
		flux_value = np.zeros_like(energy,subok=False)

		# Calculate peak energy
		e_peak = (alpha - beta)*e0

		i = energy <= e_peak
		if i.max():
			flux_value[i] = norm * np.power(energy[i]/100., alpha) * np.exp(- energy[i] / e0)
		
		i = energy > e_peak
		if i.max():
			flux_value[i] = norm * np.power((alpha - beta) * e0/100., alpha - beta) * np.exp(beta - alpha) * np.power(energy[i]/100,beta)

		return flux_value

	@staticmethod
	def fit_deriv(energy, e0, alpha, beta, norm):
		"""
		Compute derivative of the model
		"""

		# Initialize the return value
		deriv_flux_value = np.zeros_like(energy,subok=False)

		# Calculate peak energy
		e_peak = (alpha - beta)*e0

		i = energy <= e_peak
		if i.max():
			deriv_flux_value[i] = norm *(1/100.)* np.power(energy[i]/100.,alpha-1)*np.exp(-energy[i]/e0)*( (alpha) - (energy[i]/e0) )
		
		i = energy > e_peak
		if i.max():
			deriv_flux_value[i] = norm * np.power((alpha - beta) * e0/100., alpha - beta) * np.exp(beta - alpha) * (beta/100.) * np.power(energy[i]/100.,beta-1)

		return deriv_flux_value



class Blackbody(Fittable1DModel,ModelComponent):
	"""
	One dimensional Blackbody function.

	Parameters
	----------
	temp : float
		Blackbody temperature (in units of energy, i.e., k_B*T where k_B is the Boltzmann constant)
	alpha : float
		Index of the power law below temperature
	norm : float
		Model normalization
	"""
	# _name = "Blackbody"
	# _n_inputs = 1
	# _inputs = ('x',)
	# _n_outputs = 1
	# _outputs = ('y',)
	temp = Parameter(default=20,name="Temperature",description="Temperature",min=1e-99,max=1e10)
	alpha = Parameter(default=-0.4,name="alpha",description="Power law index",min=-10,max=5)
	norm = Parameter(default=1,name="Normalization",description="Normalization",min=1e-99)

	@staticmethod
	def evaluate(energy, temp, alpha, norm):
		"""
		Compute the black body spectrum at a particular energy
		"""

		# Initialize the return value
		flux_value = np.zeros_like(energy,subok=False)

		i = energy < 2e3
		if i.max():
			# If the energy is less than 2 MeV
			flux_value[i] = norm * np.power(energy[i]/temp,1.+alpha)/(np.exp(energy[i]/temp) - 1.)
		i = energy >= 2e3
		if i.max():
			flux_value[i] = 0

		return flux_value

	@staticmethod
	def fit_deriv(energy, temp, alpha, norm):
		"""
		Compute derivative of the model
		"""
		# Initialize the return value
		deriv_flux_value = np.zeros_like(energy,subok=False)

		i = energy < 2e3
		if i.max():
			# If the energy is less than 2 MeV
			deriv_flux_value[i] = np.power(energy[i]/temp,alpha) * ( (alpha+1)*temp*(np.exp(energy[i]/temp)-1) - energy[i]*np.exp(energy[i]/temp) ) / ( temp**2 *(np.exp(energy[i]/temp)-1)**2 )
		i = energy >= 2e3
		if i.max():
			deriv_flux_value[i] = 0

		return deriv_flux_value

def make_3comp(
	TH_Tp=50,TH_alpha=0.4,TH_norm=1,
	nTH1_e0=1e3,nTH1_alpha=-0.6,nTH1_beta=-2.5,nTH1_norm=1,
	nTH2_e0=1e8,nTH2_alpha=-1.1,nTH2_beta=-2.5,nTH2_norm=1):
	"""
	Function to quickly generate a three component model
	"""

	comp1 = Blackbody(temp=TH_Tp,alpha=TH_alpha,norm=TH_norm)
	comp1.color="r"
	comp1.alpha.fixed = True
	
	comp2 = Band(e0=nTH1_e0,alpha=nTH1_alpha,beta=nTH1_beta,norm=nTH1_norm)
	comp2.color="C0"
	comp2.alpha.fixed = True
	comp2.beta.fixed = True

	comp3 = Band(e0=nTH2_e0,alpha=nTH2_alpha,beta=nTH2_beta,norm=nTH2_norm)
	comp3.color="C2"
	comp3.alpha.fixed = True
	comp3.beta.fixed = True
	
	model_name = comp1 + comp2 + comp3

	return model_name



###
# Below attributes are added to the CompoundModel class of astropy
###
def print_info(self):
	"""
	Print current models and parameters of the compound model
	"""

	for i in range(self.n_submodels):
		print("Model Component: {}".format(self[i].__class__.__name__))
		for param in self[i].param_names:
			print("Parameter name: {}\n\tdescription: {}\n\tvalue = {}".format(param,getattr(type(self[i]), param)._description,getattr(type(self[i]), param).value))
		print("----------------")

	return 0; 

setattr(CompoundModel,'print_info', print_info) # Add information printing method
setattr(CompoundModel,'color', 'k') # Add color, used when plotting total spectrum
