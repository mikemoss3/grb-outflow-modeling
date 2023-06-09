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

	def model_name(self):
		"""
		Name of the compound model, automatically set as the combination of submodel names
		"""
		return self.name;

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
	name = "Band"
	# _n_inputs = 1
	# _inputs = ('x',)
	# _n_outputs = 1
	# _outputs = ('y',)
	e0 = Parameter(default=1e3,name="Break Energy", description="Break energy",min=1e-99,max=1e10)
	alpha = Parameter(default=-1,name="alpha",description="Low energy power law index",min=-1.99,max=2)
	beta = Parameter(default=-2.5,name="beta",description="High energy power law index",min=-10,max=-2)
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
		Compute the partial derivatives for each input parameter of the model
		"""
		# Initialize the return values
		d_e0 = np.zeros_like(energy,subok=False) # partial derivative of the break energy
		d_alpha = np.zeros_like(energy,subok=False)# partial derivative of the low energy power law index
		d_beta = np.zeros_like(energy,subok=False) # partial derivative of the high energy power law index
		d_norm = np.zeros_like(energy,subok=False) # partial derivative of the normalization

		# Calculate peak energy
		e_peak = (alpha - beta)*e0

		# Calculate partial derivative of the break energy
		i = energy <= e_peak
		if i.max():
			d_e0[i] = norm * np.power(100,-alpha) * np.power(energy[i],alpha+1) * np.exp(-(energy[i]/e0) - 2) * (np.log(e0) - 1)
		i = energy > e_peak
		if i.max():
			d_e0[i] = norm * np.power(100,-alpha) * np.exp(beta - alpha) * np.power(energy[i],beta) * np.power((alpha-beta)*e0,alpha-beta+1) / e0**2

		# Calculate partial derivative of the low energy power law index
		i = energy <= e_peak
		if i.max():
			d_alpha[i] = norm * np.power(100,-alpha) * np.power(energy[i],alpha) * np.exp(-energy[i]/e0) * np.log(energy[i]/100)
		i = energy > e_peak
		if i.max():
			d_alpha[i] =  norm * np.power(100,-alpha) * np.exp(beta-alpha) * np.power(energy[i],beta) * np.power((alpha-beta)*e0,alpha-beta) * np.log((alpha-beta)*e0/100)

		# Calculate partial derivative of the high energy power law index
		i = energy <= e_peak
		if i.max():
			d_beta[i] = 0
		i = energy > e_peak
		if i.max():
			d_beta[i] = norm * np.power(100,-alpha) * np.exp(beta-alpha) * np.power(energy[i],beta) * np.power((alpha-beta)*e0,alpha-beta) * (np.log(e0*(alpha - beta)) - np.log(energy[i]))

		# Calculate partial derivative of the break norm
		i = energy <= e_peak
		if i.max():
			d_norm[i] = np.power(energy[i]/100., alpha) * np.exp(- energy[i] / e0)
		i = energy > e_peak
		if i.max():
			d_norm[i] = np.power((alpha - beta) * e0/100., alpha - beta) * np.exp(beta - alpha) * np.power(energy[i]/100,beta)

		return [d_e0, d_alpha, d_beta, d_norm]


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
	name = "Blackbody"
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
		Compute the partial derivatives for each input parameter of the model
		"""
		# Initialize the return values
		d_temp = np.zeros_like(energy,subok=False) # partial derivative of the temperature
		d_alpha = np.zeros_like(energy,subok=False) # partial derivative of the power law index
		d_norm = np.zeros_like(energy,subok=False) # partial derivative of the normalization

		# Calculate partial derivative of the temperature
		i = energy < 2e3
		if i.max():
			d_temp[i] = norm * np.power(energy[i]/temp,alpha)*(energy[i] *np.exp(energy[i]/temp) - alpha*temp*(np.exp(energy[i]/temp) - 1)) / ( temp**2 * (np.exp(energy[i]/temp) -1)**2  )
		i = energy >= 2e3
		if i.max():
			d_temp[i] = 0

		# Calculate partial derivative of the power law index
		i = energy < 2e3
		if i.max():
			d_alpha[i] = norm * np.power(energy[i]/temp,alpha) * np.log(energy[i]/temp) / (np.exp(energy[i]/temp) -1)
		i = energy >= 2e3
		if i.max():
			d_alpha[i] = 0

		# Calculate partial derivative of the normalization
		i = energy < 2e3
		if i.max():
			d_norm[i] = np.power(energy[i]/temp,1.+alpha)/(np.exp(energy[i]/temp) - 1.)
		i = energy >= 2e3
		if i.max():
			d_norm[i] = 0

		return [d_temp, d_alpha, d_norm]

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

def model_name(self):
	"""
	Name of the compound model, automatically set as the combination of submodel names
	"""

	name = self[0].name
	for i in range(1,self.n_submodels):
		name += ", {}".format(self[i].name)

	return name;


setattr(CompoundModel,'print_info', print_info) # Add information printing method
setattr(CompoundModel,'color', 'k') # Add color, used when plotting total spectrum
setattr(CompoundModel,'model_name', model_name) # Add name for compound model


# Add uncertainty to Parameter class
setattr(Parameter, 'unc', np.nan)