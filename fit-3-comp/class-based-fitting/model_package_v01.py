import numpy as np
from astropy.modeling import Fittable1DModel, Parameter

class ModelComponent():
	""" 
	Definition of ModelComponent class
	"""
	def __init__(self):
		self.color = "C0"
		self.param_list = []

	def set_param(self,param_num, new_value):
		"""
		Set a parameter's value
		"""

		self.param_list[param_num].value = new_value

		return 0

	def print_params(self):
		"""
		Print current parameters
		"""

		for param in self.param_list:
			print("Name: {}\n\tdescription: {}\n\tvalue = {}\n".format(param.name,param._description,param.value))

		return 0; 

	def spectral_function(self,energy):
		"""
		Un-defined spectral function 
		"""
		print("No emission component has been set.")
		return;

	def make_spectra(self,emin,emax,num=None,spec_type=0):
		"""
		Make a spectra for the model within the defined energy range

		spectra type:
			0 == count spectra
			1 == flux spectra 
			2 == power spectra
		"""
		if num == None:
			num = int(np.log10(emax/emin)*20)
		spectra = np.zeros(shape=num,dtype=[("ENERGY",float),("RATE",float)])
		spectra['ENERGY'] = np.logspace(np.log10(emin),np.log10(emax),num=num)
		if spec_type == 0:
			spectra['RATE'] = self.spectral_function(spectra['ENERGY'])	
		if spec_type == 1:
			spectra['RATE'] = self.spectral_function(spectra['ENERGY'])	* spectra['ENERGY']
		elif spec_type == 2:
			spectra['RATE'] = self.spectral_function(spectra['ENERGY'])	* spectra['ENERGY']**2

		return spectra

class BAND(Fittable1DModel,ModelComponent):
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
	_name = "Band"
	_n_inputs = 1
	_inputs = ('x',)
	_n_outputs = 1
	_outputs = ('y',)
	e0 = Parameter(default=1,name="Break Energy", description="Break energy", min = 0)
	alpha = Parameter(default=1,name="alpha",description="Low energy power law index")
	beta = Parameter(default=1,name="beta",description="High energy power law index")
	norm = Parameter(default=1,name="Normalization",description="Normalization", min = 0)
	
	def __init__(self,e0=1e4,alpha=-1,beta=-2.5,norm=1):
		ModelComponent.__init__(self)

	# 	self._name = "Band"
	# 	self._n_inputs = 1
	# 	self._inputs = ('x',)
	# 	self._n_outputs = 1
	# 	self._outputs = ('y',)

	# 	self.param_num = 4
	# 	# Model Parameters: 
	# 	self.e0 = Parameter(default=e0,name="Break Energy", description="Break energy", min = 0)
	# 	self.alpha = Parameter(default=alpha,name="alpha",description="Low energy power law index")
	# 	self.beta = Parameter(default=beta,name="beta",description="High energy power law index")
	# 	self.norm = Parameter(default=norm,name="Normalization",description="Normalization", min = 0)

	# 	self.param_list = [self.e0,self.alpha,self.beta,self.norm]

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
	def fit_deriv(energy, temp, alpha, norm):
		"""
		Compute derivative of the model
		"""
		a = "will do"
		return a

	# def spectral_function(self,energy):
	# 	"""
	#     Compute the broken power law spectrum at a particular energy
	# 	"""

	# 	# Initialize the return value
	# 	flux_value = np.zeros_like(energy,subok=False)

	# 	# Calculate peak energy
	# 	e_peak = (self.alpha.value - self.beta.value)*self.e0.value

	# 	i = energy <= e_peak
	# 	if i.max():
	# 		flux_value[i] = self.norm.value * np.power(energy[i]/100., self.alpha.value) * np.exp(- energy[i] / self.e0.value)
		
	# 	i = energy > e_peak
	# 	if i.max():
	# 		flux_value[i] = self.norm.value * np.power((self.alpha.value - self.beta.value) * self.e0.value/100., self.alpha.value - self.beta.value) * np.exp(self.beta.value - self.alpha.value) * np.power(energy[i]/100,self.beta.value)

	# 	return flux_value

class BB(Fittable1DModel,ModelComponent):
	"""
	One dimensional Band function.

	Parameters
	----------
	temp : float
		Blackbody temperature (in units of energy, i.e., k_B*T where k_B is the Boltzmann constant)
	alpha : float
		Index of the power law below temperature
	norm : float
		Model normalization
	"""
	def __init__(self,temp=50,alpha=0.4,norm=1):
		ModelComponent.__init__(self)

		self._name = "Blackbody"

		# Model Parameters: 
		# 	peak temperature (keV)
		#	alpha
		#	normalization
		self.param_num = 3
		self.temp = Parameter(default=temp,name="Temperature",description="Temperature",min=0)
		self.alpha = Parameter(default=alpha,name="alpha",description="Power law index")
		self.norm = Parameter(default=norm,name="Normalization",description="Normalization",min=0)

		self.param_list = [self.temp,self.alpha,self.norm]

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
		a = "will do"
		return a

	# def spectral_function(self,energy):
	# 	"""
	# 	Compute the black body spectrum at a particular energy
	# 	"""

	# 	# Initialize the return value
	# 	flux_value = np.zeros_like(energy,subok=False)

	# 	i = energy < 2e3
	# 	if i.max():
	# 		# If the energy is less than 2 MeV
	# 		flux_value[i] = self.norm.value * np.power(energy[i]/self.temp.value,1.+self.alpha.value)/(np.exp(energy[i]/self.temp.value) - 1.)
	# 	i = energy >= 2e3
	# 	if i.max():
	# 		flux_value[i] = 0

	# 	return flux_value


class Model(object):
	def __init__(self,component:ModelComponent=None):
		self.num_comps = 0
		self.model_comps = np.zeros(shape=0,dtype=ModelComponent)
		self.spectral_function = None

		if component is not None:
			self.num_comps += 1
			self.model_comps = np.append(self.model_comps,component)
			self._make_spectral_function()

	def add_component(self, component:ModelComponent):
		"""
		Method to add spectral component to model
		"""
		self.num_comps += 1
		self.model_comps = np.append(self.model_comps,component)

		self._make_spectral_function()

	def _make_spectral_function(self):
		"""
		Make composite spectral function.
		Add the spectral function of each component to the total spectral function
		"""

		self.spectral_function = lambda energy: ( np.sum(self.model_comps[i].spectral_function(energy) for i in range(self.num_comps)) )

	def make_spectra(self,emin,emax,num=None,spec_type=0):
		"""
		Make a spectra for the model within the defined energy range

		spectra type:
			0 == count spectra
			1 == power spectra 
			2 == nu F nu spectra
		"""
		if num == None:
			num = int(np.log10(emax/emin)*20)
		spectra = np.zeros(shape=num,dtype=[("ENERGY",float),("RATE",float)])
		spectra['ENERGY'] = np.logspace(np.log10(emin),np.log10(emax),num=num)
		if spec_type == 0:
			spectra['RATE'] = self.spectral_function(spectra['ENERGY'])	
		if spec_type == 1:
			spectra['RATE'] = self.spectral_function(spectra['ENERGY'])	* spectra['ENERGY']
		elif spec_type == 2:
			spectra['RATE'] = self.spectral_function(spectra['ENERGY'])	* spectra['ENERGY']**2

		return spectra

class Model3Comp(Model):
	def __init__(self,TH_Tp=50,TH_alpha=0.4,TH_norm=1,nTH1_e0=1e4,nTH1_alpha=-0.6,nTH1_beta=-2.5,nTH1_norm=1,nTH2_e0=1e8,nTH2_alpha=-1.1,nTH2_beta=-2.5,nTH2_norm=1):
		## TH_ == thermal component parameters
		## nTH1_ == first non-thermal component parameters
		## nTH2_ == second non-thermal component parameters
		Model.__init__(self)

		self.add_component(BB(temp=TH_Tp,alpha=TH_alpha,norm=TH_norm))
		self.model_comps[0].color="r"
		self.add_component(BAND(e0=nTH1_e0,alpha=nTH1_alpha,beta=nTH1_beta,norm=nTH1_norm)) 
		self.model_comps[1].color="C0"
		self.add_component(BAND(e0=nTH2_e0,alpha=nTH2_alpha,beta=nTH2_beta,norm=nTH2_norm))
		self.model_comps[2].color="C2"