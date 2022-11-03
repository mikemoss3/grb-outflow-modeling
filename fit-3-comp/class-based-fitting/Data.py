import numpy as np

# File used to define Data class

class Data(object):
	def __init__(self,file_name = None):
		if file_name is None:
			self.data=None
		else:
			self.load_data(file_name)

	def load_data(self,file_name):
		self.data = np.genfromtxt(file_name,dtype=[('ENERGY',float),('RATE',float),('ERROR',float)])



	
