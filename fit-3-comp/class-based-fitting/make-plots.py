import numpy as np
import matplotlib.pyplot as plt

# Import custom classes
import Data
import Model
import Model3Comp

def plot_model(model : Model,ax=None):
	"""
	Method to plot model spectra 
	"""

	return 0;

def plot_data(data : Data,ax=None):
	"""
	Method to plot spectra data points
	"""

	return 0;


if __name__ == "__main__":

	# Load simulation data
	data_file_name = "synthetic_data.txt"
	data = Data(data_file_name)

	# Load Model
	model = Model3Comp()

	# Fit Model to Data 
	fit_package = FitPackage()
	fit_package.fit(model,data)


	# Plot Model and Data