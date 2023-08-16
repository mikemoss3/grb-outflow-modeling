# Import custom classes
# from packages-analysis.data_package import *
from packages_analysis.model_package import *
from packages_analysis.fit_package import *
from packages_analysis.plot_package import *

# from asymmetric_uncertainty import a_u
# import astropy

if __name__ == "__main__":
	"""
	#### Generating synthetic data from empirical models 
	input_model = astropy.modeling.powerlaws.PowerLaw1D(alpha=0.2,amplitude=1000)
	input_model.x_0.fixed = True
	num = 120 # Number of data points 
	energy_points = np.logspace(0,4,num) # Make energy axis 
	tmp_rate_data = input_model(energy_points) # Evaluate the model at each energy 
	tmp_rate_data = np.random.normal(loc=tmp_rate_data,scale=np.sqrt(tmp_rate_data)) # Add Gaussian fluctuations
	tmp_rate_data[tmp_rate_data<0] *=-1 # Remove negative fluctuations

	# Initialize the data array to store data 
	input_array = np.zeros(shape=num,dtype=[("ENERGY",float),("RATE",a_u)])
	input_array['ENERGY'] = energy_points # set energy axis
	input_array['RATE'] = [a_u(tmp_rate_data[i],np.sqrt(tmp_rate_data[i]),np.sqrt(tmp_rate_data[i])) for i in range(len(tmp_rate_data)) ]

	## Fitting non-linear data

	# Initialize model
	model = astropy.modeling.powerlaws.PowerLaw1D(alpha=0.5,amplitude=1000)
	model.x_0.fixed = True

	# Initialize the fitter
	fitter = FittedModel()

	# Perform the fit. By specifying "verbose=True", the best fit parameter and reduced chi-squared will be printed to the terminal
	best_fit_model, fitstat = fitter.fit(model, input_array,verbose=True)
	# best_fit_model, fitstat = fitter.fit_mc(model, input_array,verbose=True,iterations=50)

	# Make fit spectrum 
	best_fit_model_spec = np.zeros(shape=num,dtype=[("ENERGY",float),("RATE",float),('ERR',float)])
	best_fit_model_spec['ENERGY'] = energy_points
	best_fit_model_spec['RATE'] = best_fit_model.amplitude * np.power(best_fit_model_spec['ENERGY'],-best_fit_model.alpha)


	## Plot
	ax = plt.figure().gca()

	# Plot the input model used the generate the data
	# plot_model_spec(input_model,emin=input_array['ENERGY'][0],emax=input_array['ENERGY'][-1], ax=ax,inc_comps=True,spec_type=0,linestyle="dashdot",comp_linestyle="dashdot",alpha=0.6, label="Input Model")
	# Plot the data 
	plot_data_spec(input_array,ax=ax,spec_type=0,alpha=0.5)
	# Plot the best-fit model
	# plot_data_spec(best_fit_model_spec, ax=ax,spec_type=0,linewidth=3,joined=True,unc=False,label="Best Fit")

	plt.show()
	"""
	
	### Trying to fit data generated by the simulation code 


	fn = []
	fn.append('files-data/synthetic-data/synthGRB_spec_TOT_0.txt')
	fn.append('files-data/synthetic-data/synthGRB_spec_TOT_1.txt')
	fn.append('files-data/synthetic-data/synthGRB_spec_TOT_2.txt')
	
	tstart = [0,1,2]
	tend = [1,2,3]
	colors= ["C0","C1","C2"]

	data = []
	data_inst = Data()
	
	# ax = plt.figure().gca()
	fig, ax = plt.subplots(2,1,figsize=(6,8),sharex=True,gridspec_kw={'hspace':0.,'height_ratios':[2,1]})
	
	for i in range(len(fn)):
		data.append(np.genfromtxt(fn[i],dtype=[("ENERGY",float),("RATE",float),("ERR",float)]))
		# np.savetxt(fn[i],data[i])

		data_inst.load_spectrum(fn[i],tstart[i],tend[i])

		# plot_data_spec(data_inst.spectra[i]['SPECTRUM'],ax=ax,color=colors[i],spec_type=2)

		# model = Blackbody(temp=40,alpha=0.4,norm=3e3) + Band(e0=5e2,alpha=-1.1,beta=-2.05,norm=10812.635)
		model = Blackbody(temp=20,alpha=0.4,norm=3e4) + Band(e0=1e1,alpha=-1.,beta=-2.,norm=3e4)
		model[0].alpha.fixed = True
		model[0].color = colors[i]
		model[1].color = "white"
		# model[1].alpha.fixed = True
		# model[1].beta.fixed = True


		model.color=colors[i]
		fitter = FittedModel()
		# best_fit_model, fitstat = fitter.fit(model, data_inst.spectra[i]['SPECTRUM'],verbose=True)
		best_fit_model, fitstat = fitter.fit_mc(model, data_inst.spectra[i]['SPECTRUM'],verbose=True,iterations=10)

		# plot_data_spec(data_inst.spectra[i]['SPECTRUM'],ax=ax,color=colors[i],spec_type=2)
		# plot_model_spec(best_fit_model,ax=ax,inc_comps=True,spec_type=2)

		# plot_spec(spectrum_data = data_inst.spectra[i]['SPECTRUM'],spec_type=2,color=colors[i],ax=ax[0],xlabel=False)
		# plot_spec(spectrum_model=best_fit_model,spec_type=2,ax=ax[0],xlabel=False)
		plot_spec(spectrum_model=best_fit_model,spectrum_data = data_inst.spectra[i]['SPECTRUM'],plot_res=True,color=colors[i],spec_type=2,ax=ax)


	plt.show()