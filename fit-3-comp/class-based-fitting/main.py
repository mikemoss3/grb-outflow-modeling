# Import custom classes
# from data_package import *
from model_package import *
from fit_package import *
from plot_package import *

if __name__ == "__main__":

	# from astropy.modeling.powerlaws import BrokenPowerLaw1D
	# from astropy.modeling import CompoundModel


	data = np.genfromtxt('../synthGRB_spec_TOT_2.txt',dtype=[("ENERGY",float),("RATE",float),("ERR",float)])
	data = data[data['ENERGY']<1e4]
	data['RATE']/=1e42
	data['ERR'] = np.sqrt(data['RATE'])
	data['RATE'] = np.random.normal(loc=data['RATE'],scale=np.sqrt(data['RATE'])) # Add Gaussian fluctuations

	# test = make_3comp(
	# TH_Tp=25,TH_alpha=0.4,TH_norm=1e6,
	# nTH1_e0=5e3,nTH1_alpha=-0.6,nTH1_beta=-2.5,nTH1_norm=1e5,
	# nTH2_e0=1e6,nTH2_alpha=-1.1,nTH2_beta=-2.5,nTH2_norm=1e5)

	test2 = Blackbody(temp=20,alpha=0.4,norm=5e5) + Band(e0=5e3,alpha=-1.5,beta=-2,norm=5e4)
	test2[0].color = "r"
	test2[0].alpha.fixed = True
	# test2[1].alpha.fixed = True
	# test2[1].beta.fixed = True


	# fitter = fitting.DogBoxLSQFitter()
	fitter = fitting.LMLSQFitter()
	# best_fit_model = fitter(test, data['ENERGY'],data['RATE'],maxiter=50000,weights=1.0/data['ERR'])
	best_fit_model = fitter(test2, data['ENERGY'],data['RATE'],maxiter=5000,weights=1.0/data['ERR'])

	print(best_fit_model)
	cov = fitter.fit_info['param_cov']
	print(dict(zip(test2.param_names, np.diag(cov)**0.5)))

	ax = plt.figure().gca()
	plot_data(data,ax=ax,spec_type=2)
	# plot_model(test2,ax=ax,spec_type=2)
	plot_model(best_fit_model,ax=ax,spec_type=2)

	plt.show()
