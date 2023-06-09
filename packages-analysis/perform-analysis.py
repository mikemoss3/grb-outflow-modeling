# Import custom classes
# from data_package import *
from packages.model_package import *
from packages.fit_package import *
from packages.plot_package import *

if __name__ == "__main__":

	# from astropy.modeling.powerlaws import BrokenPowerLaw1D
	# from astropy.modeling import CompoundModel


	data = np.genfromtxt('data-files/synthGRB_spec_TOT_0.txt',dtype=[("ENERGY",float),("RATE",float),("ERR",float)])
	data = data[data['ENERGY']<1e4]
	data['RATE'] /= np.min(data['RATE'])
	data['ERR'] = np.sqrt(data['RATE'])*3
	data['RATE'] = np.random.normal(loc=data['RATE'],scale=np.sqrt(data['RATE'])) # Add Gaussian fluctuations
	data['RATE'][data['RATE']<0] = 0

	# test = make_3comp(
	# TH_Tp=20,TH_alpha=0.4,TH_norm=1e5,
	# nTH1_e0=5e4,nTH1_alpha=-0.7,nTH1_beta=-8.5,nTH1_norm=1e3,
	# nTH2_e0=1e7,nTH2_alpha=-1.5,nTH2_beta=-8.5,nTH2_norm=1e3)

	# test = Blackbody(temp=20,alpha=0.4,norm=5e5) + Band(e0=5e3,alpha=-1.1,beta=-2.5,norm=5e4)
	test = Blackbody(temp=40,alpha=0.4,norm=3e3) + Band(e0=5e2,alpha=-1.1,beta=-2.05,norm=10812.635)
	test[0].color = "r"
	test[0].alpha.fixed = True
	# test[1].alpha.fixed = True
	# test[1].beta.fixed = True

	# test = Blackbody(temp=20,alpha=0.4,norm=5e5)
	# test = Band(e0=5e3,alpha=-1.5,beta=-2,norm=5e4)

	best_fit = FittedModel()
	model, fitstat = best_fit.fit(test, data,verbose=True)

	ax = plt.figure().gca()
	plot_data(data,ax=ax,spec_type=0)
	# plot_model(test,ax=ax,spec_type=0)
	plot_model(model,ax=ax,spec_type=0)

	plt.show()	
