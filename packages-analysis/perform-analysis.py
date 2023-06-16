# Import custom classes
# from data_package import *
from model_package import *
from fit_package import *
from plot_package import *

if __name__ == "__main__":

	fn = []
	fn.append('../files-data/synthetic-data/synthGRB_spec_TOT_0.txt')
	fn.append('../files-data/synthetic-data/synthGRB_spec_TOT_1.txt')
	fn.append('../files-data/synthetic-data/synthGRB_spec_TOT_2.txt')
	
	tstart = [0,1,2]
	tend = [1,2,3]
	colors= ["C0","C1","C2"]

	data = []
	data_inst = Data()
	ax = plt.figure().gca()
	for i in range(len(fn)):
		
		data.append(np.genfromtxt(fn[i],dtype=[("ENERGY",float),("RATE",float),("ERR",float)]))
		data[i]['ERR'] = np.sqrt(data[i]['RATE'])
		np.savetxt(fn[i],data[i])

		data_inst.load_spectrum(fn[i],tstart[i],tend[i])

		# plot_data_spec(data_inst.spectra[i]['SPECTRUM'],ax=ax,color=colors[i],spec_type=2)

		# model = Blackbody(temp=40,alpha=0.4,norm=3e3) + Band(e0=5e2,alpha=-1.1,beta=-2.05,norm=10812.635)
		# model = Blackbody(temp=20,alpha=0.4,norm=5e5) + Band(e0=5e3,alpha=-1.1,beta=-2.5,norm=5e4)
		model = Blackbody(temp=35,alpha=0.4,norm=3e4) + Band(e0=1e3,alpha=-1.1,beta=-2.3,norm=3e4)
		model[0].alpha.fixed = True
		model[0].color = colors[i]
		model[1].color = "white"
		model.color=colors[i]
		fitter = FittedModel()
		best_fit_model, fitstat = fitter.fit(model, data_inst.spectra[i]['SPECTRUM'],verbose=True)

		plot_data_spec(data_inst.spectra[i]['SPECTRUM'],ax=ax,color=colors[i],spec_type=2)
		plot_model_spec(best_fit_model,ax=ax,inc_comps=True,spec_type=2)

	plt.show()	
