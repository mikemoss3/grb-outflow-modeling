import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import warnings
# import Fitting

#suppress warnings
warnings.filterwarnings('ignore')

def plot_aesthetics(ax,fontsize=14,fontweight='bold',xax=True,yax=True):
	"""
	This function is used to make bold and increase the font size of all plot tick markers
	"""

	if xax is True:
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize=fontsize)
			tick.label1.set_fontweight(fontweight)

			tick.label2.set_fontsize(fontsize=fontsize)
			tick.label2.set_fontweight(fontweight)

	if yax is True:
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize=fontsize)
			tick.label1.set_fontweight(fontweight)

			tick.label2.set_fontsize(fontsize=fontsize)
			tick.label2.set_fontweight(fontweight)

def thermal(energy,TH_kT,TH_norm):
	TH_alpha = 0.4

	if hasattr(energy,"__len__"):
		flux_value = np.zeros(shape=len(energy))
		
		for i in range(len(energy)):
			if energy[i] < 2*1e3:
				flux_value[i] = TH_norm * np.power(energy[i]/TH_kT, 1+TH_alpha) / (np.exp(energy[i]/TH_kT) - 1.)	
			else:
				flux_value[i] = 0.

	else:
		flux_value = 0.

		if energy < 2*1e3:
			flux_value = TH_norm * np.power(energy/TH_kT, 1+TH_alpha) / (np.exp(energy/TH_kT) - 1.)	

		else:
			flux_value = 0.

	return flux_value

def CPL(energy,cpl_ep,cpl_alpha,cpl_norm):
	if hasattr(energy,"__len__"):
		flux_value = np.zeros(shape=len(energy))
		
		for i in range(len(energy)):
			flux_value[i] = cpl_norm * np.power(energy[i]/cpl_ep, cpl_alpha) * np.exp(-energy[i]/cpl_ep)

	else:
		flux_value = 0.

		flux_value = cpl_norm * np.power(energy/cpl_ep, cpl_alpha) * np.exp(-energy/cpl_ep)

	return flux_value

def non_thermal(energy,nTH_e0,nTH_alpha,nTH_beta,nTH_norm):
	"""
	Band function
	"""

	if hasattr(energy,"__len__"):
		flux_value = np.zeros(shape=len(energy))

		for i in range(len(energy)):
			if energy[i] < (nTH_alpha - nTH_beta)*nTH_e0 / ( 2 + nTH_alpha):
				flux_value[i] += nTH_norm * np.power(energy[i]/100., nTH_alpha) * np.exp(- energy[i]*( 2 + nTH_alpha) / nTH_e0)
			else:
				flux_value[i] += nTH_norm * np.power((nTH_alpha - nTH_beta) * nTH_e0/100./ ( 2 + nTH_alpha), nTH_alpha - nTH_beta) * np.exp(nTH_beta - nTH_alpha) * np.power(energy[i] /100,nTH_beta)

	else:
		flux_value = 0.

		if energy < (nTH_alpha - nTH_beta)*nTH_e0/( 2 + nTH_alpha):
			flux_value = nTH_norm * np.power(energy/100., nTH_alpha) * np.exp(- energy*( 2 + nTH_alpha) / nTH_e0)
		else:
			flux_value = nTH_norm * np.power((nTH_alpha - nTH_beta) * nTH_e0/100./( 2 + nTH_alpha), nTH_alpha - nTH_beta) * np.exp(nTH_beta - nTH_alpha) * np.power(energy /100,nTH_beta)

	return flux_value

def model_2comp(energy,
	TH_kT = 10,
	# TH_alpha = 0.4,
	TH_norm = 1.,
	nTH1_e0 = 1e4,
	nTH1_alpha = -0.6, 
	# nTH1_beta = -2.5,
	nTH1_norm = 1,
	):
	"""
	Evaluate the three component model at a specific energy
	"""

	# If energy is an array 
	if hasattr(energy,"__len__"):
		flux_value = np.zeros(shape=len(energy))

		for i in range(len(energy)):

			tmp_en = energy[i]

			# Thermal component contribution (if less than 2 MeV)
			flux_value[i] += thermal(tmp_en,TH_kT,TH_norm)

			# First non-thermal component contribution
			# flux_value[i] += non_thermal(tmp_en,nTH1_e0,nTH1_alpha,nTH1_beta,nTH1_norm)
			flux_value[i] += CPL(tmp_en,nTH1_e0,nTH1_alpha,nTH1_norm)

	else:
		flux_value = 0. 

		# Thermal component contribution (if less than 2 MeV)
		flux_value += thermal(energy,TH_kT,TH_norm)

		# First non-thermal component contribution
		# flux_value += non_thermal(energy,nTH1_e0,nTH1_alpha,nTH1_beta,nTH1_norm)
		flux_value += CPL(energy,nTH1_e0,nTH1_alpha,nTH1_norm)

	return flux_value

# make three component model
def model_3comp(energy,
	TH_kT = 10,
	# TH_alpha = 0.4,
	TH_norm = 1.,
	nTH1_e0 = 1e4,
	nTH1_alpha = -0.6, 
	# nTH1_beta = -2.5,
	nTH1_norm = 1,
	nTH2_e0 = 1e7,
	nTH2_alpha = -1.1, 
	# nTH2_beta = -2.5,
	nTH2_norm = 1
	):
	"""
	Evaluate the three component model at a specific energy
	"""

	# If energy is an array 
	if hasattr(energy,"__len__"):
		flux_value = np.zeros(shape=len(energy))

		for i in range(len(energy)):

			tmp_en = energy[i]

			# Thermal component contribution (if less than 2 MeV)
			flux_value[i] += thermal(tmp_en,TH_kT,TH_norm)

			# First non-thermal component contribution
			flux_value[i] += CPL(tmp_en,nTH1_e0,nTH1_alpha,nTH1_norm)

			# Second non-thermal component contribution
			flux_value[i] += CPL(tmp_en,nTH2_e0,nTH2_alpha,nTH2_norm)

	else:
		flux_value = 0. 

		# Thermal component contribution (if less than 2 MeV)
		flux_value += thermal(energy,TH_kT,TH_norm)

		# First non-thermal component contribution
		flux_value += CPL(energy,nTH1_e0,nTH1_alpha,nTH1_norm)

		# Second non-thermal component contribution
		flux_value += CPL(energy,nTH2_e0,nTH2_alpha,nTH2_norm)

	return flux_value

# Plot results
fig, ax = plt.subplots(2,3,figsize=(18,6),sharex=True,gridspec_kw={'height_ratios': [3, 1]})

file_name_arr = np.array(["synthGRB_spec_TOT_0.txt", "synthGRB_spec_TOT_1.txt", "synthGRB_spec_TOT_2.txt"], dtype="U50")

fontsize=18
fontweight="bold"


for i in range(3):
	# Load synthetic data
	file_name = file_name_arr[i]
	data = np.genfromtxt(file_name,dtype=[("ENERGY",np.float32),("RATE",np.float64),("ERR",np.float64)])
	data['ENERGY'] = data['ENERGY']/2
	data = data[data['ENERGY']<8e3]


	data['RATE'] /= 5e42
	data['ERR'] = np.sqrt(data['RATE'])*5
	# data['ERR'] = data['RATE']/10

	# Make it so data can be put into float32 type
	data['RATE'].astype(np.float32)
	data['ERR'].astype(np.float32)

	# Use every 5 data points 
	data = data[::10]

	# Perform fit

	# 1-comp
	# popt1, pcov1 = opt.curve_fit(f=non_thermal ,xdata=data['ENERGY'],ydata=data['RATE'],sigma=data['ERR'],p0=np.array([ 3000,-1.6, -2.1,1e3]))
	popt1, pcov1 = opt.curve_fit(f=non_thermal ,xdata=data['ENERGY'],ydata=data['RATE'],sigma=data['ERR'],p0=np.array([ 3000,-1.6, -2.1,1e3]),
		bounds= ((0.1,-5,-8,1e-3),(1e7,0,0,1e10) ))
	resids1 = (non_thermal(data['ENERGY'],*popt1) - data['RATE']) / data['ERR']

	# 2-comp
	# Unconstrained, uses LM fitting method
	# popt2, pcov2 = opt.curve_fit(f=model_2comp ,xdata=data['ENERGY'],ydata=data['RATE'],sigma=data['ERR'],p0=np.array([ 20, 1, 3000,-1.6,1]))
	# Constrained, uses 'trf' method. 
	popt2, pcov2 = opt.curve_fit(f=model_2comp ,xdata=data['ENERGY'],ydata=data['RATE'],sigma=data['ERR'],p0=np.array([10,1e3, 3e3,-0.7,10]),
		bounds= ((0.1,1e-3,  1,-5,1e-3),(1e3,1e5,  1e7,0,1e3) ))
	resids2 =(model_2comp(data['ENERGY'],*popt2) - data['RATE']) / data['ERR']

	# popt2, pcov2 = opt.curve_fit(f=model_2comp ,xdata=data['ENERGY'],ydata=data['RATE'],sigma=data['ERR'],p0=np.array([10,1e3, 3e3,-1,-2.1,1e3]),
	# 	bounds= ((1,1e-3,  2e1,-5,-5,1e-3),(1e3,1e5,  1e7,0,0,1e5) ))
	# resids2 = (model_2comp(data['ENERGY'],*popt2) - data['RATE']) / data['ERR']


	# 3-comp 
	# popt3, pcov3 = opt.curve_fit(f=model_3comp ,xdata=data['ENERGY'],ydata=data['RATE'],sigma=data['ERR'],p0=np.array([10,1e3, 3e3,-0.7,20, 3e5,-1.,20]),maxfev=100000)
	# popt3, pcov3 = opt.curve_fit(f=model_3comp ,xdata=data['ENERGY'],ydata=data['RATE'],sigma=data['ERR'],p0=np.array([10,1e3, 3e3,-0.7,20, 3e5,-1.,20]), 
		# bounds= ((1,1e-3,  2e3,-5,1e-3,  2e3,-5,1e-3),(1e3,1e5,  1e7,0,1e3,  1e7,0,1e3) ))
	# resids3 = (model_3comp(data['ENERGY'],*popt3) - data['RATE'])


	# print(*popt1)
	# print(*popt2)
	# print(*popt3)
	# print()

	energy_axis = np.logspace(np.log10(data['ENERGY'][0]),np.log10(data['ENERGY'][-1]))


	ax[0,i].scatter(data['ENERGY'],data['ENERGY']**2 * data['RATE'],marker=".",color="k",alpha=0.3)
	ax[0,i].errorbar(data['ENERGY'],data['ENERGY']**2 * data['RATE'],data['ENERGY']**2 * data['ERR'] ,fmt=" ",color="k",alpha=0.3)

	ax[0,i].plot(energy_axis, np.power(energy_axis,2.) * thermal(energy_axis,*popt2[:2]),color="r",linestyle="dashed",alpha=0.9)
	ax[0,i].plot(energy_axis, np.power(energy_axis,2.) * CPL(energy_axis,*popt2[2:]),color="C0",linestyle="dashed",alpha=0.9)
	# ax[0,i].plot(energy_axis, np.power(energy_axis,2.) * non_thermal(energy_axis,*popt2[2:]),color="C0",linestyle="dashed",alpha=0.9)

	ax[0,i].plot(energy_axis, np.power(energy_axis,2.) * non_thermal(energy_axis,*popt1),color="purple",alpha=0.7)
	ax[0,i].plot(energy_axis, np.power(energy_axis,2.) * model_2comp(energy_axis,*popt2),color="C1")
	# ax[0].plot(energy_axis, np.power(energy_axis,2.) * model_3comp(energy_axis,*popt3),color="C2")

	# xerr = (data['ENERGY'][1:] - data['ENERGY'][0:-1])
	ax[1,i].errorbar(x=data['ENERGY'], y=resids2, yerr = 1 ,fmt=" ",color="C1",alpha=0.9,marker="+",zorder=0)
	# ax[1,i].errorbar(x=data['ENERGY'], y=resids3, yerr = 1 ,fmt=" ",color="C2",alpha=0.7,marker="+",zorder=0)
	ax[1,i].errorbar(x=data['ENERGY'], y=resids1, yerr = 1 ,fmt=" ",color="purple",alpha=0.7,marker="+",zorder=4)

	ax[1,i].hlines(y=0,xmin=energy_axis[0], xmax=energy_axis[-1],alpha=0.5)

	# Plot aesthetics 
	ax[0,i].set_yscale('log')
	ax[0,i].set_xscale('log')
	ax[1,i].set_xscale('log')

	ax[0,i].set_ylim(2e5,2e7)
	ax[1,i].set_ylim(-3,3)


	if i == 0:
		ax[0,i].set_ylabel("Rate (Arb. Units)",fontsize=fontsize,fontweight=fontweight)
		ax[1,i].set_ylabel("Sigma",fontsize=fontsize,fontweight=fontweight)
	if i == 1:
		ax[1,i].set_xlabel("Energy (keV)",fontsize=fontsize,fontweight=fontweight)

	if i > 0:
		ax[0,i].set_yticklabels([])
		ax[1,i].set_yticklabels([])

		plot_aesthetics(ax[0,i],fontsize=fontsize,fontweight=fontweight,yax=False,xax=False)
		plot_aesthetics(ax[1,i],fontsize=fontsize,fontweight=fontweight,yax=False)
	else:
		plot_aesthetics(ax[0,i],fontsize=fontsize,fontweight=fontweight,xax=False)
		plot_aesthetics(ax[1,i],fontsize=fontsize,fontweight=fontweight)


plt.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)

plt.savefig("fitted-synth-data-v02.png")

plt.show()


