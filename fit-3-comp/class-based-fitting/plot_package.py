import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import Fittable1DModel, CompoundModel

# Import custom classes
from data_package import *
from model_package import *
from fit_package import *

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

		for tick in ax.xaxis.get_minor_ticks():
			tick.label1.set_fontweight(fontweight)

			tick.label2.set_fontweight(fontweight)

	if yax is True:
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize=fontsize)
			tick.label1.set_fontweight(fontweight)

			tick.label2.set_fontsize(fontsize=fontsize)
			tick.label2.set_fontweight(fontweight)

		for tick in ax.yaxis.get_minor_ticks():
			tick.label1.set_fontweight(fontweight)

			tick.label2.set_fontweight(fontweight)
		
	ax.tick_params(direction="in",which="both")
	ax.margins(x=0,y=0)

def plot_model(model, spec_type = 0, inc_comps=True, emin=8, emax=4e4,ax=None,ymin=None,ymax=None,
	save_pref=None,xlabel=True,ylabel=True,label=None,fontsize=14,fontweight='bold',
	alpha=1,linestyle="solid",linewidth=1.5,
	comp_alpha=0.8,comp_linestyle="dashed",comp_linewidth=1.,en_steps = None,y_factor=1):
	"""
	Method to plot model spectra 

	spec_type
		0 == count spectrum
		1 == flux spectrum
		2 == power spectrum
	"""

	# Make plot instance if it doesn't exist
	if ax is None:
		ax = plt.figure().gca()
	# fig = plt.gcf()

	if en_steps is None:
		en_steps = int(np.log10(emax/emin)*20)
	energy_axis = np.logspace(np.log10(emin),np.log10(emax),num=en_steps)

	if spec_type == 1:
		y_factor *= energy_axis
	if spec_type == 2:
		y_factor *= energy_axis**2

	if isinstance(model, CompoundModel) and inc_comps:
		for i in range(model.n_submodels):
			ax.plot(energy_axis,y_factor*model[i](energy_axis),alpha=comp_alpha,color=model[i].color,linestyle=comp_linestyle,linewidth=comp_linewidth)

	ax.plot(energy_axis,y_factor*model(energy_axis), alpha=alpha,color=model.color,linestyle=linestyle,linewidth=linewidth)

	# Plot aesthetics
	ax.set_xscale('log')
	ax.set_yscale('log')
	if ymin is None:
		ymin = ax.get_ylim()[1] / 1e5
	ax.set_ylim(ymin,ymax)

	if xlabel is True:
		ax.set_xlabel('Energy (keV)',fontsize=fontsize,fontweight=fontweight)
	if ylabel is True:
		if spec_type == 0:
			ax.set_ylabel(r'N(E) (counts cm$^{-2}$ s$^{-1}$ keV$^{-1}$)',fontsize=fontsize,fontweight=fontweight)
		elif spec_type == 1:
			ax.set_ylabel(r'$F_{\nu}$ (erg$ cm$^{-2}$ s$^{-1}$ keV$^{-1}$)',fontsize=fontsize,fontweight=fontweight)
		elif spec_type == 2:
			ax.set_ylabel(r'$\nu F_{\nu}$ (erg$^2$ cm$^{-2}$ s$^{-1}$ keV$^{-1}$)',fontsize=fontsize,fontweight=fontweight)

	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)

	if label is not None:
		ax.legend(fontsize=fontsize-2)

	plt.tight_layout()

	if save_pref is not None :
		plt.savefig('{}.png'.format(save_pref))

	return 0;

def plot_data(data: Data, spec_type = 0, emin=None, emax=None,ax=None,ymin=None,ymax=None,
	save_pref=None,xlabel=True,ylabel=True,label=None,fontsize=14,fontweight='bold',
	color="C0",marker=".",alpha=1,linestyle="solid",linewidth=1.5,y_factor=1):
	"""
	Method to plot observed spectra 

	spec_type
		0 == count spectrum
		1 == flux spectrum
		2 == power spectrum
	"""

	# Make plot instance if it doesn't exist
	if ax is None:
		ax = plt.figure().gca()
	# fig = plt.gcf()

	if spec_type == 1:
		y_factor *= data['ENERGY']
	if spec_type == 2:
		y_factor *= data['ENERGY']**2

	ax.errorbar(data['ENERGY'],y_factor*data['RATE'], yerr = y_factor*data['ERR'],color=color,fmt=" ",alpha=alpha)

	# ax[1].errorbar(x=data['ENERGY'], y=resids2, yerr = 1 ,fmt=" ",color="C1",alpha=0.9,marker="+",zorder=0)

	# Plot aesthetics
	ax.set_xscale('log')
	ax.set_yscale('log')
	if ymin is None:
		ymin = ax.get_ylim()[1] / 1e5
	ax.set_ylim(ymin,ymax)
	ax.set_xlim(emin,emax)

	if xlabel is True:
		ax.set_xlabel('Energy (keV)',fontsize=fontsize,fontweight=fontweight)
	if ylabel is True:
		if spec_type == 0:
			ax.set_ylabel(r'N(E) (counts cm$^{-2}$ s$^{-1}$ keV$^{-1}$)',fontsize=fontsize,fontweight=fontweight)
		elif spec_type == 1:
			ax.set_ylabel(r'$F_{\nu}$ (erg$ cm$^{-2}$ s$^{-1}$ keV$^{-1}$)',fontsize=fontsize,fontweight=fontweight)
		elif spec_type == 2:
			ax.set_ylabel(r'$\nu F_{\nu}$ (erg$^2$ cm$^{-2}$ s$^{-1}$ keV$^{-1}$)',fontsize=fontsize,fontweight=fontweight)

	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)

	if label is not None:
		ax.legend(fontsize=fontsize-2)

	plt.tight_layout()

	if save_pref is not None :
		plt.savefig('{}.png'.format(save_pref))

	return 0;
