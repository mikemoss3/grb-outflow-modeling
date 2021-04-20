"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-02-15


Meta script to plot desired simulation results

"""
import matplotlib.pyplot as plt
import numpy as np
import os
import shelldists as sd
import radiation as rd


def plot_aesthetics(ax,fontsize=14,fontweight='bold'):

	for tick in ax.xaxis.get_major_ticks():
	    tick.label1.set_fontsize(fontsize=fontsize)
	    tick.label1.set_fontweight(fontweight)
	for tick in ax.yaxis.get_major_ticks():
	    tick.label1.set_fontsize(fontsize=fontsize)
	    tick.label1.set_fontweight(fontweight)

def plot_lor_dist(save=False):
	"""
	Plot Lorentz distributions
	"""
	# Load data
	t0=np.genfromtxt('sim_results/t0_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
	if os.path.exists('sim_results/t3e4_shells.txt'):
		t3e4=np.genfromtxt('sim_results/t3e4_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
	if os.path.exists('sim_results/t2e5_shells.txt'):
		t2e5=np.genfromtxt('sim_results/t2e5_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
	if os.path.exists('sim_results/t5e5_shells.txt'):
		t5e5=np.genfromtxt('sim_results/t5e5_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])

	fig = plt.figure()
	ax = fig.gca()
	sd.plot_lorentz_dist(ax,t0[t0['STATUS']!=0], label='t0')
	if os.path.exists('sim_results/t3e4_shells.txt'):
		sd.plot_lorentz_dist(ax,t3e4[t3e4['STATUS']!=0], label='t3e4, s#={}'.format(len(t3e4[t3e4['STATUS']==1])))
	if os.path.exists('sim_results/t2e5_shells.txt'):
		sd.plot_lorentz_dist(ax,t2e5[t2e5['STATUS']!=0], label='t2e5, s#={}'.format(len(t2e5[t2e5['STATUS']==1])))
	if os.path.exists('sim_results/t5e5_shells.txt'):
		sd.plot_lorentz_dist(ax,t5e5[t5e5['STATUS']!=0], label='t5e5, s#={}'.format(len(t5e5[t5e5['STATUS']==1])))
	plt.legend()
	
	if save == True:
		plt.savefig('sim_results/lorentz-dist.png')


def plot_spec(ax, emission,z=0,comp_num=None, plot_comps=True, nuFnu=True, Tmin=None, Tmax=None, Emin=0.1, Emax=5e5, num_bins=60,save=False):
	"""
	Plot Spectrum and components
	"""

	# fig = plt.figure()
	# ax = fig.gca()
	emission.plot_spectrum(ax,comp_num=comp_num,plot_comps=plot_comps, nuFnu=nuFnu, Tmin=Tmin, Tmax=Tmax, Emin=Emin, Emax=Emax, num_bins=num_bins,z=z)

	# ax.vlines(x=4,ymin=1e46,ymax=1e53,linestyle='dashed',color='black',label='GBM Band')
	# ax.vlines(x=4*1e4,ymin=1e46,ymax=1e53,linestyle='dashed',color='black')

	# Plot aesthetics

	ax.set_xscale('log')
	ax.set_yscale('log')

	# For axis labels
	fontsize=14
	fontweight='bold'

	ax.set_xlabel('E (keV)',fontsize=fontsize,fontweight=fontweight)

	ax.set_ylim(1e48)
	ax.set_xlim(Emin,Emax)

	plot_aesthetics(ax)
	ax.legend()
	# fig.tight_layout()
	
	if save == True:
		plt.savefig('sim_results/spectrum.png')

def plot_light_curve(ax,emission,z=0,comp_num=None, plot_comps=True, Tmin=None, Tmax=None,dt=0.05, Emin=0.1, Emax=5e5,save=False):
	"""
	Plot light curve per time bin 
	"""

	emission.plot_light_curve(ax, comp_num=comp_num,plot_comps=plot_comps, Tmin=Tmin, Tmax=Tmax, dt=dt, Emin=Emin, Emax=Emax, z=z)

	# Plot aesthetics
	ax.set_xlim(Tmin,Tmax)

	# For axis labels
	fontsize=14
	fontweight='bold'

	ax.set_ylabel('Rate (cts/sec)',fontsize=fontsize,fontweight=fontweight)
	ax.set_xlabel('Time (sec)',fontsize=fontsize,fontweight=fontweight)

	plot_aesthetics(ax)
	ax.legend()
	fig.tight_layout()
	
	if save == True:
		plt.savefig('sim_results/lightcurve.png')

def plot_param_vs_t(param,z,save=False):
	"""
	Plot Spectrum and components
	"""

	# Load spectrum data 
	spectrum = rd.Emission()
	spectrum.load_therm_spec('sim_results/ordlor_spectrum_therm.txt',z)
	spectrum.load_synch_spec('sim_results/ordlor_spectrum_synch.txt',z)

	spec_therm = spectrum.spec_therm
	spec_synch = spectrum.spec_synch


	fig = plt.figure()
	ax = fig.gca()

	if param == 'T' or param =='L':
		ax.scatter(spec_therm['te'],spec_therm[param],marker='.',label=param)
	else:
		ax.scatter(spec_synch['te'],spec_synch[param],marker='.',label=param)

	fig.legend()

	fontsize=14
	fontweight='bold'

	ax.set_ylabel(param,fontsize=fontsize,fontweight=fontweight)
	ax.set_xlabel(r't$_e$',fontsize=fontsize,fontweight=fontweight)

	for tick in ax.xaxis.get_major_ticks():
	    tick.label1.set_fontsize(fontsize=fontsize)
	    tick.label1.set_fontweight(fontweight)
	for tick in ax.yaxis.get_major_ticks():
	    tick.label1.set_fontsize(fontsize=fontsize)
	    tick.label1.set_fontweight(fontweight)

	fig.tight_layout()
		
	if save == True:
		plt.savefig('sim_results/param-{}-vs-t.png'.format(param))

def plot_evo_therm(z,save=False):
	"""
	Plot evolution of thermal emission parameters 
	"""
	# Load spectrum data 
	spectrum = rd.Emission()
	spectrum.load_therm_spec('sim_results/ordlor_spectrum_therm.txt',z)
	spec_therm = spectrum.spec_therm

	fontsize=14
	fontweight='bold'

	fig = plt.figure()
	ax = fig.gca()
	boltz= 8.617*1e-8 # keV K^-1 
	ax.scatter(spec_therm['ta'],spec_therm['T']*boltz)
	ax.set_ylabel(r'k$_B$T (KeV)',fontsize=fontsize,fontweight=fontweight)

	# ax.scatter(spec_therm['ta'],spec_therm['L'])
	# plt.set_ylabel('L (erg/s)',fontsize=fontsize,fontweight=fontweight)

	ax.set_xlabel('Arrival Time',fontsize=fontsize,fontweight=fontweight)

	ax.set_yscale('log')

	for tick in ax.xaxis.get_major_ticks():
	    tick.label1.set_fontsize(fontsize=fontsize)
	    tick.label1.set_fontweight(fontweight)
	for tick in ax.yaxis.get_major_ticks():
	    tick.label1.set_fontsize(fontsize=fontsize)
	    tick.label1.set_fontweight(fontweight)

	
	if save == True:
		plt.savefig('sim_results/thermal-evo.png')

def plot_evo_synch(z,dt=0.1,save=False):
	"""
	Plot evolution of synchrotron emission parameters 
	"""

	# Load spectrum data 
	spectrum = rd.Emission()
	spectrum.load_synch_spec('sim_results/ordlor_spectrum_synch.txt',z)
	spec_synch = spectrum.spec_synch

	fontsize=14
	fontweight='bold'


	fig0, ax = plt.subplots(2,1,figsize=(7,8),sharex=True)
	
	ax0cp = ax[0].twinx()
	ax[0].scatter(spec_synch['ta']*(1+z),spec_synch['te'],marker='^')
	ax0cp.scatter(spec_synch['ta']*(1+z),spec_synch['delt'],marker='.',color='r')
	ax[0].set_ylabel(r'$t_{e}$',fontsize=fontsize,fontweight=fontweight)
	ax0cp.set_ylabel(r'$\Delta t$',fontsize=fontsize,fontweight=fontweight)
	ax0cp.yaxis.set_label_position("right")
	ax0cp.yaxis.tick_right()


	ax1cp = ax[1].twinx()
	ax[1].scatter(spec_synch['ta']*(1+z),spec_synch['e']/spec_synch['e'][0],marker='^')

	ax1cp.scatter(spec_synch['ta']*(1+z),spec_synch['gammar']/100,marker='.',color='r')
	ax[1].set_ylabel('e (en. diss.)',fontsize=fontsize,fontweight=fontweight)
	ax1cp.set_ylabel(r'$\Gamma_{r}/100$',fontsize=fontsize,fontweight=fontweight)
	ax1cp.yaxis.set_label_position("right")
	ax1cp.yaxis.tick_right()
	
	ax[0].set_ylim(-0.25*10**5,6.25*10**5)
	ax0cp.set_ylim(-0.2,4.2)
	ax[1].set_xlim(-1)
	ax1cp.set_ylim(-0.2,4.2)

	ax[0].grid(axis='x')
	ax[1].grid(axis='x')

	for i in range(2):
		for tick in ax[i].xaxis.get_major_ticks():
		    tick.label1.set_fontsize(fontsize=fontsize)
		    tick.label1.set_fontweight(fontweight)
		for tick in ax[i].yaxis.get_major_ticks():
		    tick.label1.set_fontsize(fontsize=fontsize)
		    tick.label1.set_fontweight(fontweight)

	for twin in [ax0cp,ax1cp]:
		for tick in twin.xaxis.get_major_ticks():
		    tick.label1.set_fontsize(fontsize=fontsize)
		    tick.label1.set_fontweight(fontweight)
		for tick in twin.yaxis.get_major_ticks():
		    tick.label2.set_fontsize(fontsize=fontsize)
		    tick.label2.set_fontweight(fontweight)

	plt.tight_layout()
	plt.subplots_adjust(hspace=0)

	if save == True:
		plt.savefig('sim_results/synch-evo-fig0.png')

	fig1, ax = plt.subplots(2,2,sharex=True,figsize=(12,8))

	ax[0,0].scatter(spec_synch['ta']*(1+z),spec_synch['asyn'])
	ax[0,0].set_ylabel(r'$\alpha_{syn}$',fontsize=fontsize,fontweight=fontweight)
	# ax[0,0].set_xlabel(r't$_a$ (sec), Arrival Time',fontsize=fontsize,fontweight=fontweight)
	
	ax[0,1].scatter(spec_synch['ta']*(1+z),spec_synch['gammae']/1e4)
	ax[0,1].set_ylabel(r'$\Gamma_{e}$/1e4',fontsize=fontsize,fontweight=fontweight)
	# ax[0,1].set_xlabel(r't$_a$ (sec), Arrival Time',fontsize=fontsize,fontweight=fontweight)
	ax[0,1].yaxis.set_label_position("right")
	ax[0,1].yaxis.tick_right()

	ax[1,0].scatter(spec_synch['ta']*(1+z),spec_synch['Beq'],marker='.')
	ax[1,0].set_yscale('log')
	ax[1,0].set_ylabel(r'B$_{eq}$',fontsize=fontsize,fontweight=fontweight)
	ax[1,0].set_xlabel(r't$_a$ (sec), Arrival Time',fontsize=fontsize,fontweight=fontweight)


	ax[1,1].scatter(spec_synch['ta']*(1+z),spec_synch['Esyn'])
	ax[1,1].set_yscale('log')
	ax[1,1].set_ylabel(r'$E_{syn}$/1e3',fontsize=fontsize,fontweight=fontweight)
	ax[1,1].set_xlabel(r't$_a$ (sec), Arrival Time',fontsize=fontsize,fontweight=fontweight)
	ax[1,1].yaxis.set_label_position("right")
	ax[1,1].yaxis.tick_right()

	for i in range(2):
		for j in range(2):
			for tick in ax[i,j].xaxis.get_major_ticks():
			    tick.label1.set_fontsize(fontsize=fontsize)
			    tick.label1.set_fontweight(fontweight)
			for tick in ax[i,j].yaxis.get_major_ticks():
			    tick.label1.set_fontsize(fontsize=fontsize)
			    tick.label1.set_fontweight(fontweight)
			for tick in ax[i,j].yaxis.get_major_ticks():
			    tick.label2.set_fontsize(fontsize=fontsize)
			    tick.label2.set_fontweight(fontweight)
	plt.subplots_adjust(wspace=0,hspace=0)

	ax[0,0].set_xlim(-1)
	ax[0,0].set_ylim(-0.05,1.05)
	ax[0,1].set_ylim(-0.2,4.2)
	ax[1,0].set_ylim(0.5,2*10**6)
	ax[1,1].set_ylim(0.5,2*10**6)

	ax[0,0].grid(axis='x')
	ax[0,1].grid(axis='x')
	ax[1,0].grid(axis='x')
	ax[1,1].grid(axis='x')

	if save == True:
		plt.savefig('sim_results/synch-evo-fig1.png')




if __name__ == '__main__':
	
	# Load emission components
	emission = rd.Emission(components=None, Tmin=None, Tmax=None, Emin=0.01, Emax=5e5, dE=60)
	emission.load_spec('sim_results/ordlor_spectrum_therm.txt')
	emission.load_spec('sim_results/ordlor_spectrum_synch.txt')

	z = 0.5
	# Tmin, Tmax = 8,15

	# plot_spec(emission, z=z, comp_num=1, plot_comps=False, nuFnu=True, Tmin=None, Tmax=None, Emin=0.1, Emax=5e5, num_bins=60, save=False)


	fig = plt.figure()
	ax = fig.gca()

	# Make plots of emission
	# plot_lor_dist(save=False)
	Tmin, Tmax = None, None 
	
	# plot_spec(ax,emission, z=z, comp_num=None, plot_comps=True, nuFnu=True, Tmin=Tmin, Tmax=Tmax, Emin=0.1, Emax=5e5, num_bins=60, save=False)
	# Tmin, Tmax = 3,5
	# plot_spec(ax,emission, z=z, comp_num=1, plot_comps=False, nuFnu=True, Tmin=Tmin, Tmax=Tmax, Emin=0.1, Emax=5e5, num_bins=60, save=False)
	# Tmin, Tmax = 8,15
	# plot_spec(ax,emission, z=z, comp_num=1, plot_comps=False, nuFnu=True, Tmin=Tmin, Tmax=Tmax, Emin=0.1, Emax=5e5, num_bins=60, save=False)
	
	plot_light_curve(ax,emission, z=z, comp_num=1, plot_comps=False, Tmin=Tmin, Tmax=Tmax, dt=0.01, Emin=50, Emax=350, save=False)


	# plot_evo_therm(spectrum,z=z,save=False)
	# plot_evo_synch(spectrum,z=z,save=False)
	# plot_param_vs_t(spectrum,z=z,param='gammar',save=False)
	plt.show()






