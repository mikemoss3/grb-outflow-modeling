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


def plot_spec(save=False):
	"""
	Plot Spectrum and components
	"""

	# Load spectrum data 
	spec_therm = np.genfromtxt('sim_results/ordlor_spectrum_therm.txt',dtype=[('te',float),('ta',float),('T',float),('L',float)])
	spec_synch = np.genfromtxt('sim_results/ordlor_spectrum_synch.txt',dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float),('gammar',float),('e',float),('delt',float)])

	fig = plt.figure()
	ax = fig.gca()
	rd.plot_spectrum(ax=ax,nuFnu=True,spec_therm=spec_therm,spec_synch=spec_synch)
	# rd.plot_spectrum(ax=ax,nuFnu=True,spec_therm=spec_therm)
	# rd.plot_spectrum(ax=ax,nuFnu=True,spec_synch=spec_synch)
	ax.vlines(x=4*1e3,ymin=1e46,ymax=1e60,linestyle='dashed',color='black',label='GBM Band')
	ax.vlines(x=4*1e7,ymin=1e46,ymax=1e60,linestyle='dashed',color='black')

	# ax.set_ylim(1e48,1e53)
	ax.set_xlim(100,1e8)
	fig.tight_layout()
	fig.legend()
	
	if save == True:
		plt.savefig('sim_results/spectrum.png')

def plot_light_curve(save=False):
	"""
	Plot light curve per time bin 
	"""

	# Load spectrum data 
	spec_therm = np.genfromtxt('sim_results/ordlor_spectrum_therm.txt',dtype=[('te',float),('ta',float),('T',float),('L',float)])
	spec_synch = np.genfromtxt('sim_results/ordlor_spectrum_synch.txt',dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float),('gammar',float),('e',float),('delt',float)])

	fig = plt.figure()
	ax = fig.gca()
	rd.plot_light_curve(ax=ax,em_therm=spec_therm,em_synch=spec_synch)
	# rd.plot_light_curve(ax=ax,em_therm=spec_therm)
	# rd.plot_light_curve(ax=ax,em_synch=spec_synch)

	fig.tight_layout()
	fig.legend()
	
	if save == True:
		plt.savefig('sim_results/lightcurve.png')

def plot_param_vs_t(param,save=False):
	"""
	Plot Spectrum and components
	"""

	# Load spectrum data 
	spec_therm = np.genfromtxt('sim_results/ordlor_spectrum_therm.txt',dtype=[('te',float),('ta',float),('T',float),('L',float)])
	spec_synch = np.genfromtxt('sim_results/ordlor_spectrum_synch.txt',dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float),('gammar',float),('e',float),('delt',float)])

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

def plot_evo_therm(save=False):
	"""
	Plot evolution of thermal emission parameters 
	"""
	spec_therm = np.genfromtxt('sim_results/ordlor_spectrum_therm.txt',dtype=[('te',float),('ta',float),('T',float),('L',float)])

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



def plot_evo_synch(dt=0.01,save=False):
	"""
	Plot evolution of synchrotron emission parameters 
	"""

	spec_synch = np.genfromtxt('sim_results/ordlor_spectrum_synch.txt',dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float),('gammar',float),('e',float),('delt',float)])

	fontsize=14
	fontweight='bold'


	fig0, ax = plt.subplots(2,1,figsize=(7,8))
	
	ax0cp = ax[0].twinx()
	ax[0].scatter(spec_synch['ta'],spec_synch['te'],marker='^')
	ax0cp.scatter(spec_synch['ta'],spec_synch['delt'],marker='.',color='r')
	ax[0].set_ylabel(r'$t_{e}$',fontsize=fontsize,fontweight=fontweight)
	ax0cp.set_ylabel(r'$\Delta t$',fontsize=fontsize,fontweight=fontweight)
	ax0cp.yaxis.set_label_position("right")
	ax0cp.yaxis.tick_right()


	ax1cp = ax[1].twinx()
	ax[1].scatter(spec_synch['ta'],spec_synch['e']/spec_synch['e'][0],marker='^')

	# tstep=0.1
	# ta_arr = np.arange(start=spec_synch['ta'][0],stop=np.max(spec_synch['ta']+spec_synch['delt']),step=tstep)
	# en_diss_arr = np.zeros(shape=len(ta_arr),dtype=[('t',float),('e',float)])
	# en_diss_arr['t']=ta_arr
	# for i in range(len(spec_synch['e'])):
	# 	start_ind = np.argmax(en_diss_arr['t']>spec_synch['ta'][i])-1
	# 	stop_ind = np.argmax(en_diss_arr['t']>spec_synch['ta'][i]+spec_synch['delt'][i])-1
	# 	en_diss_arr[start_ind:stop_ind]['e'] += spec_synch['e'][i]*tstep
	# ax[1].scatter(en_diss_arr['t'],en_diss_arr['e']/en_diss_arr['e'][0],marker='^')


	ax1cp.scatter(spec_synch['ta'],spec_synch['gammar']/100,marker='.',color='r')
	ax[1].set_ylabel('e (en. diss.)',fontsize=fontsize,fontweight=fontweight)
	ax1cp.set_ylabel(r'$\Gamma_{r}/100$',fontsize=fontsize,fontweight=fontweight)
	ax1cp.yaxis.set_label_position("right")
	ax1cp.yaxis.tick_right()
	ax1cp.set_ylim(0)

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

	ax[0].set_xticklabels([])
	ax0cp.set_xticklabels([])

	plt.tight_layout()
	plt.subplots_adjust(hspace=0)

	if save == True:
		plt.savefig('sim_results/synch-evo-fig0.png')

	fig1, ax = plt.subplots(2,2,figsize=(12,8))

	ax[0,0].scatter(spec_synch['ta'],spec_synch['asyn'])
	ax[0,0].set_ylabel(r'$\alpha_{syn}$',fontsize=fontsize,fontweight=fontweight)
	# ax[0,0].set_xlabel(r't$_a$ (sec), Arrival Time',fontsize=fontsize,fontweight=fontweight)
	
	ax[0,1].scatter(spec_synch['ta'],spec_synch['gammae']/1e4)
	ax[0,1].set_yscale('log')
	ax[0,1].set_ylabel(r'$\Gamma_{e}$/1e4',fontsize=fontsize,fontweight=fontweight)
	# ax[0,1].set_xlabel(r't$_a$ (sec), Arrival Time',fontsize=fontsize,fontweight=fontweight)
	ax[0,1].yaxis.set_label_position("right")
	ax[0,1].yaxis.tick_right()

	ax[1,0].scatter(spec_synch['ta'],spec_synch['Beq'],marker='.')
	ax[1,0].set_yscale('log')
	ax[1,0].set_ylabel(r'B$_{eq}$',fontsize=fontsize,fontweight=fontweight)
	ax[1,0].set_xlabel(r't$_a$ (sec), Arrival Time',fontsize=fontsize,fontweight=fontweight)


	ax[1,1].scatter(spec_synch['ta'],spec_synch['Esyn']/1e3)
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

	ax[0,0].set_xticklabels([])
	ax[0,1].set_xticklabels([])

	if save == True:
		plt.savefig('sim_results/synch-evo-fig1.png')




if __name__ == '__main__':
	# plot_lor_dist(save=False)
	# plot_spec(save=False)
	plot_light_curve(save=False)
	# plot_param_vs_t(param='gammar',save=False)
	# plot_evo_therm(save=False)
	# plot_evo_synch(save=False)
	plt.show()






