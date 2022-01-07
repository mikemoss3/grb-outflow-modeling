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
import cosmologicalconstants as cc

import subprocess
from subprocess import STDOUT


def plot_aesthetics(ax,fontsize=14,fontweight='bold'):

	for tick in ax.xaxis.get_major_ticks():
	    tick.label1.set_fontsize(fontsize=fontsize)
	    tick.label1.set_fontweight(fontweight)
	for tick in ax.yaxis.get_major_ticks():
	    tick.label1.set_fontsize(fontsize=fontsize)
	    tick.label1.set_fontweight(fontweight)
	for tick in ax.xaxis.get_major_ticks():
	    tick.label2.set_fontsize(fontsize=fontsize)
	    tick.label2.set_fontweight(fontweight)
	for tick in ax.yaxis.get_major_ticks():
	    tick.label2.set_fontsize(fontsize=fontsize)
	    tick.label2.set_fontweight(fontweight)

def plot_lor_dist(ax=None,save=False,fontsize=14,fontweight='bold'):
	"""
	Plot Lorentz distributions
	"""
	# Load data
	t0=np.genfromtxt('sim_results/t0_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
	# if os.path.exists('sim_results/t3e4_shells.txt'):
	# 	t3e4=np.genfromtxt('sim_results/t3e4_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
	# if os.path.exists('sim_results/t2e5_shells.txt'):
	# 	t2e5=np.genfromtxt('sim_results/t2e5_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
	# if os.path.exists('sim_results/t5e5_shells.txt'):
	# 	t5e5=np.genfromtxt('sim_results/t5e5_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])

	# Make plot instance if it doesn't exist
	if ax is None:
		fig = plt.figure()
		ax = fig.gca()

	sd.plot_lorentz_dist(ax,t0[t0['STATUS']!=0], label='t0')
	# if os.path.exists('sim_results/t3e4_shells.txt'):
	# 	sd.plot_lorentz_dist(ax,t3e4[t3e4['STATUS']!=0], label='t3e4, s#={}'.format(len(t3e4[t3e4['STATUS']==1])))
	# if os.path.exists('sim_results/t2e5_shells.txt'):
	# 	sd.plot_lorentz_dist(ax,t2e5[t2e5['STATUS']!=0], label='t2e5, s#={}'.format(len(t2e5[t2e5['STATUS']==1])))
	# if os.path.exists('sim_results/t5e5_shells.txt'):
	# 	sd.plot_lorentz_dist(ax,t5e5[t5e5['STATUS']!=0], label='t5e5, s#={}'.format(len(t5e5[t5e5['STATUS']==1])))
	
	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)
	ax.legend(fontsize=fontsize)

	if save == True:
		plt.savefig('sim_results/lorentz-dist.png')


def plot_spec(emission,ax=None,z=0,comp_num=None, plot_comps=True, nuFnu=True, Tmin=None, Tmax=None, Emin=0.1, Emax=5e5, num_bins=60,save=False,fontsize=14,fontweight='bold'):
	"""
	Plot Spectrum and components
	"""
	# Make plot instance if it doesn't exist
	if ax is None:
		fig = plt.figure()
		ax = fig.gca()

	emission.plot_spectrum(ax,comp_num=comp_num,plot_comps=plot_comps, nuFnu=nuFnu, Tmin=Tmin, Tmax=Tmax, Emin=Emin, Emax=Emax, num_bins=num_bins,z=z)

	curr_ymin, curr_ymax = ax.get_ylim()
	curr_ymax *= 2
	ax.vlines(x=4,ymin=1e46,ymax=curr_ymax,linestyle='dashed',color='black',label='Fermi/GBM Band')
	ax.vlines(x=4*1e4,ymin=1e46,ymax=curr_ymax,linestyle='dashed',color='black')

	ax.vlines(x=5,ymin=1e46,ymax=curr_ymax,linestyle='dashed',color='red',label='Swift/BAT Band')
	ax.vlines(x=150,ymin=1e46,ymax=curr_ymax,linestyle='dashed',color='red')

	# Plot aesthetics

	ax.set_xscale('log')
	ax.set_yscale('log')

	# For axis labels
	ax.set_xlabel('E (keV)',fontsize=fontsize,fontweight=fontweight)

	ax.set_ylim(1e48,curr_ymax)
	ax.set_xlim(Emin,Emax)

	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)
	ax.legend(fontsize=fontsize)
	
	if save == True:
		plt.savefig('sim_results/spectrum.png')

def plot_spec_cpp(emission_types, emission_file_names,ax=None,z=0, nuFnu=True, Tmin=None, Tmax=None, Emin=0.1, Emax=5e5, num_bins=60,save=False,fontsize=14,fontweight='bold'):
	"""
	Plot Spectrum and components
	"""
	# Make plot instance if it doesn't exist
	if ax is None:
		fig = plt.figure()
		ax = fig.gca()

	# Write the nuFnu boolean as a string in order to be correctly read by c++ code.
	if nuFnu is True:
		nuFnu = 'true'
	else:
		nuFnu = 'false'

	# Check if multiple emission component files were given:
	if len(emission_file_names) > 1:
		# Initialize the total spectrum array
		spec_tot = np.zeros(shape=num_bins,dtype=[("ENERGY",float),("RATE",float)])
		# For each emission component supplied to the method
		for i in range(len(emission_file_names)):
			# Execute cpp code for each file name specified
			subprocess.run(["./c_scripts/spectrum_make {} {} {} {} {} {} {} {} {} false".format(emission_types[i], emission_file_names[i], Tmin, Tmax, Emin, Emax, num_bins, z, nuFnu)],shell=True,stderr=STDOUT)
			# Load the text file that contains the spectrum and append it
			spec_plot_points = np.genfromtxt('./sim_results/spectrum_points.txt',dtype=[("ENERGY",float),("RATE",float)])
			# Plot spectrum component
			ax.scatter(spec_plot_points['ENERGY'],spec_plot_points['RATE'],label=emission_types[i],marker='.')

			# Add spectrum to total
			spec_tot['RATE'] += spec_plot_points['RATE']
		# Grab energy axis for total spectrum
		spec_tot['ENERGY'] = spec_plot_points['ENERGY']
		# Plot total spectrum
		ax.scatter(spec_tot['ENERGY'],spec_tot['RATE'],label='Total',marker='.')

	# If only one file name was given
	else:
		# Execute the Cpp script on the file name
		subprocess.run(["./c_scripts/spectrum_make {} {} {} {} {} {} {} {} {} false".format(emission_types, emission_file_names, Tmin, Tmax, Emin, Emax, num_bins, z, nuFnu)],shell=True,stderr=STDOUT)
		# Load the text file that contains plot points:
		spec_plot_points = np.genfromtxt('./sim_results/spectrum_points.txt',dtype=[("ENERGY",float),("RATE",float)])
		# Plot the spectrum 
		ax.scatter(spec_plot_points['ENERGY'],spec_plot_points['RATE'],marker='.')


	# Grab the current ymin and ymax, this is used to set the lower and upper bounds of the vertical lines which indicate instrument observation energy range
	curr_ymin, curr_ymax = ax.get_ylim()
	curr_ymax *= 2
	# Display Fermi/GBM energy band
	ax.vlines(x=4,ymin=curr_ymin,ymax=curr_ymax,linestyle='dashed',color='black',label='Fermi/GBM Band')
	ax.vlines(x=4*1e4,ymin=curr_ymin,ymax=curr_ymax,linestyle='dashed',color='black')

	# Display Swift/BAT energy band
	ax.vlines(x=5,ymin=curr_ymin,ymax=curr_ymax,linestyle='dashed',color='red',label='Swift/BAT Band')
	ax.vlines(x=150,ymin=curr_ymin,ymax=curr_ymax,linestyle='dashed',color='red')

	# Plot aesthetics
	ax.set_xscale('log')
	ax.set_yscale('log')

	# For axis labels
	ax.set_xlabel('E (keV)',fontsize=fontsize,fontweight=fontweight)

	ax.set_ylim(1e48,curr_ymax)
	ax.set_xlim(Emin,Emax)

	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)
	ax.legend(fontsize=fontsize)
	
	if save == True:
		plt.savefig('sim_results/spectrum.png')

def plot_light_curve(emission,ax=None,z=0,comp_num=None, plot_comps=True, Tmin=None, Tmax=None,dt=0.05, Emin=0.1, Emax=5e5,save=False,fontsize=14,fontweight='bold'):
	"""
	Plot light curve per time bin 
	"""

	# Make plot instance if it doesn't exist
	if ax is None:
		fig = plt.figure()
		ax = fig.gca()

	emission.plot_light_curve(ax, comp_num=comp_num,plot_comps=plot_comps, Tmin=Tmin, Tmax=Tmax, dt=dt, Emin=Emin, Emax=Emax, z=z)

	# Plot aesthetics
	ax.set_yscale('log')
	ax.set_xlim(Tmin,Tmax)

	# For axis labels
	ax.set_ylabel('Rate (cts/sec)',fontsize=fontsize,fontweight=fontweight)
	ax.set_xlabel('Time (sec)',fontsize=fontsize,fontweight=fontweight)

	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)
	ax.legend(fontsize=fontsize)
	
	if save == True:
		plt.savefig('sim_results/lightcurve.png')

def plot_light_curve_cpp(emission_types,emission_file_names, ax=None, z=0, Tmin=None, Tmax=None, dt=0.05, Emin=0.1, Emax=5e5,save=False,fontsize=14,fontweight='bold'):
	"""
	Plot light curve per time bin 
	"""

	# Make plot instance if it doesn't exist
	if ax is None:
		fig = plt.figure()
		ax = fig.gca()

	# Number of time bins to use
	num_bins = int((Tmax-Tmin)/dt)

	# Check if multiple file names were given:
	if len(emission_file_names) > 1:
		# Initialize the total spectrum array
		lc_tot = np.zeros(shape=num_bins,dtype=[("TIME",float),("RATE",float)])
		# For each emission component supplied to the method
		for i in range(len(emission_file_names)):
			# Execute cpp code for each file name specified
			subprocess.run(["./c_scripts/light_curve_make {} {} {} {} {} {} {} {}".format(emission_types[i], emission_file_names[i], Tmin, Tmax, dt, Emin, Emax, z)],shell=True,stderr=STDOUT)
			# Load the text file that contains the spectrum and append it
			lc_plot_points = np.genfromtxt('./sim_results/light_curve_points.txt',dtype=[("TIME",float),("RATE",float)])
			# Plot spectrum component
			ax.scatter(lc_plot_points['TIME'],lc_plot_points['RATE'],label=emission_types[i],marker='.')

			# Add spectrum to total
			lc_tot['RATE'] += lc_plot_points['RATE']
		# Grab energy axis for total spectrum
		lc_tot['TIME'] = lc_plot_points['TIME']
		# Plot total spectrum
		ax.scatter(lc_tot['TIME'],lc_tot['RATE'],label='Total',marker='.')
		ax.legend(fontsize=fontsize)

	# If only one file name was given
	else:
		# Execute the Cpp script on the file name
		subprocess.run(["./c_scripts/light_curve_make {} {} {} {} {} {} {} {}".format(emission_types[0], emission_file_names[0], Tmin, Tmax, dt, Emin, Emax, z)],shell=True,stderr=STDOUT)
		# Load the text file that contains plot points:
		lc_plot_points = np.genfromtxt('./sim_results/light_curve_points.txt',dtype=[("TIME",float),("RATE",float)])
		# Plot the spectrum 
		ax.scatter(lc_plot_points['TIME'],lc_plot_points['RATE'],marker='.')


	# Plot aesthetics
	ax.set_xlim(Tmin,Tmax)

	# For axis labels
	ax.set_ylabel(r'Rate (ph cm$^{-2}$ s$^{-1}$)',fontsize=fontsize,fontweight=fontweight)
	ax.set_xlabel('Obs Time (sec)',fontsize=fontsize,fontweight=fontweight)

	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)
	
	if save == True:
		plt.savefig('sim_results/lightcurve.png')

def plot_param_vs_ta(emission_comp,param,ax=None,z=0, y_factor=1, Tmin=None, Tmax=None,save=False,fontsize=14,fontweight='bold',disp_xax=True,disp_yax=True,color='mediumblue',marker='.'):
	"""
	Plot emission parameters as a function of time (in the observer frame)
	"""

	# Make plot instance if it doesn't exist
	if ax is None:
		fig = plt.figure()
		ax = fig.gca()

	ax_time = emission_comp['ta'] * (1+z)

	# Find the indices of the start and stop time
	ind_start, ind_stop = 0,-1 
	if Tmin is not None:
		ind_start = np.argmax(ax_time>Tmin)
	if Tmax is not None:
		ind_stop = np.argmax(ax_time>Tmax)

	# Load time axis and the parameter value axis 
	ax_time = emission_comp['ta'][ind_start:ind_stop] 
	ax_param = emission_comp[param][ind_start:ind_stop] * y_factor

	ax.scatter(ax_time,ax_param,label=param,c=color,marker=marker)

	if disp_yax is True:
		ax.set_ylabel(param,fontsize=fontsize,fontweight=fontweight)
	if disp_xax is True:
		ax.set_xlabel(r't$_a$',fontsize=fontsize,fontweight=fontweight)

	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)
		
	if save == True:
		plt.savefig('sim_results/param-{}-vs-t.png'.format(param))

def plot_evo_therm(thermal_emission,ax=None,z=0,Tmin=None, Tmax=None,save=False,fontsize=14,fontweight='bold'):
	"""
	Plot evolution of thermal emission parameters 
	"""

	# Make plot instance if it doesn't exist
	if ax is None:
		fig, ax = plt.subplots(2,1,figsize=(5,8))

	# Plot temperature of the thermal component vs time (in observer frame)
	plot_param_vs_ta(thermal_emission,'T', ax=ax[0], z=z, y_factor=cc.kb_kev/(1+z),Tmin=Tmin, Tmax=Tmax, save=save,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False)
	
	ax[0].set_xlabel(r't$_{obs}$',fontsize=fontsize,fontweight=fontweight)
	ax[0].set_ylabel(r'k$_B$T (KeV)',fontsize=fontsize,fontweight=fontweight)
	ax[0].set_yscale('log')

	# Plot Rphot vs Tphot
	ax[1].scatter(thermal_emission['Rphot'],thermal_emission['T'])
	
	ax[1].set_xlabel(r'R$_{phot}$ (light sec)',fontsize=fontsize,fontweight=fontweight)
	ax[1].set_ylabel(r'T$_{phot}}$ (K)',fontsize=fontsize,fontweight=fontweight)

	
	plot_aesthetics(ax[0],fontsize=fontsize,fontweight=fontweight)
	plot_aesthetics(ax[1],fontsize=fontsize,fontweight=fontweight)
	
	if save == True:
		plt.savefig('sim_results/thermal-evo.png')

def plot_evo_synch(synch_emission,ax=None,z=0,Tmin=None, Tmax=None,save=False,fontsize=14,fontweight='bold'):
	"""
	Plot evolution of synchrotron emission parameters 
	"""

	if ax is None:
		fig, ax = plt.subplots(2,1,figsize=(7,8),sharex=True)
	
	# Make a copy of the axis in order to over plot two separate data sets
	ax0cp = ax[0].twinx()

	# Plot Arrival Time (ta) vs Emission Time (te)
	plot_param_vs_ta(synch_emission,'te', ax=ax[0], z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker='^')
	# Plot Arrival Time (ta) vs delta T
	plot_param_vs_ta(synch_emission,'delt', ax=ax0cp, z=z,Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker='.',color='r')
	
	ax[0].set_ylabel(r'$t_{e}$',fontsize=fontsize,fontweight=fontweight)
	ax0cp.set_ylabel(r'$\Delta t$',fontsize=fontsize,fontweight=fontweight)
	ax0cp.yaxis.set_label_position("right")
	ax0cp.yaxis.tick_right()

	# Make a copy of the axis in order to over plot two separate data sets
	ax1cp = ax[1].twinx()
	
	ax[1].scatter(synch_emission['ta']*(1+z),synch_emission['e']/synch_emission['e'][0],marker='^')
	ax1cp.scatter(synch_emission['ta']*(1+z),synch_emission['gammar']/100,marker='.',color='r')

	# Plot Arrival Time (ta) vs the dissipated energy (e)
	plot_param_vs_ta(synch_emission,'e', ax=ax[1], y_factor=1/synch_emission['e'][0], z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker='^')
	# Plot Arrival Time (ta) vs approximate Lorentz factor (gamma_r)
	plot_param_vs_ta(synch_emission,'gammar', ax=ax1cp, y_factor=1/100, z=z,Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker='.',color='r')

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
		plot_aesthetics(ax[i],fontsize=fontsize,fontweight=fontweight)
	for twin in [ax0cp,ax1cp]:
		plot_aesthetics(twin,fontsize=fontsize,fontweight=fontweight)

	plt.tight_layout()
	plt.subplots_adjust(hspace=0)

	if save == True:
		plt.savefig('sim_results/synch-evo-fig0.png')

	fig, ax = plt.subplots(2,2,sharex=True,figsize=(12,8))

	# Plot Arrival Time (ta) vs the energy fraction in synchrotron electron (asyn)
	plot_param_vs_ta(synch_emission,'asyn', ax=ax[0,0], z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False)
	ax[0,0].set_ylabel(r'$\alpha_{syn}$',fontsize=fontsize,fontweight=fontweight)
	
	# Plot Arrival Time (ta) vs the dissipated energy (e)
	plot_param_vs_ta(synch_emission,'gammae', ax=ax[0,1], y_factor=1/1e4,z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False)
	ax[0,1].set_ylabel(r'$\Gamma_{e}$/1e4',fontsize=fontsize,fontweight=fontweight)
	ax[0,1].yaxis.set_label_position("right")
	ax[0,1].yaxis.tick_right()

	# Plot Arrival Time (ta) vs the equipartition magnetic field (Beq)
	plot_param_vs_ta(synch_emission,'Beq', ax=ax[1,0], z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker='.')
	ax[1,0].set_yscale('log')
	ax[1,0].set_ylabel(r'B$_{eq}$',fontsize=fontsize,fontweight=fontweight)
	ax[1,0].set_xlabel(r't$_a$ (sec), Arrival Time',fontsize=fontsize,fontweight=fontweight)

	# Plot Arrival Time (ta) vs the synchrotron energy (Esyn)
	plot_param_vs_ta(synch_emission,'Esyn', ax=ax[1,1], z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False)
	ax[1,1].set_yscale('log')
	ax[1,1].set_ylabel(r'$E_{syn}$/1e3',fontsize=fontsize,fontweight=fontweight)
	ax[1,1].set_xlabel(r't$_a$ (sec), Arrival Time',fontsize=fontsize,fontweight=fontweight)
	ax[1,1].yaxis.set_label_position("right")
	ax[1,1].yaxis.tick_right()

	for i in range(2):
		for j in range(2):
			plot_aesthetics(ax[i,j],fontsize=fontsize,fontweight=fontweight)

	plt.tight_layout()
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

	z = 1
	Tmin, Tmax = None, None 

	# Make plots of emission:
	fig = plt.figure()
	ax = fig.gca()
	
	## For one comp:
	# plot_spec(emission,ax=ax, z=z, comp_num=0, plot_comps=False, nuFnu=True, Tmin=Tmin, Tmax=Tmax,Emin=0.1, Emax=1e4, num_bins=50, save=False)
	# plot_light_curve(emission,ax=ax, z=z, comp_num=1, plot_comps=False, Tmin=Tmin, Tmax=Tmax, dt=0.1, Emin=50, Emax=300, save=False)

	## For both Comps:
	# plot_spec(emission,ax=ax, z=z, plot_comps=True, nuFnu=True, Tmin=Tmin, Tmax=Tmax, num_bins=140, save=False)
	# plot_light_curve(emission,ax=ax, z=z, plot_comps=True, Tmin=Tmin, Tmax=Tmax, dt=1, Emin=50, Emax=300, save=False)

	## Other emission plots:
	plot_lor_dist(save=False,ax=ax)
	# plot_evo_therm(emission.components[0],z=z,save=False)
	# plot_evo_synch(emission.components[1],z=z,save=False)
	# plot_param_vs_ta(emission.components[1],ax=ax,z=z,param='Esyn',save=False)

	
	# plot_spec_cpp(["thermal","synchrotron"],['sim_results/ordlor_spectrum_therm.txt','sim_results/ordlor_spectrum_synch.txt'], ax=ax, z=0.5, nuFnu=True, Tmin=0, Tmax=20, Emin=0.1, Emax=1e5, num_bins=160, save=False)
	# plot_light_curve_cpp(["thermal","synchrotron"],['sim_results/ordlor_spectrum_therm.txt','sim_results/ordlor_spectrum_synch.txt'], ax=ax, z=z, Tmin=0, Tmax=30, dt=0.05, Emin=8, Emax=40000, save=False)
	# plot_light_curve_cpp(["thermal"],['sim_results/ordlor_spectrum_therm.txt'], ax=ax, z=z, Tmin=0, Tmax=20, dt=0.05, Emin=0.1, Emax=5e5, save=False)

	fig.tight_layout()
	plt.show()

	
