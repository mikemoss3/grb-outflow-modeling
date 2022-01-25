"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-11-26


Meta script to plot desired simulation results created by c++ code.

"""
import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
# import cosmologicalconstants as cc
import scipy.integrate as integrate 
# from matplotlib.widgets import RangeSlider
import time
from matplotlib.widgets import TextBox

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

##############################################################################################################################

def plot_lor_dist(file_name,ax=None,save_pref=None,xlabel=True,ylabel=True,label=None,fontsize=14,fontweight='bold',linestyle='solid'):
	"""
	Method to plot the given Lorentz factor distribution

	Attributes:
	ax = the matplotlib.pyplot.axes instance to make the plot on
	save_pref = if not left as None, the plot will be saved and the file name will have this prefix
	xlabel, ylabel = indicate whether x- and y- labels should be included (boolean)
	fontsize, fontweight = fontsize and fontweight of the plot font and labels on the plot
	linestyle = style of the plotting line 
	"""

	# Load data
	shell_dist = np.genfromtxt(file_name,dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])

	# Make plot instance if it doesn't exist
	if ax is None:
		ax = plt.figure().gca()

	# To match Daigne and Mochkovitch 1998 paper figures
	flipped_mass_arr = np.flip(shell_dist['MASS'])
	flipped_gamma_arr = np.flip(shell_dist['GAMMA'])

	# Cumulative mass
	masscum = np.cumsum(flipped_mass_arr)
	massfraccum = masscum/masscum[-1]

	# Plot distribution
	line, = ax.step(massfraccum,flipped_gamma_arr,where='pre',linestyle=linestyle,label=label)

	if xlabel is True:
		ax.set_xlabel(r'M/M$_{tot}$',fontsize=fontsize,fontweight=fontweight)
	if ylabel is True:
		ax.set_ylabel(r'$\Gamma$',fontsize=fontsize,fontweight=fontweight)

	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)
	if label is not None:
		ax.legend(fontsize=fontsize)

	if save_pref is not None :
		plt.savefig('figs/{}-lorentz-dist.png'.format(save_pref))

##############################################################################################################################

def plot_spec(file_name, z=0, joined=False, label = None, ax=None, nuFnu=True, unc=False, Emin=None, Emax=None, save_pref=None,fontsize=14,fontweight='bold'):
	"""
	Method to plot the input spectrum data files

	Attributes:
	file_name = file name which contains spectrum data points 
	joined = boolean, indicates whether the points are joined or not.
	label = optional label for the plotted spectra 
	ax = the matplotlib.pyplot.axes instance to make the plot on
	
	nuFnu = boolean, indicates whether the spectrum should be a count spectrum or energy density spectrum
	unc = boolean, indicates whether to include uncertainty bars on the data points 
	Emin, Emax = indicates the minimum and maximum energy range to plot. If None is supplied, the minimum and maximum energies of the supplied data files are used

	save_pref = if not left as None, the plot will be saved and the file name will have this prefix
	xlabel, ylabel = indicate whether x- and y- labels should be included (boolean)
	fontsize, fontweight = fontsize and fontweight of the plot font and labels on the plot
	linestyle = style of the plotting line 
	"""

	# Make plot instance if it doesn't exist
	if ax is None:
		ax = plt.figure().gca()

	# Load spectrum data
	spec_data = np.genfromtxt(file_name,dtype=[("ENERG",float),("RATE",float),('UNC',float)])

	spec_data['ENERG'] /= (1+z)
	spec_data['RATE'] /= (1+z)
	spec_data['UNC'] /= (1+z)

	if joined is True:
		# Plot spectrum data
		if nuFnu is True:
			if unc is True:
				line = ax.errorbar(x=spec_data['ENERG'],y=spec_data['RATE']*(spec_data['ENERG']**2),yerr=spec_data['UNC']*(spec_data['ENERG']**2),label=label)
			else:
				line, = ax.plot(spec_data['ENERG'],spec_data['RATE']*(spec_data['ENERG']**2),label=label)
		else:
			if unc is True:
				line = ax.errorbar(x=spec_data['ENERG'],y=spec_data['RATE'],yerr=spec_data['UNC'],label=label)
			else:
				line = ax.plot(spec_data['ENERG'],spec_data['RATE'],label=label)
	else:
		# Plot spectrum data
		if nuFnu is True:
			if unc is True:
				line = ax.errorbar(x=spec_data['ENERG'],y=spec_data['RATE']*(spec_data['ENERG']**2),yerr=spec_data['UNC']*(spec_data['ENERG']**2),label=label,fmt=" ",marker="+")
			else:
				line = ax.errorbar(x=spec_data['ENERG'],y=spec_data['RATE']*(spec_data['ENERG']**2),label=label,fmt=" ",marker="+")
		else:
			if unc is True:
				line = ax.errorbar(x=spec_data['ENERG'],y=spec_data['RATE'],yerr=spec_data['UNC'],label=label,fmt=" ",marker="+")
			else:
				line = ax.errorbar(x=spec_data['ENERG'],y=spec_data['RATE'],label=label,fmt=" ",marker="+")

	# Plot aesthetics
	ax.set_xscale('log')
	ax.set_yscale('log')

	# Force lower bound
	ax.set_ylim(1e48,1e53)

	# For axis labels
	ax.set_xlabel('E (keV)',fontsize=fontsize,fontweight=fontweight)

	if nuFnu is True:
		ax.set_ylabel(r'$\nu$F$_\nu$ erg sec$^{-1}$ cm$^{-2}$ keV$^{-2}$',fontsize=fontsize,fontweight=fontweight)
	else:
		ax.set_ylabel( "N(E) counts sec$^{-1}$ keV$^{-1}$",fontsize=fontsize,fontweight=fontweight)

	# curr_ymin, curr_ymax = ax.get_ylim()
	# ax.set_ylim(curr_ymin,curr_ymax)
	# ax.set_xlim(Emin,Emax)

	# Add label names to plot if supplied
	if label is not None:
		ax.legend(fontsize=fontsize-2)

	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)
	
	# plt.tight_layout()
	if save_pref is not None:
		plt.savefig('figs/{}-spectrum.png'.format(save_pref))	

	return line

##############################################################################################################################

def add_FermiGBM_band(ax,fontsize=12):
	# Grab the current ymin and ymax, this is used to set the lower and upper bounds of the vertical lines which indicate instrument observation energy range
	curr_ymin, curr_ymax = ax.get_ylim()
	curr_xmin, curr_xmax = ax.get_xlim()

	# Display Fermi/GBM - NAI energy band
	ax.vlines(x=8,ymin=curr_ymin,ymax=curr_ymax,linestyle='dashed',color='black',label='Fermi/GBM-NAI Band')
	ax.vlines(x=1e3,ymin=curr_ymin,ymax=curr_ymax,linestyle='dashed',color='black')

	# Display Fermi/GBM - BGO energy band
	ax.vlines(x=150,ymin=curr_ymin,ymax=curr_ymax,linestyle='dashed',color='orange',label='Fermi/GBM-BGO Band')
	ax.vlines(x=3*1e4,ymin=curr_ymin,ymax=curr_ymax,linestyle='dashed',color='orange')

	# Add to legend	
	ax.legend(fontsize=fontsize)

	# We don't want the plotting window to change if either of the energy band edges do not overlap with the plotted energy spectra
	ax.set_ylim(curr_ymin,curr_ymax)
	ax.set_xlim(curr_xmin,curr_xmax)

##############################################################################################################################

def add_SwiftBAT_band(ax,fontsize=12):
	# Grab the current ymin and ymax, this is used to set the lower and upper bounds of the vertical lines which indicate instrument observation energy range
	curr_ymin, curr_ymax = ax.get_ylim()
	curr_xmin, curr_xmax = ax.get_xlim()

	# Display Swift/BAT energy band
	ax.vlines(x=5,ymin=curr_ymin,ymax=curr_ymax,linestyle='dashed',color='red',label='Swift/BAT Band')
	ax.vlines(x=350,ymin=curr_ymin,ymax=curr_ymax,linestyle='dashed',color='red')

	# Add to legend	
	ax.legend(fontsize=fontsize)

	# We don't want the plotting window to change if either of the energy band edges do not overlap with the plotted energy spectra
	ax.set_ylim(curr_ymin,curr_ymax)
	ax.set_xlim(curr_xmin,curr_xmax)

##############################################################################################################################

def plot_light_curve(file_name, z=0, label=None, ax=None, Tmin=None, Tmax=None, save_pref=None, fontsize=14,fontweight='bold'):
	"""
	Method to plot the input light curve data files

	Attributes:
	file_name = file name which contains spectrum data points 
	label = optional label for the plotted light curve 
	ax = the matplotlib.pyplot.axes instance to make the plot on
	
	Tmin, Tmax = indicates the minimum and maximum time range to plot. If None is supplied, the minimum and maximum times of the supplied data files are used

	save_pref = if not left as None, the plot will be saved and the file name will have this prefix
	xlabel, ylabel = indicate whether x- and y- labels should be included (boolean)
	fontsize, fontweight = fontsize and fontweight of the plot font and labels on the plot
	linestyle = style of the plotting line 
	"""

	if(z<0):
		print("Please provide a non-negative redshift.")
		return;
	else:
		# Make plot instance if it doesn't exist
		if ax is None:
			ax = plt.figure().gca()

		# Load light curve data
		light_curve_data = np.genfromtxt(file_name,dtype=[("TIME",float),("RATE",float)])

		# Plot light curve data

		if(z>0):
			# ax.scatter(light_curve_data['TIME']*(1+z),light_curve_data['RATE']/(4*np.pi*lum_dis(z)**2),label=label,marker=".")
			ax.step(light_curve_data['TIME']*(1+z),light_curve_data['RATE']/(4*np.pi*lum_dis(z)**2),label=label,marker=" ",where="mid")
		else: 
			# If z = 0, return luminosity
			# ax.scatter(light_curve_data['TIME'],light_curve_data['RATE'],label=label,marker=".")
			ax.step(light_curve_data['TIME'],light_curve_data['RATE'],label=label,marker=" ",where="mid")

		# Plot aesthetics
		# For axis labels
		ax.set_ylabel(r'Rate (ph cm$^{-2}$ s$^{-1}$)',fontsize=fontsize,fontweight=fontweight)
		ax.set_xlabel('Obs Time (sec)',fontsize=fontsize,fontweight=fontweight)

		# Add label names to plot if supplied
		if label is not None:
			plt.legend(fontsize=fontsize-2)

		plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)
		
		plt.tight_layout()
		if save_pref is not None:
			plt.savefig('figs/{}-light-curve.png'.format(save_pref))

##############################################################################################################################

def plot_light_curve_interactive(file_name, z=0, label=None, ax=None, Tmin=None, Tmax=None, save_pref=None, fontsize=14,fontweight='bold'):
	"""
	Method to plot the input light curve data files

	Attributes:
	file_name = file name which contains spectrum data points 
	label = optional label for the plotted light curve 
	ax = the matplotlib.pyplot.axes instance to make the plot on
	
	Tmin, Tmax = indicates the minimum and maximum time range to plot. If None is supplied, the minimum and maximum times of the supplied data files are used

	save_pref = if not left as None, the plot will be saved and the file name will have this prefix
	xlabel, ylabel = indicate whether x- and y- labels should be included (boolean)
	fontsize, fontweight = fontsize and fontweight of the plot font and labels on the plot
	linestyle = style of the plotting line 
	"""

	if(z<0):
		print("Please provide a non-negative redshift.")
		return;
	else:
		# Make plot instance if it doesn't exist
		fig, ax = plt.subplots(1,2,figsize=(12,6))	
		
		# Load light curve data
		light_curve_data = np.genfromtxt(file_name,dtype=[("TIME",float),("RATE",float)])

		# Plot light curve data

		if(z>0):
			# ax.scatter(light_curve_data['TIME']*(1+z),light_curve_data['RATE']/(4*np.pi*lum_dis(z)**2),label=label,marker=".")
			ax[0].step(light_curve_data['TIME']*(1+z),light_curve_data['RATE']/(4*np.pi*lum_dis(z)**2),label=label,marker=" ",where="mid")
		else: 
			# If z = 0, return luminosity
			# ax.scatter(light_curve_data['TIME'],light_curve_data['RATE'],label=label,marker=".")
			ax[0].step(light_curve_data['TIME'],light_curve_data['RATE'],label=label,marker=" ",where="mid")

		tmin = np.max(light_curve_data['TIME']*(1+z))/4
		tmax = np.max(light_curve_data['TIME']*(1+z))*3/4
		subprocess.run(["./main","quickplot", "{}".format(tmin), "{}".format(tmax)])

		line = plot_spec("data-file-dir/quickplot_spectrum.txt",ax=ax[1],z=z,joined=True)
		add_FermiGBM_band(ax[1])

		fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)
		axbox = plt.axes([0.2, 0.08, 0.22, 0.04])
		text_box = TextBox(ax=axbox, label='Time Min.', initial="{0:.2f}, {1:.2f}".format(tmin,tmax) )

		# Create the Vertical lines on the plot 
		lower_limit_line = ax[0].axvline(np.max(light_curve_data['TIME']*(1+z))/4, color='k')
		upper_limit_line = ax[0].axvline(np.max(light_curve_data['TIME']*(1+z))*3/4, color='k')

		# Plot aesthetics
		# For axis labels
		ax[0].set_ylabel(r'Rate (ph cm$^{-2}$ s$^{-1}$)',fontsize=fontsize,fontweight=fontweight)
		ax[0].set_xlabel('Obs Time (sec)',fontsize=fontsize,fontweight=fontweight)

		# Add label names to plot if supplied
		if label is not None:
			ax[0].legend(fontsize=fontsize-2)

		plot_aesthetics(ax[0],fontsize=fontsize,fontweight=fontweight)
		fig.canvas.draw()
		
		def submit(val):
			# val will be separated into two values, split by the comma
			x = val.split(", ")

			tmin = float(x[0])
			tmax = float(x[1])

			# Update the position of the vertical lines
			lower_limit_line.set_xdata(tmin)
			upper_limit_line.set_xdata(tmax)

			subprocess.run(["./main","quickplot", "{}".format(tmin), "{}".format(tmax)])

			spec_data = np.genfromtxt("data-file-dir/quickplot_spectrum.txt",dtype=[("ENERG",float),("RATE",float),("UNC",float)])
			spec_data["ENERG"] /= (1+z)
			spec_data["RATE"] /= (1+z)
			line.set_xdata(spec_data["ENERG"])
			line.set_ydata( spec_data["RATE"]*(spec_data['ENERG']**2) )
			ax[1].redraw_in_frame()

			# Redraw the figure to ensure it updates
			fig.canvas.draw_idle()
			
		text_box.on_submit(submit)

		fig.subplots_adjust(bottom=0.25,wspace=0.4)

		# plt.tight_layout()
		if save_pref is not None:
			plt.savefig('figs/{}-light-curve.png'.format(save_pref))

		return text_box

##############################################################################################################################

def load_therm_emission(file_name):
	"""
	Method to load thermal emission data from the given file name
	"""

	dtype = np.dtype([('TE',float),('TA',float),('DELT',float),('T',float),('FLUX',float),('RPHOT',float)])

	return np.genfromtxt(file_name,dtype=dtype)

##############################################################################################################################

def load_is_emission(file_name):
	"""
	Method to load internal shock emission data from the given file name
	"""

	dtype = np.dtype([('TE',float),('TA',float),('DELT',float),('BEQ',float),('GAMMAE',float),('ESYN',float),('GAMMAR',float),('EDISS',float),('ASYN',float),('TAU',float),('RELVEL',float)])

	return np.genfromtxt(file_name,dtype=dtype)

##############################################################################################################################

def load_fs_emission(file_name):
	"""
	Method to load forward shock emission data from the given file name
	"""

	dtype = np.dtype([('TE',float),('TA',float),('DELT',float),('BEQ',float),('GAMMAE',float),('ESYN',float),('GAMMAR',float),('EDISS',float)])

	return np.genfromtxt(file_name,dtype=dtype)

##############################################################################################################################

def load_rs_emission(file_name):
	"""
	Method to load reverse shock emission data from the given file name
	"""

	dtype = np.dtype([('TE',float),('TA',float),('DELT',float),('BEQ',float),('GAMMAE',float),('ESYN',float),('GAMMAR',float),('EDISS',float)])

	return np.genfromtxt(file_name,dtype=dtype)

##############################################################################################################################

def plot_param_vs_ta(emission_comp,param,ax=None,z=0, y_factor=1, label=None, Tmin=None, Tmax=None,save_pref=None,fontsize=14,fontweight='bold',disp_xax=True,disp_yax=True,color='C0',marker='.'):
	"""
	Plot emission parameters as a function of time (in the observer frame)
	"""

	# Make plot instance if it doesn't exist
	if ax is None:
		ax = plt.figure().gca()

	# Multiply by 1+z for the time axis and apply the supplied factor on the y-axis 
	ax_time = emission_comp['TA'] * (1+z)
	emission_comp[param]*=y_factor

	# Find the indices of the start and stop time
	ind_start, ind_stop = 0,-1 
	if Tmin is not None:
		ind_start = np.argmax(ax_time>Tmin)
	if Tmax is not None:
		ind_stop = np.argmax(ax_time>Tmax)

	# Load time axis and the parameter value axis 
	ax_time = ax_time[ind_start:ind_stop] 
	ax_param = emission_comp[param][ind_start:ind_stop]

	if label is None:
		label = param

	ax.scatter(ax_time,ax_param,label=label,c=color,marker=marker)

	if disp_yax is True:
		ax.set_ylabel(param,fontsize=fontsize,fontweight=fontweight)
	if disp_xax is True:
		ax.set_xlabel(r't$_a$',fontsize=fontsize,fontweight=fontweight)

	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)
		
	if save_pref is not None:
		plt.savefig('figs/{}-param-{}-vs-t.png'.format(save_pref,param))

##############################################################################################################################

def plot_evo_therm(thermal_emission,ax=None,z=0,Tmin=None, Tmax=None,save_pref=None,fontsize=14,fontweight='bold'):
	"""
	Plot evolution of thermal emission parameters 
	"""

	# Make plot instance if it doesn't exist
	if ax is None:
		fig, ax = plt.subplots(2,1,figsize=(5,8))

	# Plot temperature of the thermal component vs time (in observer frame)
	kb_kev = 8.617*1e-8
	plot_param_vs_ta(thermal_emission,'T', ax=ax[0], z=z, y_factor=kb_kev/(1+z),Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False)
	
	ax[0].set_xlabel(r't$_{obs}$',fontsize=fontsize,fontweight=fontweight)
	ax[0].set_ylabel(r'k$_B$T (KeV)',fontsize=fontsize,fontweight=fontweight)
	ax[0].set_yscale('log')

	# Plot Rphot vs Tphot
	ax[1].scatter(thermal_emission['RPHOT'],thermal_emission['T'])
	
	ax[1].set_xlabel(r'R$_{phot}$ (light sec)',fontsize=fontsize,fontweight=fontweight)
	ax[1].set_ylabel(r'T$_{phot}}$ (K)',fontsize=fontsize,fontweight=fontweight)

	plot_aesthetics(ax[0],fontsize=fontsize,fontweight=fontweight)
	plot_aesthetics(ax[1],fontsize=fontsize,fontweight=fontweight)
	
	plt.tight_layout()
	if save_pref is not None:
		plt.savefig('figs/{}-thermal-evo.png'.format(save_pref))

##############################################################################################################################

def plot_evo_int_shock(is_emission,ax=None,z=0,Tmin=None, Tmax=None,save_pref=None,fontsize=14,fontweight='bold'):
	"""
	Plot evolution of internal shock emission parameters 
	"""

	if ax is None:
		fig, ax = plt.subplots(2,1,figsize=(7,8),sharex=True)
	
	# Make a copy of the axis in order to over plot two separate data sets
	ax0cp = ax[0].twinx()

	# Plot Arrival Time (ta) vs Emission Time (te)
	plot_param_vs_ta(is_emission,'TE', ax=ax[0], z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker='^')
	# Plot Arrival Time (ta) vs delta T
	plot_param_vs_ta(is_emission,'DELT', ax=ax0cp, z=z,Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker='.',color='C1')
	
	ax[0].set_ylabel(r'$t_{e}$',fontsize=fontsize,fontweight=fontweight)
	ax0cp.set_ylabel(r'$\Delta t$',fontsize=fontsize,fontweight=fontweight)
	ax0cp.yaxis.set_label_position("right")
	ax0cp.yaxis.tick_right()

	# Make a copy of the axis in order to over plot two separate data sets
	ax1cp = ax[1].twinx()
	
	# Plot Arrival Time (ta) vs the dissipated energy (e)
	plot_param_vs_ta(is_emission,'EDISS', ax=ax[1], y_factor=is_emission['DELT'], z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker='^')
	# Plot Arrival Time (ta) vs approximate Lorentz factor (gamma_r)
	plot_param_vs_ta(is_emission,'GAMMAR', ax=ax1cp, y_factor=1/100, z=z,Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker='.',color='C1')

	ax[1].set_ylabel('e (en. diss.)',fontsize=fontsize,fontweight=fontweight)
	ax[1].set_xlabel(r't$_a$ (sec), Arrival Time',fontsize=fontsize,fontweight=fontweight)
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

	if save_pref is not None:
		plt.savefig('figs/{}-int-shock-evo-fig0.png',format(save_pref))

	fig, ax = plt.subplots(2,2,sharex=True,figsize=(12,8))

	# Plot Arrival Time (ta) vs the energy fraction in synchrotron electron (asyn)
	plot_param_vs_ta(is_emission,'ASYN', ax=ax[0,0], z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False)
	ax[0,0].set_ylabel(r'$\alpha_{syn}$',fontsize=fontsize,fontweight=fontweight)
	
	# Plot Arrival Time (ta) vs the dissipated energy (e)
	plot_param_vs_ta(is_emission,'GAMMAE', ax=ax[0,1], y_factor=1/1e4,z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False)
	ax[0,1].set_ylabel(r'$\Gamma_{e}$/1e4',fontsize=fontsize,fontweight=fontweight)
	ax[0,1].yaxis.set_label_position("right")
	ax[0,1].yaxis.tick_right()

	# Plot Arrival Time (ta) vs the equipartition magnetic field (Beq)
	plot_param_vs_ta(is_emission,'BEQ', ax=ax[1,0], z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker='.')
	ax[1,0].set_yscale('log')
	ax[1,0].set_ylabel(r'B$_{eq}$',fontsize=fontsize,fontweight=fontweight)
	ax[1,0].set_xlabel(r't$_a$ (sec), Arrival Time',fontsize=fontsize,fontweight=fontweight)

	# Plot Arrival Time (ta) vs the synchrotron energy (Esyn)
	plot_param_vs_ta(is_emission,'ESYN', ax=ax[1,1], z=z, Tmin=Tmin, Tmax=Tmax,
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

	if save_pref == True:
		plt.savefig('figs/{}-int-shock-evo-fig1.png'.format(save_pref))

##############################################################################################################################

def plot_evo_ext_shock(fs_data=None, rs_data=None, z=0, Tmin=None, Tmax=None,save_pref=None,fontsize=14,fontweight='bold'):

	fig, ax = plt.subplots(2,2,sharex=True,figsize=(12,8))

	if( (fs_data is None) and (rs_data is None)):
		print("Please provide either forward shock or reverse shock data.")
		return
	if(fs_data is not None):
		# T_a vs T_e
		plot_param_vs_ta(fs_data,'TE', ax=ax[0,0], z=z, Tmin=Tmin, Tmax=Tmax,
				fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker=',',color="C1")

		# T_a vs del T_a
		
		# T_a vs Gamma_r
		rhowindline = lambda t, t0, norm: norm*np.power(t/t0,-1./4.)
		rhoconstline = lambda t, t0, norm: norm*np.power(t/t0,-3./8.)
		tstart = 1e2
		tstop = 1e8
		tnum = 1e2
		g_norm_wind = rhowindline(tstart,fs_data['TA'][np.argmax(fs_data['TA']>tstart)],fs_data['GAMMAR'][np.argmax(fs_data['TA']>tstart)])
		g_norm_const = rhoconstline(tstart,fs_data['TA'][np.argmax(fs_data['TA']>tstart)],fs_data['GAMMAR'][np.argmax(fs_data['TA']>tstart)])
		t = np.linspace(tstart,tstop,num=int(tnum))	
		ax[0,1].plot(t,rhowindline(t,tstart,g_norm_wind),label=r"$t^{-1/4}$",color='r')
		ax[0,1].plot(t,rhoconstline(t,tstart,g_norm_const),label=r"$t^{-3/8}$",color='k')

		plot_param_vs_ta(fs_data,'GAMMAR', ax=ax[0,1], z=z, Tmin=Tmin, Tmax=Tmax,
				fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker=',',color="C1")

		# T_a vs e_diss
		plot_param_vs_ta(fs_data,'EDISS', ax=ax[1,0], z=z, Tmin=Tmin, Tmax=Tmax,
				fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker=',',color="C1")

		# T_a vs E_syn
		plot_param_vs_ta(fs_data,'ESYN', ax=ax[1,1], z=z, Tmin=Tmin, Tmax=Tmax,
				fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker=',',color="C1")

	if(rs_data is not None):
		# T_a vs T_e
		plot_param_vs_ta(rs_data,'TE', ax=ax[0,0], z=z, Tmin=Tmin, Tmax=Tmax,
			fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker='.')	
		# T_a vs del T_a

		# T_a vs Gamma_r
		plot_param_vs_ta(rs_data,'GAMMAR', ax=ax[0,1], z=z, Tmin=Tmin, Tmax=Tmax,
			fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker='.')
		# T_a vs e_diss
		plot_param_vs_ta(rs_data,'EDISS', ax=ax[1,0], z=z, Tmin=Tmin, Tmax=Tmax,
			fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker='.')
		# T_a vs E_syn
		plot_param_vs_ta(rs_data,'ESYN', ax=ax[1,1], z=z, Tmin=Tmin, Tmax=Tmax,
			fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker='.')

	# Plot Aesthetics
	fontsize = 14
	fontweight = "bold"

	# Format Top Left plot, T_a vs T_e
	
	ax[0,0].set_ylabel(r"$t_e$",fontsize=fontsize,fontweight=fontweight)
	ax[0,0].set_xlabel(r"$t_a$",fontsize=fontsize,fontweight=fontweight)

	ax[0,0].set_yscale("log")
	ax[0,0].set_xscale("log")

	# Format Top Right plot, T_a vs Gamma_r 
	ax[0,1].set_ylabel(r"$\Gamma_r$",fontsize=fontsize,fontweight=fontweight)
	ax[0,1].set_xlabel(r"$t_a$",fontsize=fontsize,fontweight=fontweight)

	ax[0,1].yaxis.set_label_position("right")
	ax[0,1].yaxis.tick_right()
	ax[0,1].set_yscale("log")
	ax[0,1].set_xscale("log")
	ax[0,1].legend([r"$\rho_{wind}$, $t^{-1/4}$",r"$\rho_{const}$, $t^{-3/8}$","FS","RS"],fontsize=fontsize-4)


	# Format Bottom Left plot, T_a vs E_diss
	ax[1,0].set_ylabel(r"E$_{diss}$",fontsize=fontsize,fontweight=fontweight)
	ax[1,0].set_xlabel(r"$t_a$",fontsize=fontsize,fontweight=fontweight)
	ax[1,0].set_yscale("log")
	ax[1,0].set_xscale("log")

	# Format Bottom Right plot, T_a vs E_synch
	ax[1,1].set_ylabel(r"E$_{syn}$",fontsize=fontsize,fontweight=fontweight)
	ax[1,1].set_xlabel(r"$t_a$",fontsize=fontsize,fontweight=fontweight)
	ax[1,1].yaxis.set_label_position("right")
	ax[1,1].yaxis.tick_right()
	ax[1,1].set_yscale("log")
	ax[1,1].set_xscale("log")

	# Make plots look good
	for i in range(2):
		for j in range(2):
			plot_aesthetics(ax[i,j],fontsize=fontsize,fontweight=fontweight)

	plt.tight_layout()
	plt.subplots_adjust(wspace=0,hspace=0)

	if save_pref is not None:
		plt.savefig('figs/{}-ext-shock-evo-fig0.png'.format(save_pref))


	fig, ax = plt.subplots(2,1,sharex=True,figsize=(6,6))

	if(fs_data is not None):
		# # T_a vs EPS
		# plot_param_vs_ta(fs_data,'EPS', ax=ax[0,0], z=z, Tmin=Tmin, Tmax=Tmax,
		# 		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker=',',color="C1")

		# # T_a vs rho
		# plot_param_vs_ta(fs_data,'RHO', ax=ax[0,1], z=z, Tmin=Tmin, Tmax=Tmax,
		# 		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker=',',color="C1")

		# T_a vs B_eq
		plot_param_vs_ta(fs_data,'BEQ', ax=ax[0], z=z, Tmin=Tmin, Tmax=Tmax,
				fontsize=fontsize, fontweight=fontweight, disp_xax=False, marker=',',color="C1")

		# T_a vs Gamma_e
		plot_param_vs_ta(fs_data,'GAMMAE', ax=ax[1], z=z, Tmin=Tmin, Tmax=Tmax,
				fontsize=fontsize, fontweight=fontweight, marker=',',color="C1")

	if(rs_data is not None):
		# # T_a vs EPS
		# plot_param_vs_ta(rs_data,'EPS', ax=ax[0,0], z=z, Tmin=Tmin, Tmax=Tmax,
		# 	fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker='.')	

		# # T_a vs rho
		# plot_param_vs_ta(rs_data,'RHO', ax=ax[0,1], z=z, Tmin=Tmin, Tmax=Tmax,
		# 	fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker='.')
		
		# T_a vs B_eq
		plot_param_vs_ta(rs_data,'BEQ', ax=ax[0], z=z, Tmin=Tmin, Tmax=Tmax,
			fontsize=fontsize, fontweight=fontweight, disp_xax=False, marker='.')
		# T_a vs Gamma_e
		plot_param_vs_ta(rs_data,'GAMMAE', ax=ax[1], z=z, Tmin=Tmin, Tmax=Tmax,
			fontsize=fontsize, fontweight=fontweight, marker='.')

	# Plot Aesthetics
	fontsize = 14
	fontweight = "bold"

	# Format Top Left plot, T_a vs Eps
	
	# ax[0,0].set_ylabel(r"$\epsilon$",fontsize=fontsize,fontweight=fontweight)
	# ax[0,0].set_xlabel(r"$t_a$",fontsize=fontsize,fontweight=fontweight)

	# ax[0,0].set_yscale("log")
	# ax[0,0].set_xscale("log")

	# # Format Top Right plot, T_a vs rho
	# ax[0,1].set_ylabel(r"$\rho$",fontsize=fontsize,fontweight=fontweight)
	# ax[0,1].set_xlabel(r"$t_a$",fontsize=fontsize,fontweight=fontweight)

	# ax[0,1].yaxis.set_label_position("right")
	# ax[0,1].yaxis.tick_right()
	# ax[0,1].set_yscale("log")
	# ax[0,1].set_xscale("log")
	# ax[0,1].legend(["FS","RS"],fontsize=fontsize-4)

	# Format Bottom Left plot, T_a vs B_eq
	ax[0].set_ylabel(r"B$_{EQ}$",fontsize=fontsize,fontweight=fontweight)
	ax[0].set_xlabel(r"$t_a$",fontsize=fontsize,fontweight=fontweight)
	ax[0].set_yscale("log")
	ax[0].set_xscale("log")

	# Format Bottom Right plot, T_a vs Gamma_e
	ax[1].set_ylabel(r"$\Gamma_{e}$",fontsize=fontsize,fontweight=fontweight)
	ax[1].set_xlabel(r"$t_a$",fontsize=fontsize,fontweight=fontweight)
	# ax[1,1].yaxis.set_label_position("right")
	# ax[1,1].yaxis.tick_right()
	ax[1].set_yscale("log")
	ax[1].set_xscale("log")

	# Make plots look good
	for i in range(2):
		plot_aesthetics(ax[i],fontsize=fontsize,fontweight=fontweight)

	plt.tight_layout()
	plt.subplots_adjust(wspace=0,hspace=0)

	if save_pref is not None:
		plt.savefig('figs/{}-ext-shock-evo-fig1.png'.format(save_pref))

##############################################################################################################################

def plot_gamma_r(is_data = None,fs_data=None, rs_data=None, ax=None, z=0, Tmin=None, Tmax=None,save_pref=None,fontsize=14,fontweight='bold'):

	# Make plot instance if it doesn't exist
	if ax is None:
		ax = plt.figure().gca()

	if is_data is not None:
		plot_param_vs_ta(is_data,'GAMMAR', ax=ax, y_factor=1/100, z=z,Tmin=Tmin, Tmax=Tmax,fontsize=fontsize, fontweight=fontweight,marker='.',color='C0',label="IS")
	if fs_data is not None:
		plot_param_vs_ta(fs_data,'GAMMAR', ax=ax, y_factor=1/100, z=z,Tmin=Tmin, Tmax=Tmax,fontsize=fontsize, fontweight=fontweight,marker='.',color='C1',label="FS")
	if rs_data is not None:
		plot_param_vs_ta(rs_data,'GAMMAR', ax=ax, y_factor=1/100, z=z,Tmin=Tmin, Tmax=Tmax,fontsize=fontsize, fontweight=fontweight,marker='.',color='C2',label="RS")

	
	ax.set_xscale("log")
	ax.set_yscale("log")

	ax.set_xlabel(r't$_{a}$ (sec)',fontsize=fontsize,fontweight=fontweight)
	ax.set_ylabel(r'Gamma$_{r}$/100',fontsize=fontsize,fontweight=fontweight)

	ax.legend(fontsize=fontsize-4)

	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)
	
	plt.tight_layout()

	if save_pref == True:
		plt.savefig('figs/{}-gamma_r.png'.format(save_pref))

##############################################################################################################################

def lum_dis(z: float):
	""" 
	Caclulate luminosity distance for redshift z
	"""
	if(z == 0):
		return 1
	else:
		# bol_lum = [1,100000] # bolumetric luminosity range
		c = 3*np.power(10,10) # speed of light, cm/s
		omega_m = 0.3 # matter density of the universe
		omega_lam = 0.7 # dark energy density of the universe
		H0 = 67.4*np.power(10,5) # Hubbles Constant cm/s/Mpc

		lum_dis_Mpc = ((1+z)*c/(H0) ) * integrate.quad(lambda zi: 1/np.sqrt( ((omega_m*np.power(1+zi,3) )+omega_lam) ),0,z)[0]
		lum_dis_cm = lum_dis_Mpc * 3.086e24 # Mpc -> cm
		return lum_dis_cm

##############################################################################################################################


if __name__ == '__main__':

	z = 1

	"""
	Shell Lorentz Distribution
	"""
	"""
	ax_sd = plt.figure().gca()
	plot_lor_dist('data-file-dir/shell_dist.txt', ax=ax_sd)
	"""

	"""
	Synthetic spectrum 
	"""
	
	
	ax_spec = plt.figure(figsize=(9,8)).gca()

	plot_spec("data-file-dir/test_spec_therm.txt",ax=ax_spec,z=z,label="Thermal")
	plot_spec("data-file-dir/test_spec_is.txt",ax=ax_spec,z=z,label="IS")
	plot_spec("data-file-dir/test_spec_fs.txt",ax=ax_spec,z=z,label="FS")
	plot_spec("data-file-dir/test_spec_rs.txt",ax=ax_spec,z=z,label="RS")
	plot_spec("data-file-dir/test_spec_total.txt",ax=ax_spec,z=z,label="Total")

	## Synthetic spectrum before convolusion
	# plot_spec("data-file-dir/spec_source.txt",ax=ax_spec,unc=False,label="Source")
	# plot_spec("data-file-dir/spec_source_fluc.txt",ax=ax_spec,unc=True,label="Pre-Conv")
	# plot_spec("data-file-dir/spec_model.txt",ax=ax_spec,unc=False,label="Model",joined=True)
	
	## Synthetic spectrum after convolusion
	# plot_spec("data-file-dir/spec_obs.txt",ax=ax_spec,unc=True,label="Obs")
	# plot_spec("data-file-dir/spec_model_conv.txt",ax=ax_spec,unc=False,label="Model",joined=True)

	add_FermiGBM_band(ax_spec)
	# add_SwiftBAT_band(ax_spec)

	# ax_spec.set_xlim(0.1,1e5)
	# ax_spec.set_ylim(1e48,1e52)
	
	

	"""
	Synthetic light curve
	"""
	
	ax_lc = plt.figure().gca()
	plot_light_curve("data-file-dir/test_light_curve.txt",ax=ax_lc,z=z,label="Total")
	# plot_light_curve("data-file-dir/test_light_curve_therm.txt",ax=ax_lc,z=0.5,label="Therm")
	# plot_light_curve("data-file-dir/test_light_curve_is.txt",ax=ax_lc,z=0.5,label="IS")

	# Interactive light curve
	# tbox = plot_light_curve_interactive("data-file-dir/test_light_curve.txt",z=z,label="Total")
	

	"""
	Jet dynamics plots 
	"""
	
	
	# therm_emission = load_therm_emission("data-file-dir/synthGRB_jet_params_therm.txt")
	# plot_evo_therm(therm_emission,z=1)

	is_emission = load_is_emission("data-file-dir/synthGRB_jet_params_is.txt")
	plot_evo_int_shock(is_emission,z=0.5)

	fs_data = load_fs_emission("data-file-dir/synthGRB_jet_params_fs.txt")
	rs_data = load_rs_emission("data-file-dir/synthGRB_jet_params_rs.txt")

	plot_evo_ext_shock(fs_data=fs_data,rs_data=rs_data,z=0)
	# plot_evo_ext_shock(fs_data=fs_data)		
	
	
	"""
	Display real observed data
	"""
	"""
	ax_spec = plt.figure(figsize=(8,4)).gca()
	plot_spec("data-file-dir/190114C_n4_tte_spec_bak.txt",ax=ax_spec,unc=False,label="BGD")
	plot_spec("data-file-dir/190114C_n4_tte_spec_rise.txt",ax=ax_spec,unc=False,label="Rise")
	plot_spec("data-file-dir/190114C_n4_tte_spec_peak.txt",ax=ax_spec,unc=False,label="Peak")
	plot_spec("data-file-dir/190114C_n4_tte_spec_fall.txt",ax=ax_spec,unc=False,label="Fall")

	ax_lc = plt.figure().gca()
	plot_light_curve("data-file-dir/190114C_n4_tte_light_curve.txt",ax=ax_lc)
	"""

	"""
	Multiple Synthetic spectrum 
	"""
	
	"""
	## Model spectrum before convolusion
	ax_spec_preconv = plt.figure(figsize=(8,4)).gca()
	plot_spec("data-file-dir/spec_model_1.txt",ax=ax_spec_preconv,unc=True,label="NaI")
	plot_spec("data-file-dir/spec_model_2.txt",ax=ax_spec_preconv,unc=True,label="BGO")

	add_FermiGBM_band(ax_spec_preconv)

	## Synthetic spectrum after convolusion
	ax_spec_postconv = plt.figure(figsize=(8,4)).gca()
	plot_spec("data-file-dir/spec_obs_1.txt",ax=ax_spec_postconv,unc=True,label="NaI")
	plot_spec("data-file-dir/spec_obs_2.txt",ax=ax_spec_postconv,unc=True,label="BGO")
	plot_spec("data-file-dir/spec_model_conv_1.txt",ax=ax_spec_postconv,unc=False,label="Model NaI",joined=True)
	plot_spec("data-file-dir/spec_model_conv_2.txt",ax=ax_spec_postconv,unc=False,label="Model BGO",joined=True)

	add_FermiGBM_band(ax_spec_postconv)
	"""


	"""
	Testing
	"""









	plt.show()
