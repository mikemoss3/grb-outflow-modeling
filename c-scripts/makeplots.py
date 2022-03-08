"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-11-26


Meta script to plot desired simulation results created by c++ code.

"""
import matplotlib.pyplot as plt
import numpy as np
import subprocess
# import cosmologicalconstants as cc
import scipy.integrate as integrate 
from matplotlib.widgets import TextBox
from matplotlib.lines import Line2D
from matplotlib.animation import FuncAnimation, PillowWriter 

kb_kev = 8.617*1e-8
kev_to_erg = 1.6022*np.power(10.,-9.)
planck_kev = 4.136 * np.power(10,-18.) # keV Hz^-1

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
	lor_dist_list = load_lor_dist(file_name)

	# Make plot instance if it doesn't exist
	if ax is None:
		fig, ax = plt.subplots(1,2,figsize=(16, 6))	

	ind = np.array([0])
	def on_press(event):
		if event.key == 'right':
			ind[0] +=1
		elif event.key == 'left':
			ind[0] -=1
		# Update shell index plot 
		line_shell_ind.set_xdata(np.linspace(start=0,stop=len(lor_dist_list[ind[0]]['GAMMA']),num=len(lor_dist_list[ind[0]]['GAMMA'])))
		line_shell_ind.set_ydata(lor_dist_list[ind[0]]['GAMMA'])

		# Update cumulative mass fraction plot 
		flipped_mass_arr = np.flip(lor_dist_list[ind[0]]['MASS'])
		flipped_gamma_arr = np.flip(lor_dist_list[ind[0]]['GAMMA'])
		# Cumulative mass
		masscum = np.cumsum(flipped_mass_arr)
		massfraccum = masscum/masscum[-1]
		line_mass_frac.set_xdata(massfraccum)
		line_mass_frac.set_ydata(flipped_gamma_arr)

		line_mass_frac_zoom.set_xdata(massfraccum)
		line_mass_frac_zoom.set_ydata(flipped_gamma_arr)

		# Redraw the figure to implement updates
		ax[0].redraw_in_frame()
		ax[1].redraw_in_frame()
		fig.canvas.draw_idle()

	fig.canvas.mpl_connect('key_press_event', on_press)

	# Plot distribution as a function of the shell number 
	line_shell_ind, = ax[0].step(np.linspace(start=0,stop=len(lor_dist_list[0]['GAMMA']),num=len(lor_dist_list[0]['GAMMA'])),lor_dist_list[0]['GAMMA'],where='pre',linestyle=linestyle,label=label)
	ax[0].invert_xaxis()
	ax[0].set_ylim(0)

	# Plot the Lorentz distribution as a function of the mass fraction
	flipped_mass_arr = np.flip(lor_dist_list[0]['MASS'])
	flipped_gamma_arr = np.flip(lor_dist_list[0]['GAMMA'])

	# Cumulative mass
	masscum = np.cumsum(flipped_mass_arr)
	massfraccum = masscum/masscum[-1]

	# Plot distribution as a function of the mass fraction
	line_mass_frac, = ax[1].step(massfraccum,flipped_gamma_arr,where='pre',linestyle=linestyle,label=label)
	
	ax[1].set_ylim(0)

	# Zoomed in version of the distribution as a function of the mass fraction
	ax_zoom = ax[1].twinx()
	line_mass_frac_zoom, = ax_zoom.step(massfraccum,flipped_gamma_arr,where='pre',linestyle=linestyle,label=label,color="C1")
	ax_zoom.set_ylim(0,50)

	if xlabel is True:
		ax[0].set_xlabel('Shell Ind.',fontsize=fontsize,fontweight=fontweight)
		ax[1].set_xlabel(r'M/M$_{tot}$',fontsize=fontsize,fontweight=fontweight)
	if ylabel is True:
		ax[0].set_ylabel(r'$\Gamma$',fontsize=fontsize,fontweight=fontweight)
		ax[1].set_ylabel(r'$\Gamma$',fontsize=fontsize,fontweight=fontweight)

	plot_aesthetics(ax[0],fontsize=fontsize,fontweight=fontweight)
	plot_aesthetics(ax[1],fontsize=fontsize,fontweight=fontweight)
	plot_aesthetics(ax_zoom,fontsize=fontsize,fontweight=fontweight)
	if label is not None:
		ax[0].legend(fontsize=fontsize)

	if save_pref is not None :
		plt.savefig('figs/{}-lorentz-dist.png'.format(save_pref))

##############################################################################################################################

def plot_lor_dist_anim(file_name,ax=None,save_pref=None,xlabel=True,ylabel=True,label=None,fontsize=14,fontweight='bold',linestyle='solid'):
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
	lor_dist_list = load_lor_dist(file_name)

	# Make plot instance if it doesn't exist
	fig, ax = plt.subplots(1,2,figsize=(16, 6))	

	# Plot distribution as a function of the shell number 
	line_shell_ind, = ax[0].step(np.linspace(start=0,stop=len(lor_dist_list[0]['GAMMA']),num=len(lor_dist_list[0]['GAMMA'])),lor_dist_list[0]['GAMMA'],where='pre',linestyle=linestyle,label=label)
	ax[0].invert_xaxis()
	ax[0].set_ylim(0)

	# Plot the Lorentz distribution as a function of the mass fraction
	flipped_mass_arr = np.flip(lor_dist_list[0]['MASS'])
	flipped_gamma_arr = np.flip(lor_dist_list[0]['GAMMA'])

	# Cumulative mass
	masscum = np.cumsum(flipped_mass_arr)
	massfraccum = masscum/masscum[-1]

	# Plot distribution as a function of the mass fraction
	line_mass_frac, = ax[1].step(massfraccum,flipped_gamma_arr,where='pre',linestyle=linestyle,label=label)
	ax[1].set_ylim(0)

	# Zoomed in version of the distribution as a function of the mass fraction
	ax_zoom = ax[1].twinx()
	line_mass_frac_zoom, = ax_zoom.step(massfraccum,flipped_gamma_arr,where='pre',linestyle=linestyle,label=label,color="C1")
	ax_zoom.set_ylim(0,50)

	if xlabel is True:
		ax[0].set_xlabel('Shell Ind.',fontsize=fontsize,fontweight=fontweight)
		ax[1].set_xlabel(r'M/M$_{tot}$',fontsize=fontsize,fontweight=fontweight)
	if ylabel is True:
		ax[0].set_ylabel(r'$\Gamma$',fontsize=fontsize,fontweight=fontweight)
		ax[1].set_ylabel(r'$\Gamma$',fontsize=fontsize,fontweight=fontweight)

	plot_aesthetics(ax[0],fontsize=fontsize,fontweight=fontweight)
	plot_aesthetics(ax[1],fontsize=fontsize,fontweight=fontweight)
	plot_aesthetics(ax_zoom,fontsize=fontsize,fontweight=fontweight)
	if label is not None:
		ax[0].legend(fontsize=fontsize)

	ind = np.array([0])

	def update(i):

		ind[0] = i
		
		# Update shell index plot 
		line_shell_ind.set_xdata(np.linspace(start=0,stop=len(lor_dist_list[ind[0]]['GAMMA']),num=len(lor_dist_list[ind[0]]['GAMMA'])))
		line_shell_ind.set_ydata(lor_dist_list[ind[0]]['GAMMA'])

		# Update cumulative mass fraction plot 
		flipped_mass_arr = np.flip(lor_dist_list[ind[0]]['MASS'])
		flipped_gamma_arr = np.flip(lor_dist_list[ind[0]]['GAMMA'])
		# Cumulative mass
		masscum = np.cumsum(flipped_mass_arr)
		massfraccum = masscum/masscum[-1]
		line_mass_frac.set_xdata(massfraccum)
		line_mass_frac.set_ydata(flipped_gamma_arr)

		line_mass_frac_zoom.set_xdata(massfraccum)
		line_mass_frac_zoom.set_ydata(flipped_gamma_arr)

		# Redraw the figure to implement updates
		# ax[0].redraw_in_frame()
		# ax[1].redraw_in_frame()
		fig.canvas.draw_idle()

	ani = FuncAnimation(fig=fig, func=update, frames=len(lor_dist_list) ,interval=400)

	if save_pref is not None :
		ani.save('figs/{}-lorentz-dist-anim.gif'.format(save_pref))

	return ani


##############################################################################################################################

def load_lor_dist(file_name, string_match = "// Next step\n"):
	"""
	Method to load Lorentz factor distribution
	"""
	
	# Function to find all line numbers that match the "string_match" input argument
	def lines_that_equal(string, fp):
		line_num_list = [] # List to store line numbers
		line_num = 0 # Temporary line number holder
		for line in fp:
			line_num +=1
			if line == string:
				# If the line matches, append the line number to the line number list
				line_num_list.append(line_num)
		line_num_list.append(line_num) # Append last line, this is used to indicate end of file when loading data
		return line_num_list

	# With the input file open, use the lines_that_equal function to find line numbers that match the "string_match" argument
	line_num_list = 0
	with open(file_name, "r") as fp:
		line_num_list = lines_that_equal(string_match,fp)

	# For all the data between each line number in the line_num_list, load the data and append it to the data_list
	data_list = []
	for i in range(len(line_num_list)-1):
		tmp_data = np.genfromtxt(file_name,skip_header=line_num_list[i],max_rows=line_num_list[i+1]-1-line_num_list[i], dtype=[("RADIUS",float),("GAMMA",float),("MASS",float),("TE",float),("STATUS",float)])
		data_list.append(tmp_data[tmp_data['STATUS']>0])
		# data_list.append(tmp_data)

	return data_list

##############################################################################################################################

def plot_spec(file_name, z=0, joined=False, label = None, color="C0", ax=None, nuFnu=True, unc=False, Emin=None, Emax=None, save_pref=None,fontsize=14,fontweight='bold'):
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
				line = ax.errorbar(x=spec_data['ENERG'],y=spec_data['RATE']*(spec_data['ENERG']**2),yerr=spec_data['UNC']*(spec_data['ENERG']**2),label=label,color=color)
			else:
				line, = ax.plot(spec_data['ENERG'],spec_data['RATE']*(spec_data['ENERG']**2),label=label,color=color)
		else:
			if unc is True:
				line = ax.errorbar(x=spec_data['ENERG'],y=spec_data['RATE'],yerr=spec_data['UNC'],label=label,color=color)
			else:
				line, = ax.plot(spec_data['ENERG'],spec_data['RATE'],label=label,color=color)
	else:
		# Plot spectrum data
		if nuFnu is True:
			if unc is True:
				line = ax.errorbar(x=spec_data['ENERG'],y=spec_data['RATE']*(spec_data['ENERG']**2),yerr=spec_data['UNC']*(spec_data['ENERG']**2),label=label,fmt=" ",marker="+",color=color)
			else:
				line = ax.errorbar(x=spec_data['ENERG'],y=spec_data['RATE']*(spec_data['ENERG']**2),label=label,fmt=" ",marker="+",color=color)
		else:
			if unc is True:
				line = ax.errorbar(x=spec_data['ENERG'],y=spec_data['RATE'],yerr=spec_data['UNC'],label=label,fmt=" ",marker="+",color=color)
			else:
				line = ax.errorbar(x=spec_data['ENERG'],y=spec_data['RATE'],label=label,fmt=" ",marker="+",color=color)

	# Plot aesthetics
	ax.set_xscale('log')
	ax.set_yscale('log')

	# Force lower bound
	ax.set_ylim(1e48,1e51)

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

def add_FermiGBM_band(ax,fontsize=12,axis="x"):
	# Grab the current ymin and ymax, this is used to set the lower and upper bounds of the vertical lines which indicate instrument observation energy range
	curr_ymin, curr_ymax = ax.get_ylim()
	curr_xmin, curr_xmax = ax.get_xlim()

	# Vertical axis
	if(axis == "x"):
		# Display Fermi/GBM - NAI energy band
		ax.axvspan(xmin=8,xmax=1e3,ymin=0.5,alpha=0.4,facecolor='grey',label='Fermi/GBM-NAI')

		# Display Fermi/GBM - BGO energy band
		ax.axvspan(xmin=150,xmax=3*1e4,ymin=0.5,alpha=0.4,facecolor='orange',label='Fermi/GBM-BGO')

	# Horizontal axis
	elif(axis == "y"):
		# Display Fermi/GBM - NAI energy band
		ax.axhspan(ymin=8,ymax=1e3,alpha=0.4,facecolor='grey',label='Fermi/GBM-NAI')

		# Display Fermi/GBM - BGO energy band
		ax.axhspan(ymin=150,ymax=3*1e4,alpha=0.4,facecolor='orange',label='Fermi/GBM-BGO')


	# Add to legend	
	ax.legend(fontsize=fontsize)

	# We don't want the plotting window to change if either of the energy band edges do not overlap with the plotted energy spectra
	ax.set_ylim(curr_ymin,curr_ymax)
	ax.set_xlim(curr_xmin,curr_xmax)

##############################################################################################################################

def add_SwiftBAT_band(ax,fontsize=12,axis="x"):
	# Grab the current ymin and ymax, this is used to set the lower and upper bounds of the vertical lines which indicate instrument observation energy range
	curr_ymin, curr_ymax = ax.get_ylim()
	curr_xmin, curr_xmax = ax.get_xlim()

	# Display Swift/BAT energy band

	# Vertical axis
	if(axis == "x"):
		ax.axvspan(xmin=5,xmax=350,ymin=0.5,alpha=0.4,facecolor='blue',label='Swift/BAT')

	# Horizontal axis
	elif(axis == "y"):
		ax.axhspan(ymin=5,ymax=350,alpha=0.4,facecolor='blue',label='Swift/BAT')

	# Add to legend	
	ax.legend(fontsize=fontsize)

	# We don't want the plotting window to change if either of the energy band edges do not overlap with the plotted energy spectra
	ax.set_ylim(curr_ymin,curr_ymax)
	ax.set_xlim(curr_xmin,curr_xmax)

##############################################################################################################################

def plot_light_curve(file_name, z=0, label=None, ax=None, Tmin=None, Tmax=None, save_pref=None,color="C0", fontsize=14,fontweight='bold', logscale=False):
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
			ax.step(light_curve_data['TIME']*(1+z),light_curve_data['RATE']/(4*np.pi*lum_dis(z)**2),label=label,marker=" ",where="mid",color=color	)
		else: 
			# If z = 0, return luminosity
			# ax.scatter(light_curve_data['TIME'],light_curve_data['RATE'],label=label,marker=".")
			ax.step(light_curve_data['TIME'],light_curve_data['RATE'],label=label,marker=" ",where="mid",color=color	)

		if(logscale == True):
			ax.set_yscale('log')
			ax.set_xscale('log')

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

def plot_light_curve_interactive(init_Tmin, init_Tmax, init_Emin, init_Emax, z=0, 
	with_comps=False, label=None, ax=None, save_pref=None, fontsize=14,fontweight='bold',logscale=False):
	"""
	Method to plot the input light curve data files (interactively!)

	init_Tmin, init_Tmax 
			= Sets the initial time selection interval for calculating the spectrum and also sets the initial time window for the light curve.
			In addition, this will set the largest time interval the light curve will be displayed in.
	init_Emin, init_Emax
			= Sets the initial energy selection interval for calculating the light curve and also sets the initial energy window for the spectrum.
			In addition, this will set the largest energy interval the spectrum will be displayed in.

	z = float > 0, redshift
	with_comps = boolean, indicates whether the individual components should be displayed
	"""

	comp_indicator = "false"
	if(with_comps == True):
		comp_indicator = "true"
	num_comps = 4
	comp_suff = ["TH","IS","FS","RS"]
	comp_colors = ['r','C0','C1','C2']

	if(z<0):
		print("Please provide a non-negative redshift.")
		return;
	
	# Make initial data
	subprocess.run(["./main","timechange", comp_indicator, "{}".format(init_Tmin), "{}".format(init_Tmax), "{}".format(init_Emin), "{}".format(init_Emax) ])
	subprocess.run(["./main","energychange", comp_indicator, "{}".format(init_Tmin), "{}".format(init_Tmax), "{}".format(init_Emin), "{}".format(init_Emax) ])

	# Make plot instance if it doesn't exist
	fig, ax = plt.subplots(1,2,figsize=(16, 8))	
	
	# Load initial light curve data
	lc_data_tot = np.genfromtxt("data-file-dir/quickplot_light_curve.txt",dtype=[("TIME",float),("RATE",float)])

	# Plot initial light curve
	if(z>0):
		lc_tot_line, = ax[0].step(lc_data_tot['TIME']*(1+z),lc_data_tot['RATE']/(4*np.pi*lum_dis(z)**2),label=label,marker=" ",where="mid",color="k")
		ax[0].set_ylabel(r'Rate (ph cm$^{-2}$ s$^{-1}$)',fontsize=fontsize,fontweight=fontweight)
	else: 
		# If z = 0, return luminosity
		lc_tot_line, = ax[0].step(lc_data_tot['TIME'],lc_data_tot['RATE'],label=label,marker=" ",where="mid",color="k")
		ax[0].set_ylabel(r'Rate (ph s$^{-1}$)',fontsize=fontsize,fontweight=fontweight)

	ax[0].set_xlabel('Obs Time (sec)',fontsize=fontsize,fontweight=fontweight)
	plot_aesthetics(ax[0],fontsize=fontsize,fontweight=fontweight)

	# If components are desired
	if(with_comps == True):
		# Load component data
		lc_comp_data = [0] * num_comps
		for i in range(num_comps):
			lc_comp_data[i] = np.genfromtxt("data-file-dir/quickplot_light_curve_{}.txt".format(comp_suff[i]),dtype=[("TIME",float),("RATE",float)])

		# Plot component light curves
		lc_comp_lines = [0] * num_comps
		if(z>0):
			for i in range(num_comps):
				lc_comp_lines[i], = ax[0].step(lc_comp_data[i]['TIME']*(1+z),lc_comp_data[i]['RATE']/(4*np.pi*lum_dis(z)**2),color=comp_colors[i],marker=" ",where="mid")
		else: 
			# If z = 0, return luminosity
			for i in range(num_comps):
				lc_comp_lines[i], = ax[0].step(lc_comp_data[i]['TIME'],lc_comp_data[i]['RATE'],color=comp_colors[i],marker=" ",where="mid")



	# Plot initial spectrum 
	spec_tot_line = plot_spec("data-file-dir/quickplot_spectrum.txt",ax=ax[1],z=z,joined=True,color='k',label='TOT')
	# If components are desired, plot initial component spectra 
	if(with_comps == True):
		spec_comp_lines = [0] * num_comps
		for i in range(num_comps):
			spec_comp_lines[i] = plot_spec("data-file-dir/quickplot_spectrum_{}.txt".format(comp_suff[i]),ax=ax[1],z=z,joined=True,color=comp_colors[i],label=comp_suff[i])
	add_FermiGBM_band(ax[1])

	# Display time selection interval
	selected_Tmin_line = ax[0].axvline(x=init_Tmin,color='k')
	selected_Tmax_line = ax[0].axvline(x=init_Tmax,color='k')

	# Display energy selection interval on plot
	selected_Emin_line = ax[1].axvline(x=init_Emin,color='k')
	selected_Emax_line = ax[1].axvline(x=init_Emax,color='k')

	# Make time selection and time window text boxes
	fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)
	
	axbox_time_select = plt.axes([0.2, 0.08, 0.2, 0.04])
	txtbox_time_select = TextBox(ax=axbox_time_select, label='Time Selection', initial="{0:.2f}, {1:.2f}".format(init_Tmin, init_Tmax) )

	axbox_time_window = plt.axes([0.2, 0.035, 0.2, 0.04])
	txtbox_time_window = TextBox(ax=axbox_time_window, label='Time Window', initial="{0:.2f}, {1:.2f}".format(init_Tmin, init_Tmax) )

	# Make energy selection and energy window text boxes
	axbox_energy_select = plt.axes([0.65, 0.08, 0.2, 0.04])
	txtbox_energy_select = TextBox(ax=axbox_energy_select, label='Energy ', initial="{0:.2f}, {1:.2f}".format(init_Emin,init_Emax) )

	axbox_energy_window = plt.axes([0.65, 0.035, 0.2, 0.04])
	txtbox_energy_window = TextBox(ax=axbox_energy_window, label='Energy ', initial="{0:.2f}, {1:.2f}".format(init_Emin,init_Emax) )


	# Make update functions
	
	def submit_time_select(val):
		"""
		Function call when a time selection is given, this will alter the displayed spectra
		"""
		# val will be separated into two values, split by the comma
		x = val.split(", ")

		selected_Tmin = float(x[0])
		selected_Tmax = float(x[1])

		# Update the time selection indicating lines
		selected_Tmin_line.set_xdata(selected_Tmin)
		selected_Tmax_line.set_xdata(selected_Tmax)

		# Compute the new spectrum for this time selection
		subprocess.run(["./main","timechange", comp_indicator, "{}".format(selected_Tmin), "{}".format(selected_Tmax), "{}".format(init_Emin), "{}".format(init_Emax) ])

		# If components are desired, plot initial component spectra 
		if(with_comps == True):
			for i in range(num_comps):
				# Load component spectra data
				tmp_comp_data = np.genfromtxt("data-file-dir/quickplot_spectrum_{}.txt".format(comp_suff[i]),dtype=[("ENERGY",float),("RATE",float),("UNC",float)])
				# Update plotted data 
				spec_comp_lines[i].set_xdata(tmp_comp_data['ENERGY']/(1+z))
				spec_comp_lines[i].set_ydata(tmp_comp_data['RATE'] * tmp_comp_data['ENERGY']**2 )

		# Load spectrum data
		tmp_comp_data = np.genfromtxt("data-file-dir/quickplot_spectrum.txt",dtype=[("ENERGY",float),("RATE",float),("UNC",float)])
		# Update plotted spectrum data
		spec_tot_line.set_xdata(tmp_comp_data['ENERGY']/(1+z))
		spec_tot_line.set_ydata(tmp_comp_data['RATE']* tmp_comp_data['ENERGY']**2)

		# Redraw the figure to implement updates
		ax[0].redraw_in_frame()
		ax[1].redraw_in_frame()
		fig.canvas.draw_idle()

		
	def submit_energy_select(val):
		"""
		Function call when a energy selection is given, this will edit the observed light curve 
		"""
		# val will be separated into two values, split by the comma
		x = val.split(", ")

		selected_Emin = float(x[0])
		selected_Emax = float(x[1])

		# Update the energy selection indicating lines
		selected_Emin_line.set_xdata(selected_Emin)
		selected_Emax_line.set_xdata(selected_Emax)

		# Compute the new light curve for this energy selection
		subprocess.run(["./main","energychange", comp_indicator, "{}".format(init_Tmin), "{}".format(init_Tmax), "{}".format(selected_Emin), "{}".format(selected_Emax) ])

		# If components are desired, plot initial component spectra 
		if(with_comps == True):
			for i in range(num_comps):
				# Load component light curve data
				tmp_comp_data = np.genfromtxt("data-file-dir/quickplot_light_curve_{}.txt".format(comp_suff[i]),dtype=[("TIME",float),("RATE",float)])
				# Update plotted data 
				if(z>0):
					lc_comp_lines[i].set_xdata(tmp_comp_data['TIME']*(1+z))
					lc_comp_lines[i].set_ydata(tmp_comp_data['RATE']/(4*np.pi*lum_dis(z)**2))
				else:
					lc_comp_lines[i].set_xdata(tmp_comp_data['TIME'])
					lc_comp_lines[i].set_ydata(tmp_comp_data['RATE'])
		# Load light curve data
		tmp_comp_data = np.genfromtxt("data-file-dir/quickplot_light_curve.txt",dtype=[("TIME",float),("RATE",float)])
		# Update plotted spectrum data
		if(z>0):
			lc_tot_line.set_xdata(tmp_comp_data['TIME']*(1+z))
			lc_tot_line.set_ydata(tmp_comp_data['RATE']/(4*np.pi*lum_dis(z)**2))
		else:
			lc_tot_line.set_xdata(tmp_comp_data['TIME'])
			lc_tot_line.set_ydata(tmp_comp_data['RATE'])


		# Redraw the figure to implement updates
		ax[0].redraw_in_frame()
		ax[1].redraw_in_frame()
		fig.canvas.draw_idle()

	def submit_time_window(val):
		"""
		Function call when a time window is given
		"""
		# val will be separated into two values, split by the comma
		x = val.split(", ")

		window_Tmin = float(x[0])
		window_Tmax = float(x[1])

		# Update the time window
		ax[0].set_xlim(window_Tmin,window_Tmax)

		# Redraw the figure to implement updates
		ax[0].redraw_in_frame()
		ax[1].redraw_in_frame()
		fig.canvas.draw_idle()
		
	def submit_energy_window(val):
		"""
		Function call when a energy window is given
		"""
		# val will be separated into two values, split by the comma
		x = val.split(", ")

		window_Emin = float(x[0])
		window_Emax = float(x[1])

		# Update the energy window
		ax[1].set_xlim(window_Emin,window_Emax)

		# Redraw the figure to implement updates
		ax[0].redraw_in_frame()
		ax[1].redraw_in_frame()
		fig.canvas.draw_idle()

		

	# Call functions to update selection intervals
	txtbox_time_select.on_submit(submit_time_select)
	txtbox_energy_select.on_submit(submit_energy_select)
	
	# Call functions to update windows
	txtbox_time_window.on_submit(submit_time_window)
	txtbox_energy_window.on_submit(submit_energy_window)

	fig.subplots_adjust(bottom=0.25,wspace=0.4)

	# plt.tight_layout()
	if save_pref is not None:
		plt.savefig('figs/{}-light-curve.png'.format(save_pref))

	return txtbox_time_select, txtbox_time_window, txtbox_energy_select, txtbox_energy_window


##############################################################################################################################

def load_therm_emission(file_name):
	"""
	Method to load thermal emission data from the given file name
	"""

	dtype = np.dtype([('TE',float),('TA',float),('DELT',float),('TEMP',float),('FLUX',float),('RPHOT',float)])

	return np.genfromtxt(file_name,dtype=dtype)

##############################################################################################################################

def load_is_emission(file_name):
	"""
	Method to load internal shock emission data from the given file name
	"""

	dtype = np.dtype([('TE',float),('TA',float),('DELT',float),('BEQ',float),('GAMMAE',float),('ESYN',float),('GAMMAR',float),('EDISS',float),("NUC",float),("NUM",float),("SHIND",int),('ASYN',float),('TAU',float),('RELVEL',float)])

	return np.genfromtxt(file_name,dtype=dtype)

##############################################################################################################################

def load_fs_emission(file_name):
	"""
	Method to load forward shock emission data from the given file name
	"""

	dtype = np.dtype([('TE',float),('TA',float),('DELT',float),('BEQ',float),('GAMMAE',float),('ESYN',float),('GAMMAR',float),('EDISS',float),("NUC",float),("NUM",float),("SHIND",int)])

	return np.genfromtxt(file_name,dtype=dtype)

##############################################################################################################################

def load_rs_emission(file_name):
	"""
	Method to load reverse shock emission data from the given file name
	"""

	dtype = np.dtype([('TE',float),('TA',float),('DELT',float),('BEQ',float),('GAMMAE',float),('ESYN',float),('GAMMAR',float),('EDISS',float),("NUC",float),("NUM",float),("SHIND",int)])

	return np.genfromtxt(file_name,dtype=dtype)

##############################################################################################################################

def plot_param_vs_time(emission_comp,param,ax=None,z=0, y_factor=1, label=None, Tmin=None, Tmax=None,save_pref=None,fontsize=14,fontweight='bold',disp_xax=True,disp_yax=True,
	color='C0',marker='.',alpha=1,frame="obs",markersize=7):
	"""
	Plot emission parameters as a function of time (in the observer frame)
	"""

	time_str = "TA"
	if(frame=="source"):
		time_str = "TE"

	# Make plot instance if it doesn't exist
	if ax is None:
		ax = plt.figure().gca()

	# Multiply by 1+z for the time axis and apply the supplied factor on the y-axis 
	ax_time = emission_comp[time_str] * (1+z)

	# Multiply by input factor
	ax_param = emission_comp[param]*y_factor

	# Find the indices of the start and stop time
	ind_start, ind_stop = 0,-1
	if Tmin is not None:
		ind_start = np.argmax(ax_time>Tmin)
	if(Tmax is not None):
		if ( Tmax < np.max(emission_comp[time_str]) ):
			ind_stop = np.argmax(ax_time>Tmax)

			# Load time axis and the parameter value axis 
			ax_time = ax_time[ind_start:ind_stop]
			ax_param = ax_param[ind_start:ind_stop]
		else:
			# Load time axis and the parameter value axis 
			ax_time = ax_time[ind_start:] 
			ax_param = ax_param[ind_start:]
	else:
		# Load time axis and the parameter value axis 
		ax_time = ax_time[ind_start:] 
		ax_param = ax_param[ind_start:]

	ax.scatter(x=ax_time,y=ax_param,label=label,c=color,marker=marker,s=markersize, picker=True,alpha=alpha)

	if disp_yax is True:
		ax.set_ylabel(param,fontsize=fontsize,fontweight=fontweight)
	if disp_xax is True:
		if(time_str == "TA"):
			ax.set_xlabel(r't$_a$',fontsize=fontsize,fontweight=fontweight)
		elif(time_str == "TE"):
			ax.set_xlabel(r't$_e$',fontsize=fontsize,fontweight=fontweight)

	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)
		
	if save_pref is not None:
		plt.savefig('figs/{}-param-{}-vs-t.png'.format(save_pref,param))

##############################################################################################################################

def plot_evo_therm(thermal_emission,ax=None,z=0,Tmin=None, Tmax=None,save_pref=None,fontsize=14,fontweight='bold',frame="obs"):
	"""
	Plot evolution of thermal emission parameters 
	"""

	# Make plot instance if it doesn't exist
	if ax is None:
		fig, ax = plt.subplots(2,1,figsize=(5,8))

	# Plot temperature of the thermal component vs time (in observer frame)
	plot_param_vs_time(thermal_emission,'TEMP', ax=ax[0], z=z, y_factor=kb_kev/(1+z),Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,frame=frame)
	
	ax[0].set_xlabel(r't$_{obs}$',fontsize=fontsize,fontweight=fontweight)
	ax[0].set_ylabel(r'k$_B$T (KeV)',fontsize=fontsize,fontweight=fontweight)
	ax[0].set_yscale('log')

	# Plot Rphot vs Tphot
	# ax[1].scatter(thermal_emission['RPHOT'],thermal_emission['TEMP'])
	plot_param_vs_time(thermal_emission,'RPHOT', y_factor=3*1e10, ax=ax[1], z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,frame=frame)

	ax[1].set_yscale('log')

	if(frame == "obs"):
		ax[1].set_xlabel(r't$_{obs}$ (sec)',fontsize=fontsize,fontweight=fontweight)
	if(frame == "source"):
		ax[1].set_xlabel(r't$_{e} (sec)$',fontsize=fontsize,fontweight=fontweight)
	ax[1].set_ylabel(r'R$_{phot}$ (cm)',fontsize=fontsize,fontweight=fontweight)
	# ax[1].set_ylabel(r'T$_{phot}}$ (K)',fontsize=fontsize,fontweight=fontweight)

	plot_aesthetics(ax[0],fontsize=fontsize,fontweight=fontweight)
	plot_aesthetics(ax[1],fontsize=fontsize,fontweight=fontweight)
	
	plt.tight_layout()
	if save_pref is not None:
		plt.savefig('figs/{}-thermal-evo.png'.format(save_pref))

##############################################################################################################################

def plot_evo_int_shock(is_emission,ax=None,z=0,Tmin=None, Tmax=None,save_pref=None,fontsize=14,fontweight='bold',frame="obs"):
	"""
	Plot evolution of internal shock emission parameters 
	"""

	if ax is None:
		fig, ax = plt.subplots(2,1,figsize=(7,8),sharex=True)
	
	# Make a copy of the axis in order to over plot two separate data sets
	ax0cp = ax[0].twinx()

	# Plot Arrival Time (ta) vs Emission Time (te)
	plot_param_vs_time(is_emission,'TE', ax=ax[0], z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker='^',frame=frame)
	# Plot Arrival Time (ta) vs delta T
	plot_param_vs_time(is_emission,'DELT', ax=ax0cp, z=z,Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker='.',color='C1',frame=frame)
	
	ax[0].set_ylabel(r'$t_{e}$',fontsize=fontsize,fontweight=fontweight)
	ax0cp.set_ylabel(r'$\Delta t$',fontsize=fontsize,fontweight=fontweight)
	ax0cp.yaxis.set_label_position("right")
	ax0cp.yaxis.tick_right()

	# Make a copy of the axis in order to over plot two separate data sets
	ax1cp = ax[1].twinx()
	
	# Plot Arrival Time (ta) vs the dissipated energy (e)
	plot_param_vs_time(is_emission,'EDISS', ax=ax[1], y_factor=is_emission['DELT'], z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker='^',frame=frame)
	# Plot Arrival Time (ta) vs approximate Lorentz factor (gamma_r)
	plot_param_vs_time(is_emission,'GAMMAR', ax=ax1cp, y_factor=1/100, z=z,Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker='.',color='C1',frame=frame)

	ax[1].set_ylabel(r'E$_{diss}/\Delta$t$_e$',fontsize=fontsize,fontweight=fontweight)
	if(frame == "obs"):
		ax[1].set_xlabel(r't$_{obs}$ (sec)',fontsize=fontsize,fontweight=fontweight)
	if(frame == "source"):
		ax[1].set_xlabel(r't$_{e} (sec)$',fontsize=fontsize,fontweight=fontweight)
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
	plot_param_vs_time(is_emission,'ASYN', ax=ax[0,0], z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,frame=frame)
	ax[0,0].set_ylabel(r'$\alpha_{syn}$',fontsize=fontsize,fontweight=fontweight)
	
	# Plot Arrival Time (ta) vs the dissipated energy (e)
	plot_param_vs_time(is_emission,'GAMMAE', ax=ax[0,1], y_factor=1/1e4,z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,frame=frame)
	ax[0,1].set_ylabel(r'$\Gamma_{e}$/1e4',fontsize=fontsize,fontweight=fontweight)
	ax[0,1].yaxis.set_label_position("right")
	ax[0,1].yaxis.tick_right()

	# Plot Arrival Time (ta) vs the equipartition magnetic field (Beq)
	plot_param_vs_time(is_emission,'BEQ', ax=ax[1,0], z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker='.',frame=frame)
	ax[1,0].set_yscale('log')
	ax[1,0].set_ylabel(r'B$_{eq}$',fontsize=fontsize,fontweight=fontweight)
	if(frame == "obs"):
		ax[1,0].set_xlabel(r't$_{obs}$ (sec)',fontsize=fontsize,fontweight=fontweight)
	if(frame == "source"):
		ax[1,0].set_xlabel(r't$_{e} (sec)$',fontsize=fontsize,fontweight=fontweight)

	# Plot Arrival Time (ta) vs the synchrotron energy (Esyn)
	plot_param_vs_time(is_emission,'ESYN', ax=ax[1,1], z=z, Tmin=Tmin, Tmax=Tmax,
		fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,frame=frame)
	ax[1,1].set_yscale('log')
	ax[1,1].set_ylabel(r'$E_{syn}$/1e3',fontsize=fontsize,fontweight=fontweight)
	if(frame == "obs"):
		ax[1,1].set_xlabel(r't$_{obs}$ (sec)',fontsize=fontsize,fontweight=fontweight)
	if(frame == "source"):
		ax[1,1].set_xlabel(r't$_{e} (sec)$',fontsize=fontsize,fontweight=fontweight)
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

	if save_pref is not None:
		plt.savefig('figs/{}-int-shock-evo-fig1.png'.format(save_pref))

##############################################################################################################################

def make_together_plots(shock_data, ax0, ax1, label=None, color="C1",marker=".", z=0, Tmin=None, Tmax=None,fontsize=14,fontweight='bold',guidelines=False,save_pref=None,frame="obs",markersize=7):


	### First Plot ###

	# T_a vs B_eq
	plot_param_vs_time(shock_data,'BEQ', ax=ax0[0,0], z=z, Tmin=Tmin, Tmax=Tmax,
			fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker=marker,color=color,label=label,frame=frame,markersize=markersize)

	# T_a vs Gamma_r
	if guidelines == True:
		rhowindline = lambda t, t0, norm: norm*np.power(t/t0,-1./4.)
		rhoconstline = lambda t, t0, norm: norm*np.power(t/t0,-3./8.)
		tstart = 1e2
		tstop = 1e8
		tnum = 1e2
		g_norm_wind = rhowindline(tstart,shock_data['TA'][np.argmax(shock_data['TA']>tstart)],shock_data['GAMMAR'][np.argmax(shock_data['TA']>tstart)])
		g_norm_const = rhoconstline(tstart,shock_data['TA'][np.argmax(shock_data['TA']>tstart)],shock_data['GAMMAR'][np.argmax(shock_data['TA']>tstart)])
		t = np.linspace(tstart,tstop,num=int(tnum))	
		ax0[0,1].plot(t,rhowindline(t,tstart,g_norm_wind),label=r"$t^{-1/4}$",color='r')
		ax0[0,1].plot(t,rhoconstline(t,tstart,g_norm_const),label=r"$t^{-3/8}$",color='k')

	plot_param_vs_time(shock_data,'GAMMAR', ax=ax0[0,1], z=z, Tmin=Tmin, Tmax=Tmax,
			fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker=marker,color=color,frame=frame,markersize=markersize)

	# T_a vs e_diss
	plot_param_vs_time(shock_data,'EDISS', ax=ax0[1,0], z=z, Tmin=Tmin, Tmax=Tmax,
			fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker=marker,color=color,frame=frame,markersize=markersize)

	# T_a vs E_syn
	plot_param_vs_time(shock_data,'ESYN', ax=ax0[1,1], z=z, Tmin=Tmin, Tmax=Tmax,
			fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker=marker,color=color,frame=frame,markersize=markersize)


	# Plot Aesthetics

	# Format Top Left plot, T_a vs B_eq
	ax0[0,0].set_ylabel(r"B$_{EQ}$",fontsize=fontsize,fontweight=fontweight)
	# ax0[0,0].set_xlabel(r"$t_a$",fontsize=fontsize,fontweight=fontweight)
	ax0[0,0].set_yscale("log")
	ax0[0,0].set_xscale("log")
	ax0[0,0].legend(fontsize=fontsize)

	# Format Top Right plot, T_a vs Gamma_r 
	ax0[0,1].set_ylabel(r"$\Gamma_r$",fontsize=fontsize,fontweight=fontweight)
	# ax0[0,1].set_xlabel(r"$t_a$",fontsize=fontsize,fontweight=fontweight)
	ax0[0,1].yaxis.set_label_position("right")
	ax0[0,1].yaxis.tick_right()
	ax0[0,1].set_yscale("log")
	ax0[0,1].set_xscale("log")

	# Format Bottom Left plot, T_a vs E_diss
	ax0[1,0].set_ylabel(r"E$_{diss}/\Delta$t$_{e}$",fontsize=fontsize,fontweight=fontweight)
	if(frame == "obs"):
		ax0[1,0].set_xlabel(r't$_{obs}$ (sec)',fontsize=fontsize,fontweight=fontweight)
	if(frame == "source"):
		ax0[1,0].set_xlabel(r't$_{e} (sec)$',fontsize=fontsize,fontweight=fontweight)
	ax0[1,0].set_yscale("log")
	ax0[1,0].set_xscale("log")

	# Format Bottom Right plot, T_a vs E_synch
	ax0[1,1].set_ylabel(r"E$_{syn}$",fontsize=fontsize,fontweight=fontweight)
	if(frame == "obs"):
		ax0[1,1].set_xlabel(r't$_{obs}$ (sec)',fontsize=fontsize,fontweight=fontweight)
	if(frame == "source"):
		ax0[1,1].set_xlabel(r't$_{e} (sec)$',fontsize=fontsize,fontweight=fontweight)
	ax0[1,1].yaxis.set_label_position("right")
	ax0[1,1].yaxis.tick_right()
	ax0[1,1].set_yscale("log")
	ax0[1,1].set_xscale("log")

	# Make plots look good
	for i in range(2):
		for j in range(2):
			plot_aesthetics(ax0[i,j],fontsize=fontsize,fontweight=fontweight)

	plt.tight_layout()
	plt.subplots_adjust(wspace=0,hspace=0)

	if save_pref is not None:
		plt.savefig('figs/{}-all-shock-evo-fig0.png'.format(save_pref))

	### Second Plot ###

	# T_a vs T_e
	plot_param_vs_time(shock_data,'TE', ax=ax1[0], z=z, Tmin=Tmin, Tmax=Tmax,
			fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False, marker=marker,color=color, label=label,frame=frame,markersize=markersize)

	# T_a vs del T_a
	# Make a copy of the axis in order to over plot two separate data sets
	# ax1cp = ax1[0].twinx()
	# plot_param_vs_time(shock_data,'DELT', ax=ax1cp, z=z,Tmin=Tmin, Tmax=Tmax,
	# 	fontsize=fontsize, fontweight=fontweight, disp_xax=False, disp_yax=False,marker=marker,color=color,frame=frame)


	# T_a vs Gamma_e
	plot_param_vs_time(shock_data,'GAMMAE', ax=ax1[1], z=z, Tmin=Tmin, Tmax=Tmax,
			fontsize=fontsize, fontweight=fontweight, marker=marker,color=color,frame=frame,markersize=markersize)

	# Plot Aesthetics
	# Format Top plot, T_a vs T_e
	
	ax1[0].set_ylabel(r"$t_e$",fontsize=fontsize,fontweight=fontweight)
	ax1[0].set_xlabel(r"$t_a$",fontsize=fontsize,fontweight=fontweight)

	# ax1cp.set_ylabel(r'$\Delta t$',fontsize=fontsize,fontweight=fontweight)
	# ax1cp.yaxis.set_label_position("right")
	# ax1cp.yaxis.tick_right()

	ax1[0].set_yscale("log")
	ax1[0].set_xscale("log")
	# ax1cp.set_yscale("log")
	# ax1cp.set_xscale("log")
	ax1[0].legend(fontsize=fontsize)

	# Format Bottom plot, T_a vs Gamma_e
	ax1[1].set_ylabel(r"$\Gamma_{e}$",fontsize=fontsize,fontweight=fontweight)
	if(frame == "obs"):
		ax1[1].set_xlabel(r't$_{obs}$ (sec)',fontsize=fontsize,fontweight=fontweight)
	if(frame == "source"):
		ax1[1].set_xlabel(r't$_{e}$ (sec)',fontsize=fontsize,fontweight=fontweight)
	ax1[1].set_yscale("log")
	ax1[1].set_xscale("log")

	# Make plots look good
	for i in range(2):
		plot_aesthetics(ax1[i],fontsize=fontsize,fontweight=fontweight)
	# plot_aesthetics(ax1cp,fontsize=fontsize,fontweight=fontweight)


	plt.tight_layout()
	plt.subplots_adjust(hspace=0)

	if save_pref is not None:
		plt.savefig('figs/{}-all-shock-evo-fig1.png'.format(save_pref))

##############################################################################################################################

def plot_together(is_data = None,fs_data=None, rs_data=None, z=0, Tmin=None, Tmax=None,save_pref=None,fontsize=14,fontweight='bold',frame="obs",markregime=True,markersize=10):

	fig0, ax0 = plt.subplots(2,2,sharex=True,figsize=(12,8))
	fig1, ax1 = plt.subplots(2,1,sharex=True,figsize=(6,6))

	annot0 = []
	annot0.append( ax0[0,0].annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points", bbox=dict(boxstyle="round", fc="w"),arrowprops=dict(arrowstyle="->")) )
	annot0.append( ax0[0,1].annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points", bbox=dict(boxstyle="round", fc="w"),arrowprops=dict(arrowstyle="->")) )
	annot0.append( ax0[1,0].annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points", bbox=dict(boxstyle="round", fc="w"),arrowprops=dict(arrowstyle="->")) )
	annot0.append( ax0[1,1].annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points", bbox=dict(boxstyle="round", fc="w"),arrowprops=dict(arrowstyle="->")) )
	
	annot1 = []
	annot1.append( ax1[0].annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points", bbox=dict(boxstyle="round", fc="w"),arrowprops=dict(arrowstyle="->")) )
	annot1.append( ax1[1].annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points", bbox=dict(boxstyle="round", fc="w"),arrowprops=dict(arrowstyle="->")) )

	data_holder = []

	def onpick3(event, fig, ax, data_holder):

		artist = event.artist
		for i in range( len(ax.get_children() ) ):
			if(ax.get_children()[i] == artist ):
				
				facecolor = ax.get_children()[i].get_facecolor()

				# This is a list of indices for each point located underneath of the mouse click
				ind_list = event.ind
				# Take central one in the list  
				ind = ind_list[int(len(ind_list)/2)]

				# For fig 0
				for j in range(2):
					for k in range(2):
						# For fig 1submit_time
						# Access the selected data point position 
						d = ax0[j,k].collections[i]
						pos = d.get_offsets()[ind]

						# Set the position of the annotation 
						annot0[k+(j*2)].xy = pos

						# Display the point index
						text = "SH. IND. = {}".format(data_holder[i]["SHIND"][ind])
						annot0[k+(j*2)].set_text(text)

						# Set color and alpha of text box for clarity
						annot0[k+(j*2)].get_bbox_patch().set_facecolor(facecolor)
						annot0[k+(j*2)].get_bbox_patch().set_alpha(0.4)

						# Redraw all axes to update text box position and information
						ax0[j,k].redraw_in_frame()

				# For fig 1
				for j in range(2):
					# For fig 1
					# Access the selected data point position 
					d = ax1[j].collections[i]
					pos = d.get_offsets()[ind]

					# Set the position of the annotation 
					annot1[j].xy = pos

					# Display the point index
					text = "SH. IND. = {}".format(data_holder[i]["SHIND"][ind])
					annot1[j].set_text(text)

					# Set color and alpha of text box for clarity
					annot1[j].get_bbox_patch().set_facecolor(facecolor)
					annot1[j].get_bbox_patch().set_alpha(0.4)

					# Redraw all axes to update text box position and information
					ax1[j].redraw_in_frame()

			# Redraw all figures to update text box position and information
			fig0.canvas.draw_idle()
			fig1.canvas.draw_idle()



	if is_data is not None:
		if (markregime == True):
			fastcool_data = is_data[is_data['NUM']>is_data['NUC']]
			make_together_plots(shock_data=fastcool_data,label="IS", color="C0", ax0=ax0, ax1=ax1, z=z, Tmin=Tmin, Tmax=Tmax, fontsize=fontsize,fontweight=fontweight,frame=frame,markersize=markersize)
			data_holder.append(fastcool_data)

			slowcool_data = is_data[is_data['NUM']<=is_data['NUC']]
			make_together_plots(shock_data=slowcool_data, color="C0", marker='x', ax0=ax0, ax1=ax1, z=z, Tmin=Tmin, Tmax=Tmax, fontsize=fontsize,fontweight=fontweight,frame=frame,markersize=markersize)
			data_holder.append(slowcool_data)

		else:
			make_together_plots(shock_data=is_data,label="IS", color="C0", ax0=ax0, ax1=ax1, z=z, Tmin=Tmin, Tmax=Tmax, fontsize=fontsize,fontweight=fontweight,frame=frame,markersize=markersize)
			data_holder.append(is_data)

	if fs_data is not None:
		if (markregime == True):
			fastcool_data = fs_data[fs_data['NUM']>fs_data['NUC']]
			make_together_plots(shock_data=fastcool_data,label="FS", color="C1", ax0=ax0, ax1=ax1, z=z, Tmin=Tmin, Tmax=Tmax, fontsize=fontsize,fontweight=fontweight,frame=frame,markersize=markersize)
			data_holder.append(fastcool_data)
			
			slowcool_data = fs_data[fs_data['NUM']<=fs_data['NUC']]
			make_together_plots(shock_data=slowcool_data, color="C1", marker='x', ax0=ax0, ax1=ax1, z=z, Tmin=Tmin, Tmax=Tmax, fontsize=fontsize,fontweight=fontweight,frame=frame,markersize=markersize)
			data_holder.append(slowcool_data)

		else:
			make_together_plots(shock_data=fs_data,label="FS", color="C1", ax0=ax0, ax1=ax1, z=z, Tmin=Tmin, Tmax=Tmax, fontsize=fontsize,fontweight=fontweight,frame=frame,markersize=markersize)
			data_holder.append(fs_data)

	if rs_data is not None:
		if (markregime == True):
			fastcool_data = rs_data[rs_data['NUM']>rs_data['NUC']]
			make_together_plots(shock_data=fastcool_data,label="RS", color="C2", ax0=ax0, ax1=ax1, z=z, Tmin=Tmin, Tmax=Tmax, fontsize=fontsize,fontweight=fontweight,frame=frame,markersize=markersize)
			data_holder.append(fastcool_data)

			slowcool_data = rs_data[rs_data['NUM']<=rs_data['NUC']]
			make_together_plots(shock_data=slowcool_data, color="C2", marker='x', ax0=ax0, ax1=ax1, z=z, Tmin=Tmin, Tmax=Tmax, fontsize=fontsize,fontweight=fontweight,frame=frame,markersize=markersize)
			data_holder.append(slowcool_data)

		else:
			make_together_plots(shock_data=rs_data,label="RS", color="C2", ax0=ax0, ax1=ax1, z=z, Tmin=Tmin, Tmax=Tmax, fontsize=fontsize,fontweight=fontweight,frame=frame,markersize=markersize)
			data_holder.append(rs_data)


	fig0.canvas.mpl_connect("pick_event", lambda event: onpick3(event, fig0, ax0[0,0],data_holder))
	fig0.canvas.mpl_connect("pick_event", lambda event: onpick3(event, fig0, ax0[0,1],data_holder))
	fig0.canvas.mpl_connect("pick_event", lambda event: onpick3(event, fig0, ax0[1,0],data_holder))
	fig0.canvas.mpl_connect("pick_event", lambda event: onpick3(event, fig0, ax0[1,1],data_holder))

	fig1.canvas.mpl_connect('pick_event', lambda event: onpick3(event, fig1, ax1[0],data_holder))
	fig1.canvas.mpl_connect('pick_event', lambda event: onpick3(event, fig1, ax1[1],data_holder))


	if save_pref is not None:
		fig0.savefig('figs/{}-all-shock-evo-fig0.png'.format(save_pref))
		fig1.savefig('figs/{}-all-shock-evo-fig1.png'.format(save_pref))

	return fig0, fig1

##############################################################################################################################

def plot_observables(is_emission,th_emission,ax=None,z=0, save_pref=None,fontsize=14,fontweight='bold',frame="obs"):
	"""
	Plot observables
	"""

	if ax is None:
		fig, ax = plt.subplots(2,2,figsize=(8,8))

	# Epeak vs Time of the non-thermal component
	plot_param_vs_time(is_emission, "ESYN", ax=ax[0,0], z=z, fontsize=fontsize, fontweight=fontweight, frame=frame)

	if(frame == "obs"):
		ax[0,0].set_xlabel(r't$_{obs}$ (sec)',fontsize=fontsize,fontweight=fontweight)
	if(frame == "source"):
		ax[0,0].set_xlabel(r't$_{e}$ (sec)',fontsize=fontsize,fontweight=fontweight)
	ax[0,0].set_ylabel(r"E$_{peak}$ (keV)",fontsize=fontsize,fontweight=fontweight)

	# Flux vs Epeak of the non-thermal component
	ax[0,1].scatter(x=is_emission['ESYN'],y=is_emission['EDISS'],marker=".")

	ax[0,1].set_yscale('log')
	ax[0,1].set_xscale('log')

	ax[0,1].set_xlabel(r"E$_{peak}$ (keV)",fontsize=fontsize,fontweight=fontweight)
	ax[0,1].set_ylabel(r'F (erg/s)',fontsize=fontsize,fontweight=fontweight)
	
	# alpha (low energy power law index) vs Time of the non-thermal component



	# Flux vs Temperature of the thermal component
	ax[1,1].scatter(x=th_emission['TEMP']*kb_kev,y=th_emission['FLUX'],marker=".")

	ax[1,1].set_xlabel(r"k$_B$T (keV)",fontsize=fontsize,fontweight=fontweight)
	ax[1,1].set_ylabel(r'F (erg/s)',fontsize=fontsize,fontweight=fontweight)

	# Make plots look good
	for i in range(2):
		for j in range(2):
			plot_aesthetics(ax[i,j],fontsize=fontsize,fontweight=fontweight)
			ax[i,j].grid()

	plt.tight_layout()

	if save_pref is not None:
		plt.savefig('figs/{}-observables.png'.format(save_pref))

##############################################################################################################################

def plot_synch_cooling_regime(emission,ax=None,z=0,label=None, color="C0", markers=[".","^"], markersize=None, alpha=1, save_pref=None, Tmin=None, Tmax=None, fontsize=14,fontweight='bold',frame="obs"):
	"""
	Plot nu_c and nu_m vs time
	"""

	# Make plot instance if it doesn't exist
	if ax is None:
		ax = plt.figure(figsize=(10,8)).gca()

	# Plot temperature of the thermal component vs time (in observer frame)
	plot_param_vs_time(emission,'NUC', ax=ax, z=z, y_factor=planck_kev*emission['GAMMAR'],Tmin=Tmin, Tmax=Tmax, marker=markers[0],
		markersize=markersize, color=color, label=label, fontsize=fontsize, fontweight=fontweight,frame=frame,alpha=alpha)
	plot_param_vs_time(emission,'NUM', ax=ax, z=z, y_factor=planck_kev*emission['GAMMAR'],Tmin=Tmin, Tmax=Tmax, marker=markers[1],
		markersize=markersize, color=color, fontsize=fontsize, fontweight=fontweight,frame=frame,alpha=alpha)

	ax.set_xlabel(r't$_{obs}$',fontsize=fontsize,fontweight=fontweight)
	ax.set_ylabel(r'E (KeV)',fontsize=fontsize,fontweight=fontweight)
	ax.set_title(r"$\nu_c$ = {}, $\nu_m$ = {}".format(repr(markers[0]),repr(markers[1])),fontsize=fontsize,fontweight=fontweight)
	ax.set_yscale('log')


	# Make plots look good
	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)

	ax.set_yscale('log')

	if label is not None:
		ax.legend(fontsize=fontsize)
	

	ax.grid(True)
	plt.tight_layout()

	if save_pref is not None:
		plt.savefig('figs/{}-synch-cooling-regime.png'.format(save_pref))


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

	z = 0

	"""
	Shell Lorentz Distribution
	"""
	
	# plot_lor_dist('data-file-dir/synthGRB_shell_dist.txt')
	# ani = plot_lor_dist_anim('data-file-dir/synthGRB_shell_dist.txt')

	"""
	Synthetic spectrum 
	"""
	"""
	ax_spec = plt.figure(figsize=(9,8)).gca()

	## Synthetic spectra with each component
	plot_spec("data-file-dir/synthGRB_spec_therm.txt",ax=ax_spec,z=z,label="Thermal",color="r")
	plot_spec("data-file-dir/synthGRB_spec_is.txt",ax=ax_spec,z=z,label="IS",color="C0")
	plot_spec("data-file-dir/synthGRB_spec_fs.txt",ax=ax_spec,z=z,label="FS",color="C1")
	plot_spec("data-file-dir/synthGRB_spec_rs.txt",ax=ax_spec,z=z,label="RS",color="C2")
	plot_spec("data-file-dir/synthGRB_spec_total.txt",ax=ax_spec,z=z,label="Total",color="k")

	## Synthetic spectrum before convolusion
	# plot_spec("data-file-dir/spec_source.txt",ax=ax_spec,unc=False,label="Source")
	# plot_spec("data-file-dir/spec_source_fluc.txt",ax=ax_spec,unc=True,label="Pre-Conv")
	# plot_spec("data-file-dir/spec_model.txt",ax=ax_spec,unc=False,label="Model",joined=True)
	
	## Synthetic spectrum after convolusion
	# plot_spec("data-file-dir/spec_obs.txt",ax=ax_spec,unc=True,label="Obs")
	# plot_spec("data-file-dir/spec_model_conv.txt",ax=ax_spec,unc=False,label="Model",joined=True)
	# plot_spec("data-file-dir/spec_mod_emp.txt",ax=ax_spec,unc=False,label="Model",joined=True)

	add_FermiGBM_band(ax_spec)
	# add_SwiftBAT_band(ax_spec)

	# ax_spec.set_xlim(0.1,1e5)
	# ax_spec.set_ylim(1e48,1e52)
	"""
	

	"""
	Synthetic light curve
	"""
	"""
	ax_lc = plt.figure().gca()
	plot_light_curve("data-file-dir/synthGRB_light_curve_TH.txt",ax=ax_lc,z=z,label="TH",color="r")
	plot_light_curve("data-file-dir/synthGRB_light_curve_IS.txt",ax=ax_lc,z=z,label="IS",color="C0")
	plot_light_curve("data-file-dir/synthGRB_light_curve_FS.txt",ax=ax_lc,z=z,label="FS",color="C1")
	plot_light_curve("data-file-dir/synthGRB_light_curve_RS.txt",ax=ax_lc,z=z,label="RS",color="C2")
	plot_light_curve("data-file-dir/synthGRB_light_curve.txt",ax=ax_lc,z=z,label="Total",logscale=False,color="k")

	# Interactive light curve
	# tbox = plot_light_curve_interactive(init_Tmin = 0, init_Tmax = 20, init_Emin = 8, init_Emax = 1e4,z=z,label="Total",with_comps=True)
	

	# Afterglow light curve
	ax_afg_lc = plt.figure().gca()
	plot_light_curve("data-file-dir/synthGRB_light_curve.txt",ax=ax_afg_lc,z=z,label="Prompt",logscale=True,color="k")
	plot_light_curve("data-file-dir/synthGRB_light_curve_afterglow_gbm.txt",ax=ax_afg_lc,z=z,label="AG: GBM",logscale=True,color="C1")
	# plot_light_curve("data-file-dir/synthGRB_light_curve_afterglow_xrt.txt",ax=ax_afg_lc,z=z,label="AG: XRT",logscale=True,color="C4")
	# plot_light_curve("data-file-dir/synthGRB_light_curve_afterglow_opt.txt",ax=ax_afg_lc,z=z,label="AG: OPT",logscale=True,color="C5")
	ax_afg_lc.set_ylim(1e43,1e49)
	ax_afg_lc.set_xlim(0.1)
	"""
	"""
	Jet dynamics plots 
	
	"""
	"""
	# therm_emission = load_therm_emission("data-file-dir/synthGRB_jet_params_therm.txt")
	# plot_evo_therm(therm_emission,z=1)
	
	is_data = load_is_emission("data-file-dir/synthGRB_jet_params_is.txt")
	# plot_evo_int_shock(is_data,z=0.5)

	fs_data = load_fs_emission("data-file-dir/synthGRB_jet_params_fs.txt")
	rs_data = load_rs_emission("data-file-dir/synthGRB_jet_params_rs.txt")

	# plot_evo_ext_shock(fs_data=fs_data,rs_data=rs_data,z=0)
	# plot_evo_ext_shock(fs_data=fs_data)		
	
	# Plot everything together:
	fig0, fig1 = plot_together(is_data=is_data,fs_data=fs_data,rs_data=rs_data)
	# fig0, fig1 = plot_together(fs_data=fs_data)

	
	# Plot nu_c and nu_m: 
	# ax_synch_reg = plt.figure(figsize=(10,8)).gca()
	# fig, ax_synch_reg = plt.subplots()
	# markers = [".","^"]
	# plot_synch_cooling_regime(is_data,ax=ax_synch_reg,Tmin=0,Tmax=20,label="IS",color="C0",markers=markers,alpha=0.8,markersize=16)
	# # plot_synch_cooling_regime(fs_data,ax=ax_synch_reg,label="FS",color="C1")
	# plot_synch_cooling_regime(rs_data,ax=ax_synch_reg,Tmin=0,Tmax=20,label="RS",color="C2",markers=markers,alpha=0.8,markersize=16)
	# add_FermiGBM_band(ax_synch_reg,axis="y")
	"""
	
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
	Observables
	"""
	"""
	th_emission = load_therm_emission("data-file-dir/synthGRB_jet_params_therm.txt")
	is_emission = load_is_emission("data-file-dir/synthGRB_jet_params_is.txt")

	plot_observables(is_emission, th_emission)
	"""


	"""
	Testing
	"""



	plt.show()

