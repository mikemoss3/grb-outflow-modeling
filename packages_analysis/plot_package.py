"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2023-06-13


Package that contains all plotting functions for plotting synthetic and observed data.

"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from asymmetric_uncertainty import neg_errors, pos_errors
# import subprocess

# from matplotlib.widgets import TextBox

# Import custom classes
from packages_analysis.data_package import *
from packages_analysis.model_package import *
from packages_analysis.fit_package import *

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

##############################################################################################################################

def _load_lor_dist(file_name, string_match = "// Next step\n"):
	"""
	Method to load Lorentz factor distribution from the file specified by "file_name"

	The Lorentz distribution file must contain the columns: 
	RADIUS - Radius of the shell
	GAMMA - Lorentz factor of the shell
	MASS - Mass of the shell
	TE - Time of emission of the shell
	STATUS - Status of the shell, this is used by the simulation code to indicate if a shell is still active or not. 
	
	Attributes:
	separator_string = the string between each snapshot of the Lorentz distribution
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
	data_time = []
	data_list = []
	for i in range(len(line_num_list)-1):
		
		# Read the time at which this data occurred in the simulation
		tmp_time = np.genfromtxt(file_name,skip_header=line_num_list[i],max_rows=1)
		data_time.append(tmp_time) 
		
		# Read data
		tmp_data = np.genfromtxt(file_name,skip_header=line_num_list[i]+1,max_rows=line_num_list[i+1]-2-line_num_list[i], dtype=[("RADIUS",float),("GAMMA",float),("MASS",float),("TE",float),("STATUS",float)])
		data_list.append(tmp_data[tmp_data['STATUS']>0])

	return data_time, data_list

##############################################################################################################################

def plot_lor_prof(file_name,ax=None,color="C0",save_pref=None,xlabel=True,ylabel=True,label=None,fontsize=14,fontweight='bold',marker='.', separator_string = "// Next step\n",joined=True,title=True,zoom_inset=False,ind_start=0,alpha=0.7,linestyle="solid",linewidth=1.5):
	"""
	Method to plot the given Lorentz factor distribution saved in the text file with path name "file_name".
	Multiple snapshots of the Lorentz distribution can be given in a single file. Each Lorentz distribution must be separated by a line with the string indicated by "separator_string".
	If more than one snapshot is provided, the snapshots can be scrolled through with the left and right arrow keys. 

	The Lorentz distribution file must contain the columns: 
	RADIUS - Radius of the shell
	GAMMA - Lorentz factor of the shell
	MASS - Mass of the shell
	TE - Time of emission of the shell
	STATUS - Status of the shell, this is used by the simulation code to indicate if a shell is still active or not. 

	Attributes:
	ax = 			the matplotlib.pyplot.axes instance to make the plot on	
	save_pref = 	if not left as None, the plot will be saved and the file name will have this prefix
	xlabel = 		boolean, indicates whether x-labels should be displayed
	ylabel = 		boolean, indicates whether y-labels should be displayed
	label = 		optional label of plotted line
	fontsize = 		size of the plot and label font
	fontweight = 	weight of the plot and label font
	alpha = 		opacity of plot lines
	linestyle = 	style of the plot lines
	linewidth = 	width of plot lines
	separator_string = the string between each snapshot of the Lorentz distribution
	"""

	# Load data
	lor_time_list, lor_dist_list = _load_lor_dist(file_name,separator_string)

	# Make plot instance if it doesn't exist
	if ax is None:
		ax = plt.figure().gca()
	fig = plt.gcf()

	ind = np.array([ind_start])
	def on_press(event):
		# if event.key == 'right':
		if event.key == 'down':
			ind[0] +=1
		# elif event.key == 'left':
		elif event.key == 'up':
			ind[0] -=1

		# Loop back to beginning 
		if(ind[0] == len(lor_time_list)):
			ind[0] = 0
		# Loop to the end
		if(ind[0] == -1):
			ind[0] = len(lor_time_list)-1

		# Update ejection time plot 
		if(joined == True):
			line_shell_ind.set_xdata(lor_dist_list[ind[0]]['TE'])
			line_shell_ind.set_ydata(lor_dist_list[ind[0]]['GAMMA'])
		if(joined == False):
			line_shell_ind.set_offsets(np.transpose([lor_dist_list[ind[0]]['TE'],lor_dist_list[ind[0]]['GAMMA']]))
		# ax[0].set_xlim(np.min(lor_dist_list[ind[0]]['TE']),np.max(lor_dist_list[ind[0]]['GAMMA']))

		if title == True:
			fig.suptitle("Emission Time = {0:.1e} sec".format(lor_time_list[ind[0]]),fontsize=fontsize,fontweight=fontweight)

		# Include zoom within an inset if specified 	
		if zoom_inset is True:
			if(joined == True):
				line_inset.set_xdata(lor_dist_list[ind[0]]['TE'])
				line_inset.set_ydata(lor_dist_list[ind[0]]['GAMMA'])
			elif(joined == False):
				line_inset.set_offsets(np.transpose([lor_dist_list[ind[0]]['TE'],lor_dist_list[ind[0]]['GAMMA']]))


		# Redraw the figure to implement updates
		ax.redraw_in_frame()
		fig.canvas.draw_idle()

	fig.canvas.mpl_connect('key_press_event', on_press)

	## Plot distribution as a function of the shell number 
	if(joined == True):
		line_shell_ind, = ax.step(lor_dist_list[ind_start]['TE'],lor_dist_list[ind_start]['GAMMA'],where="pre",color=color,alpha=alpha,linestyle=linestyle,label=label,linewidth=linewidth)
	elif(joined == False):
		line_shell_ind = ax.scatter(lor_dist_list[ind_start]['TE'],lor_dist_list[ind_start]['GAMMA'],marker=marker,color=color,alpha=alpha,label=label)

	# Include zoom within an inset if specified 	
	if zoom_inset is True:
		axins = ax.inset_axes([0.1, 0.5, 0.47, 0.47])

		if(joined == True):
			line_inset, = axins.step(lor_dist_list[0]['TE'],lor_dist_list[0]['GAMMA'],where="pre",color="k",alpha=0.5,linestyle="dashed",linewidth=linewidth)
		elif(joined == False):
			line_inset = axins.scatter(lor_dist_list[0]['TE'],lor_dist_list[0]['GAMMA'],marker=marker,color="k",alpha=0.5)

		if(joined == True):
			line_inset, = axins.step(lor_dist_list[ind_start]['TE'],lor_dist_list[ind_start]['GAMMA'],where="pre",color=color,alpha=alpha,linestyle=linestyle,linewidth=linewidth)
		elif(joined == False):
			line_inset = axins.scatter(lor_dist_list[ind_start]['TE'],lor_dist_list[ind_start]['GAMMA'],marker=marker,color=color,alpha=alpha)

		# sub region of the original image
		x1, x2, y1, y2 = 12, np.max(lor_dist_list[ind_start]['TE']), 0.1, 20
		axins.set_xlim(x1, x2)
		axins.set_ylim(y1, y2)

		axins.invert_xaxis()

		plot_aesthetics(axins,fontsize=fontsize-2,fontweight=fontweight)

	ax.set_ylim(0,np.max(lor_dist_list[0]['GAMMA']+10))
	ax.set_xlim(0,np.max(lor_dist_list[0]['TE']))
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	ax.invert_xaxis()

	if xlabel is True:
		ax.set_xlabel('Initial Ejection Time (sec)',fontsize=fontsize,fontweight=fontweight)
	if ylabel is True:
		ax.set_ylabel(r'$\Gamma$',fontsize=fontsize,fontweight=fontweight)


	if title == True:
		fig.suptitle("Emission Time = {0:.1e} sec".format(lor_time_list[ind_start]),fontsize=fontsize,fontweight=fontweight)


	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)

	if label is not None:
		ax.legend(fontsize=fontsize-2)

	plt.tight_layout()

	if save_pref is not None :
		plt.savefig('figs/{}-lorentz-profile.png'.format(save_pref))

##############################################################################################################################

def plot_lor_prof_evo(file_name,indices,ax=None,save_pref=None, xlabel=True,ylabel=True,label=None,fontsize=14,fontweight='bold', separator_string = "// Next step\n",title=True,zoom_inset=False,zoom_inset_range=[13,25,0,30],alpha=0.7,linestyle="solid",linewidth=1.5,color="C0",color_map=cm.PuBuGn,cm_norm_min=0.1,cm_norm_max=0.9,color_bar = True):

	if column == False:
		_lor_prof_single(file_name=file_name,indices=indices,ax=ax,save_pref=save_pref,xlabel=xlabel,ylabel=ylabel,label=label,fontsize=fontsize,fontweight=fontweight, separator_string = separator_string,title=title,zoom_inset=zoom_inset,zoom_inset_range=zoom_inset_range,alpha=alpha,linestyle=linestyle,linewidth=linewidth,color_map=color_map,cm_norm_min=cm_norm_min,cm_norm_max=cm_norm_max,color_bar = color_bar)
	else:
		_lor_prof_column(file_name=file_name,indices=indices,ax=ax,save_pref=save_pref,xlabel=xlabel,ylabel=ylabel,label=label,fontsize=fontsize,fontweight=fontweight, separator_string = separator_string,title=title,zoom_inset=zoom_inset,zoom_inset_range=zoom_inset_range,alpha=alpha,linestyle=linestyle,linewidth=linewidth,color=color)

##############################################################################################################################

def _lor_prof_single(file_name,indices,ax=None,save_pref=None,xlabel=True,ylabel=True,label=None,fontsize=14,fontweight='bold', separator_string = "// Next step\n",title=True,zoom_inset=False,zoom_inset_range=[13,25,0,30],alpha=0.7,linestyle="solid",linewidth=1.5,color_map=cm.PuBuGn,cm_norm_min=0.1,cm_norm_max=0.9,color_bar = True):
	"""
	Method to plot the given Lorentz factor distribution saved in the text file with path name "file_name".
	Multiple snapshots of the Lorentz distribution can be given in a single file. Each Lorentz distribution must be separated by a line with the string indicated by "separator_string".
	If more than one snapshot is provided, the snapshots will be overplotted in the same frame

	The Lorentz distribution file must contain the columns: 
	RADIUS - Radius of the shell
	GAMMA - Lorentz factor of the shell
	MASS - Mass of the shell
	TE - Time of emission of the shell
	STATUS - Status of the shell, this is used by the simulation code to indicate if a shell is still active or not. 

	Attributes:
	ax = 			the matplotlib.pyplot.axes instance to make the plot on	
	save_pref = 	if not left as None, the plot will be saved and the file name will have this prefix
	xlabel = 		boolean, indicates whether x-labels should be displayed
	ylabel = 		boolean, indicates whether y-labels should be displayed
	label = 		optional label of plotted line
	fontsize = 		size of the plot and label font
	fontweight = 	weight of the plot and label font
	alpha = 		opacity of plot lines
	linestyle = 	style of the plot lines
	linewidth = 	width of plot lines
	separator_string = the string between each snapshot of the Lorentz distribution
	"""

	# Load data
	lor_time_list, lor_dist_list = _load_lor_dist(file_name,separator_string)

	# Make plot instance if it doesn't exist
	if ax is None:
		ax = plt.figure().gca()

	colors = color_map( np.linspace(cm_norm_min,cm_norm_max,len(indices)) )

	for i in range(len(indices)):
		## Plot distribution as a function of the shell number 
		line_shell_ind, = ax.step(lor_dist_list[indices[i]]['TE'],lor_dist_list[indices[i]]['GAMMA'],where="pre",alpha=alpha,linestyle=linestyle,linewidth=linewidth,color=colors[i])


	# Include zoom within an inset if specified 	
	if zoom_inset is True:
		axins = ax.inset_axes([0.1, 0.5, 0.47, 0.47])

		for i in range(len(indices)):
			line_inset, = axins.step(lor_dist_list[indices[i]]['TE'],lor_dist_list[indices[i]]['GAMMA'],where="pre",alpha=alpha,linestyle=linestyle,linewidth=linewidth,color=colors[i])

		# sub region of the original image
		# x1, x2, y1, y2 = 12, np.max(lor_dist_list[indices[0]]['TE']), 5, 15
		x1, x2, y1, y2 = zoom_inset_range[0], zoom_inset_range[1], zoom_inset_range[2], zoom_inset_range[3]
		axins.set_xlim(x1, x2)
		axins.set_ylim(y1, y2)

		axins.invert_xaxis()

		plot_aesthetics(axins,fontsize=fontsize-2,fontweight=fontweight)


	if color_bar is True:
		fig = plt.gcf()
		fig.colorbar(cm.ScalarMappable(norm=plt.Normalize(vmin=cm_norm_min,vmax=cm_norm_max), cmap=color_map),location='top',orientation="horizontal",ticks=None)


	ax.set_ylim(0,np.max(lor_dist_list[indices[0]]['GAMMA']+10))
	ax.set_xlim(0,np.max(lor_dist_list[indices[0]]['TE']))
	ax.invert_xaxis()

	if xlabel is True:
		ax.set_xlabel('Initial Ejection Time (sec)',fontsize=fontsize,fontweight=fontweight)
	if ylabel is True:
		ax.set_ylabel(r'$\Gamma$',fontsize=fontsize,fontweight=fontweight)

	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)

	if label is not None:
		ax.legend(fontsize=fontsize)

	plt.tight_layout()

	if save_pref is not None :
		plt.savefig('figs/{}-lorentz-profile-simple.png'.format(save_pref))

##############################################################################################################################

def _lor_prof_column(file_name,indices,ax=None,save_pref=None,xlabel=True,ylabel=True,label=None,fontsize=14,fontweight='bold', separator_string = "// Next step\n",title=True,zoom_inset=False,zoom_inset_range=[13,25,0,30],alpha=0.7,linestyle="solid",linewidth=1.5,color="C0"):
	"""
	Method to plot the given Lorentz factor distribution saved in the text file with path name "file_name".
	Multiple snapshots of the Lorentz distribution can be given in a single file. Each Lorentz distribution must be separated by a line with the string indicated by "separator_string".
	If more than one snapshot is provided, the snapshots will be plotted in a column with each snapshot in a separate panel

	The Lorentz distribution file must contain the columns: 
	RADIUS - Radius of the shell
	GAMMA - Lorentz factor of the shell
	MASS - Mass of the shell
	TE - Time of emission of the shell
	STATUS - Status of the shell, this is used by the simulation code to indicate if a shell is still active or not. 

	Attributes:
	ax = 			the matplotlib.pyplot.axes instance to make the plot on	
	save_pref = 	if not left as None, the plot will be saved and the file name will have this prefix
	xlabel = 		boolean, indicates whether x-labels should be displayed
	ylabel = 		boolean, indicates whether y-labels should be displayed
	label = 		optional label of plotted line
	fontsize = 		size of the plot and label font
	fontweight = 	weight of the plot and label font
	alpha = 		opacity of plot lines
	linestyle = 	style of the plot lines
	linewidth = 	width of plot lines
	separator_string = the string between each snapshot of the Lorentz distribution
	"""


	# Load data
	lor_time_list, lor_dist_list = _load_lor_dist(file_name,separator_string)

	fig, ax = plt.subplots(len(indices),1,figsize=(6,int(len(indices*2))),sharex=True,sharey=True,gridspec_kw={'hspace':0.1})

	for i in range(len(indices)):
		ax[i].step(lor_dist_list[indices[i]]['TE'],lor_dist_list[indices[i]]['GAMMA'],where="pre",color=color,alpha=alpha,linestyle=linestyle,linewidth=linewidth)

		if zoom_inset is True:
			axins = ax[i].inset_axes([0.1, 0.5, 0.37, 0.37])

			line_inset, = axins.step(lor_dist_list[indices[i]]['TE'],lor_dist_list[indices[i]]['GAMMA'],where="pre",color=color,alpha=alpha,linestyle=linestyle,linewidth=linewidth)

			# sub region of the original image
			# x1, x2, y1, y2 = 12, np.max(lor_dist_list[indices[0]]['TE']), 5, 15
			x1, x2, y1, y2 = zoom_inset_range[0], zoom_inset_range[1], zoom_inset_range[2], zoom_inset_range[3]
			axins.set_xlim(x1, x2)
			axins.set_ylim(y1, y2)

			axins.invert_xaxis()

			plot_aesthetics(axins,fontsize=fontsize-2,fontweight=fontweight)


		ax[i].set_ylim(0,np.max(lor_dist_list[indices[0]]['GAMMA']+10))
		ax[i].set_xlim(0,np.max(lor_dist_list[indices[0]]['TE']))
		ax[i].invert_xaxis()

		if (xlabel is True) and ( i == len(indices)-1):
			ax[i].set_xlabel('Initial Ejection Time (sec)',fontsize=fontsize,fontweight=fontweight)
		if ylabel is True:
			ax[i].set_ylabel(r'$\Gamma$',fontsize=fontsize,fontweight=fontweight)

		plot_aesthetics(ax[i],fontsize=fontsize,fontweight=fontweight)

	if label is not None:
		ax[0].legend(fontsize=fontsize)

	# plt.tight_layout()

	if save_pref is not None :
		plt.savefig('figs/{}-lorentz-profile-column.png'.format(save_pref))

##############################################################################################################################

def add_FermiGBM_band(ax,fontsize=12,axis="x",plt_ratio_min=0.5,plt_ratio_max=1,inc_legend=True):
	"""
	Method to add two shaded boxes to indicate the Fermi/GBM (NaI and BGO) observation bands to a matplotlib.pyplot.axes instance

	Attributes:
	axis = Defines which axis the energy axis (in keV)
	"""

	# Grab the current ymin and ymax, this is used to set the lower and upper bounds of the vertical lines which indicate instrument observation energy range
	curr_ymin, curr_ymax = ax.get_ylim()
	curr_xmin, curr_xmax = ax.get_xlim()

	# Vertical axis
	if(axis == "x"):
		# Display Fermi/GBM - NAI energy band
		ax.axvspan(xmin=8,xmax=1e3,ymin=plt_ratio_min,ymax=plt_ratio_max,alpha=0.4,facecolor='grey',label='Fermi/GBM-NAI')

		# Display Fermi/GBM - BGO energy band
		ax.axvspan(xmin=150,xmax=3*1e4,ymin=plt_ratio_min,ymax=plt_ratio_max,alpha=0.4,facecolor='orange',label='Fermi/GBM-BGO')

	# Horizontal axis
	elif(axis == "y"):
		# Display Fermi/GBM - NAI energy band
		ax.axhspan(ymin=8,ymax=1e3,xmin=plt_ratio_min,xmax=plt_ratio_max,alpha=0.4,facecolor='grey',label='Fermi/GBM-NAI')

		# Display Fermi/GBM - BGO energy band
		ax.axhspan(ymin=150,ymax=3*1e4,xmin=plt_ratio_min,xmax=plt_ratio_max,alpha=0.4,facecolor='orange',label='Fermi/GBM-BGO')


	# Add to legend	
	if inc_legend is True:
		ax.legend(fontsize=fontsize)

	# We don't want the plotting window to change if either of the energy band edges do not overlap with the plotted energy spectra
	ax.set_ylim(curr_ymin,curr_ymax)
	ax.set_xlim(curr_xmin,curr_xmax)

##############################################################################################################################

def add_SwiftBAT_band(ax,fontsize=12,axis="x",plt_ratio_min=0.5,plt_ratio_max=1,inc_legend=True):
	"""
	Method to add two shaded boxes to indicate the Swift/BAT observation band to a matplotlib.pyplot.axes instance

	Attributes:
	axis = Defines which axis the energy axis (in keV)
	"""

	# Grab the current ymin and ymax, this is used to set the lower and upper bounds of the vertical lines which indicate instrument observation energy range
	curr_ymin, curr_ymax = ax.get_ylim()
	curr_xmin, curr_xmax = ax.get_xlim()

	# Display Swift/BAT energy band

	# Vertical axis
	if(axis == "x"):
		ax.axvspan(xmin=5,xmax=350,ymin=plt_ratio_min,ymax=plt_ratio_max,alpha=0.4,facecolor='blue',label='Swift/BAT')

	# Horizontal axis
	elif(axis == "y"):
		ax.axhspan(ymin=5,ymax=350,xmin=plt_ratio_min,xmax=plt_ratio_max,alpha=0.4,facecolor='blue',label='Swift/BAT')

	# Add to legend	
	if inc_legend is True:
		ax.legend(fontsize=fontsize)

	# We don't want the plotting window to change if either of the energy band edges do not overlap with the plotted energy spectra
	ax.set_ylim(curr_ymin,curr_ymax)
	ax.set_xlim(curr_xmin,curr_xmax)


##############################################################################################################################

def plot_model_spec(model, spec_type = 0,y_factor=1, inc_comps=True, emin=8, emax=4e4,ax=None,ymin=None,ymax=None,save_pref=None,xlabel=True,ylabel=True,label=None,fontsize=14,fontweight='bold', alpha=1,linestyle="solid",linewidth=1.5,comp_alpha=0.8,comp_linestyle="dashed",comp_linewidth=1.,en_steps = None, **kwargs):
	"""
	Method to plot the spectra stored in a model class 

	Attributes:
	model = 		model class instance to be plotted
	spec_type = 		spectrum type to be plotted
					0 == count spectrum
					1 == flux spectrum
					2 == power spectrum

	emin, emax = indicates the minimum and maximum energy range to plot. If None is supplied, the minimum and maximum energies of the supplied data files are used

	ax = 			the matplotlib.pyplot.axes instance to make the plot on	
	save_pref = 	if not left as None, the plot will be saved and the file name will have this prefix
	xlabel = 		boolean, indicates whether x-labels should be displayed
	ylabel = 		boolean, indicates whether y-labels should be displayed
	label = 		optional label of plotted line
	fontsize = 		size of the plot and label font
	fontweight = 	weight of the plot and label font
	alpha = 		opacity of plot lines
	linestyle = 	style of the plot lines
	linewidth = 	width of plot lines

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
			ax.plot(energy_axis,y_factor*model[i](energy_axis),alpha=comp_alpha,color=model[i].color,linestyle=comp_linestyle,linewidth=comp_linewidth,**kwargs)

	ax.plot(energy_axis,y_factor*model(energy_axis), alpha=alpha,color=model.color,linestyle=linestyle,linewidth=linewidth,label=label,**kwargs)

	# Plot aesthetics
	ax.set_xscale('log')
	ax.set_yscale('log')

	ylims = ax.get_ylim()
	# Force lower bound
	if np.max(model(energy_axis)*y_factor) < ylims[1]:
		ax.set_ylim(np.min(model(energy_axis)*y_factor)/10,ylims[1])
	else:
		ax.set_ylim(np.min(model(energy_axis)*y_factor)/10, np.max(model(energy_axis)*y_factor) * 10)
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

##############################################################################################################################

def plot_data_spec(spec_data, spec_type = 0,y_factor=1 , joined=False, unc=True, emin=None, emax=None,ax=None,ymin=None,ymax=None,save_pref=None,xlabel=True,ylabel=True,label=None,fontsize=14,fontweight='bold',alpha=1,linestyle="solid",linewidth=1.5,color="C0",marker=".", **kwargs):
	"""
	Method to plot synthetic or observed spectra 

	spec_type
		0 == count spectrum
		1 == flux spectrum
		2 == power spectrum

	ax = 			the matplotlib.pyplot.axes instance to make the plot on	
	save_pref = 	if not left as None, the plot will be saved and the file name will have this prefix
	xlabel = 		boolean, indicates whether x-labels should be displayed
	ylabel = 		boolean, indicates whether y-labels should be displayed
	label = 		optional label of plotted line
	fontsize = 		size of the plot and label font
	fontweight = 	weight of the plot and label font
	alpha = 		opacity of plot lines
	linestyle = 	style of the plot lines
	linewidth = 	width of plot lines
	"""

	# Make plot instance if it doesn't exist
	if ax is None:
		ax = plt.figure().gca()
	# fig = plt.gcf()

	if spec_type == 1:
		y_factor *= spec_data['ENERGY']
	if spec_type == 2:
		y_factor *= spec_data['ENERGY']**2

	# ax.set_ylim(1e47)
	# ylims = ax.get_ylim()

	# Plot spectrum data
	if joined is True:
		if unc is True:
			line = ax.errorbar(x=spec_data['ENERGY'],y=spec_data['RATE']*y_factor,yerr=(neg_errors(spec_data['RATE']),pos_errors(spec_data['RATE']))*y_factor,label=label,color=color,alpha=alpha,**kwargs)
		else:
			line, = ax.plot(spec_data['ENERGY'],spec_data['RATE']*y_factor,label=label,color=color,linestyle=linestyle,alpha=alpha,**kwargs)
	else:
		# Plot spectrum data
		if unc is True:
			line = ax.errorbar(x=spec_data['ENERGY'],y=spec_data['RATE']*y_factor,yerr=(neg_errors(spec_data['RATE']),pos_errors(spec_data['RATE']))*y_factor,label=label,fmt=" ",color=color,alpha=alpha,**kwargs)
		else:
			line = ax.errorbar(x=spec_data['ENERGY'],y=spec_data['RATE']*y_factor,label=label,fmt=" ",marker="+",color=color,alpha=alpha,**kwargs)

	# ax[1].errorbar(x=data['ENERGY'], y=resids2, yerr = 1 ,fmt=" ",color="C1",alpha=0.9,marker="+",zorder=0)

	# Custom power law index annotation
	annot = ax.annotate("", xy=(-100,-100), xytext=(20,20),textcoords="offset points", bbox=dict(boxstyle="round", fc="w"),arrowprops=dict(arrowstyle="->"))

	# Plot aesthetics
	ax.set_xscale('log')
	ax.set_yscale('log')

	# Force lower bound
	# if np.max(spec_data['RATE']*y_factor) < ylims[1]:
	# 	ax.set_ylim(ylims[1]/(10**2),ylims[1])
	# else:
	# 	ax.set_ylim(ylims[1]/(10**2), np.max(spec_data['RATE']*y_factor) * 10)
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

	def onclick(event, points):
		"""
		Function used to making a line and calculating the power law index of the line
		"""

		if event.dblclick:
			# If this is the first point being clicked on, add it to the list of points
			if len(points) == 0:
				points.append([event.xdata, event.ydata])
			# If this is the second point being clicked on, add it to the list of points, make the connecting line, and show the power law index
			elif len(points) == 1:
				# Append points to list
				points.append([event.xdata, event.ydata])

				# Make connecting line
				ax.plot([points[0][0] , points[1][0]] ,[points[0][1], points[1][1]],color="m")

				# Calculate the power law index
				ratio_F = points[1][1]/points[0][1] # Ratio of flux 
				ratio_E = points[1][0]/points[0][0] # Ratio of energy
				alpha = (np.log(ratio_F) / np.log(ratio_E) ) - 2.

				# Set the position of the annotation and make it visible
				annot.xy = points[1]
				annot.set(visible = True)

				# Display the point index
				text = r"$\alpha$ = {}".format(alpha)
				annot.set_text(text)

			# If this is the third point selected, remove the line and reset the points
			elif len(points) == 2:
				# Remove the line
				ax.lines[-1].remove()

				# Hide the annotation
				annot.set(visible = False)

				# Reset points list
				points.clear()


			# Redraw the figure to implement updates
			ax.redraw_in_frame()
			plt.gcf().canvas.draw_idle()

	# Call function for making line between two points
	points = []
	plt.gcf().canvas.mpl_connect('button_press_event', lambda event: onclick(event, points) )

	plt.tight_layout()

	if save_pref is not None :
		plt.savefig('{}.png'.format(save_pref))

	return line;

##############################################################################################################################

def plot_spec(spectrum_model=None,spectrum_data=None, plot_res=False, spec_type =0, y_factor = 1, joined=False, unc=True, emin=8, emax=4e4 ,ax=None,ymin=None,ymax=None,save_pref=None,xlabel=True,ylabel=True,label=None,fontsize=14,fontweight='bold',alpha=1,linestyle="solid",linewidth=1.5,color="C0",marker=".",inc_comps=True,comp_alpha=0.8,comp_linestyle="dashed",comp_linewidth=1.,en_steps = None, **kwargs):
	"""
	Method to plot the spectra stored in a model class or a data class. 

	Unique attribute:
		plot_res = boolean, indicates whether the residuals between a model and data should be displayed.
	"""	

	if (spectrum_model is None) and (spectrum_data is None):
		print("Please supply either a model or data.")
		return;

	if (plot_res is True):
		if (spectrum_data is None) or (spectrum_model is None):
			print("To plot residuals, both a model spectrum and spectral data must be supplied.")
			return;

	# Make plot instance if it doesn't exist
	if ax is not None:
		if isinstance(ax,matplotlib.axes._axes.Axes):
			ax_spec = ax
		if isinstance(ax,np.ndarray):
			ax_spec = ax[0]
	if ax is None:
		if plot_res is True:
			fig, ax = plt.subplots(2,1,figsize=(6,8),sharex=True,gridspec_kw={'hspace':0.,'height_ratios':[2,1]})
			ax_spec = ax[0]
		else:
			ax_spec = plt.figure().gca()

	if spectrum_data is not None:
		plot_data_spec(spectrum_data, spec_type = spec_type,y_factor=y_factor , joined=joined, unc=unc, emin=None, emax=None,ax=ax_spec,ymin=ymin,ymax=ymax,xlabel=xlabel,ylabel=ylabel,label=label,fontsize=fontsize,fontweight=fontweight,alpha=alpha,linewidth=linewidth,color=color,marker=".", **kwargs)
	if spectrum_model is not None:
		plot_model_spec(spectrum_model, spec_type = spec_type,y_factor=y_factor, inc_comps=inc_comps, emin=emin, emax=emax,ax=ax_spec,ymin=ymin,ymax=ymax,xlabel=xlabel,ylabel=ylabel,label=label,fontsize=fontsize,fontweight=fontweight, alpha=alpha,linewidth=linewidth,comp_alpha=0.8,comp_linestyle="dashed",comp_linewidth=1.,en_steps = None, **kwargs)

	if plot_res is True:
		res_arr = spectrum_model.calc_residuals(spectrum_data)
		ax[1].errorbar(res_arr['ENERGY'],res_arr['RES'],yerr = 1,fmt = " ",color=color,alpha=alpha,**kwargs)

		# Plot aesthetics
		ax[1].set_xscale('log')

		if xlabel is True:
			ax[1].set_xlabel('Energy (keV)',fontsize=fontsize,fontweight=fontweight)
		if ylabel is True:
			ax[1].set_ylabel('Residuals',fontsize=fontsize,fontweight=fontweight)

		plot_aesthetics(ax[1],fontsize=fontsize,fontweight=fontweight,**kwargs)

	if save_pref is not None :
		plt.savefig('{}.png'.format(save_pref))

##############################################################################################################################

def plot_light_curve(light_curve_data,xax_units="s",y_factor=1, Tmin=None, Tmax=None, ax=None, save_pref=None, xlabel=True, ylabel=True ,label=None, fontsize=14,fontweight='bold', logscale=False,color="C0", alpha=1,step=True,linestyle="solid",linewidth=1.5,**kwargs):
	"""
	Method to plot the input light curve data files

	Attributes:
	file_name = file name which contains spectrum data points 
	z = redshift to shift the light curve to
	label = optional label for the plotted light curve 
	ax = the matplotlib.pyplot.axes instance to make the plot on
	
	Tmin, Tmax = indicates the minimum and maximum time range to plot. If None is supplied, the minimum and maximum times of the supplied data files are used

	save_pref = if not left as None, the plot will be saved and the file name will have this prefix
	xlabel, ylabel = indicate whether x- and y- labels should be included (boolean)
	fontsize, fontweight = fontsize and fontweight of the plot font and labels on the plot
	linestyle = style of the plotting line 
	logscale = boolean, Indicates whether the time x- and y- axes should be in log scale 

	xax_units = "s", indicates the units of the x axis on the light curve. ("s"==seconds, "m"==minutes, "h"==hours, d=="days")
	"""
	# Make plot instance if it doesn't exist
	if ax is None:
		ax = plt.figure.gca()
	fig = plt.gcf()

	# Unit conversion
	x_conv = 1 # if xax_unit == "s", then we don't have to change anything
	if(xax_units == "m"):
		x_conv = 60 # if xax_unit == "m"
	if(xax_units == "h"):
		x_conv = 60*60 # if xax_unit == "h" 
	if(xax_units == "d"):
		x_conv = 60*60*24 # if xax_unit == "d"

	light_curve_data['TIME']/=x_conv

	# Plot light curve data

	if step is True:
		ax.step(light_curve_data['TIME'],light_curve_data['RATE']*y_factor,label=label,marker=" ",where="mid",color=color,alpha=alpha,linestyle=linestyle,linewidth=linewidth,**kwargs)
	elif step is False:
		ax.plot(light_curve_data['TIME'],light_curve_data['RATE']*y_factor,label=label,marker=" ",color=color,alpha=alpha,linestyle=linestyle,linewidth=linewidth,**kwargs)

	# Custom power law index annotation
	annot = ax.annotate("", xy=(-100,-100), xytext=(20,20),textcoords="offset points", bbox=dict(boxstyle="round", fc="w"),arrowprops=dict(arrowstyle="->"))

	if(logscale == True):
		ax.set_yscale('log')
		ax.set_xscale('log')

	# Plot aesthetics
	# For axis labels
	ax.set_ylabel(r'Rate (ph s$^{-1}$)',fontsize=fontsize,fontweight=fontweight)

	if(xax_units == "s"):
		ax.set_xlabel('Obs Time (sec)',fontsize=fontsize,fontweight=fontweight)
	if(xax_units == "m"):
		ax.set_xlabel('Obs Time (minutes)',fontsize=fontsize,fontweight=fontweight)
	if(xax_units == "h"):
		ax.set_xlabel('Obs Time (hours)',fontsize=fontsize,fontweight=fontweight)
	if(xax_units == "d"):
		ax.set_xlabel('Obs Time (days)',fontsize=fontsize,fontweight=fontweight)
	

	# Add label names to plot if supplied
	if label is not None:
		plt.legend(fontsize=fontsize-2)

	plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)
	
	plt.tight_layout()


	def onclick(event, points):
		"""
		Function used to making a line and calculating the power law index of the line
		"""
		
		if event.dblclick:
			# If this is the first point being clicked on, add it to the list of points
			if len(points) == 0:
				points.append([event.xdata, event.ydata])
			# If this is the second point being clicked on, add it to the list of points, make the connecting line, and show the power law index
			elif len(points) == 1:
				# Append points to list
				points.append([event.xdata, event.ydata])

				# Make connecting line
				ax.plot([points[0][0] , points[1][0]] ,[points[0][1], points[1][1]],color="m")

				# Calculate the power law index
				ratio_F = points[1][1]/points[0][1] # Ratio of flux 
				ratio_E = points[1][0]/points[0][0] # Ratio of energy
				alpha = (np.log(ratio_F) / np.log(ratio_E) )

				# Set the position of the annotation and make it visible
				annot.xy = points[1]
				annot.set(visible = True)

				# Display the point index
				text = r"$\alpha$ = {}".format(alpha)
				annot.set_text(text)

			# If this is the third point selected, remove the line and reset the points
			elif len(points) == 2:
				# Remove the line
				ax.lines[-1].remove()

				# Hide the annotation
				annot.set(visible = False)

				# Reset points list
				points.clear()


			# Redraw the figure to implement updates
			ax.redraw_in_frame()
			fig.canvas.draw_idle()


	# Call function for making line between two points
	points = []
	fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, points) )

	if save_pref is not None:
		plt.savefig('figs/{}-light-curve.png'.format(save_pref))

	return fig

##############################################################################################################################
