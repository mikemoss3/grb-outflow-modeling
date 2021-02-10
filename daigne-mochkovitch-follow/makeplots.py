import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter


import shelldists as sd



# Script to make the figures from Daigne and Mochkovitch 1998 


def beautify_plot(spec_ax):
	"""
	Method to customize plot appearance 
	"""

	fontsize=14
	for tick in spec_ax.xaxis.get_major_ticks():
	    tick.label1.set_fontsize(fontsize)
	    tick.label1.set_fontweight('bold')
	    tick.label2.set_fontsize(fontsize)
	    tick.label2.set_fontweight('bold')

	for tick in spec_ax.yaxis.get_major_ticks():
	    tick.label1.set_fontsize(fontsize)
	    tick.label1.set_fontweight('bold')
	    tick.label2.set_fontsize(fontsize)
	    tick.label2.set_fontweight('bold')


## Global variables
fontsize=14
fontweight='bold'

## Figure 1

# Load shell distribution at various times
t0=np.genfromtxt('sim_results/n5000_const_te/t0_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
t3e4=np.genfromtxt('sim_results/n5000_const_te/t3e4_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
t2e5=np.genfromtxt('sim_results/n5000_const_te/t2e5_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
t5e5=np.genfromtxt('sim_results/n5000_const_te/t5e5_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])

# Plot shell distribution
fig, ax = plt.subplots(2,2,sharex=True,sharey=True)
sd.plot_lorentz_dist(ax[0,0],t0[t0['STATUS']!=0], label='T=0 s',xlabel=False,fontsize=fontsize,fontweight=fontweight)

sd.plot_lorentz_dist(ax[0,1],t3e4[t3e4['STATUS']!=0], label='T = 3e4\n s#={}'.format(len(t3e4[t3e4['STATUS']==1])),xlabel=False,ylabel=False)
sd.plot_lorentz_dist(ax[0,1],t0[t0['STATUS']!=0],fontsize=fontsize,fontweight=fontweight,linestyle='dashed',xlabel=False,ylabel=False)

sd.plot_lorentz_dist(ax[1,0],t2e5[t2e5['STATUS']!=0], label='T = 2e5\n s#={}'.format(len(t2e5[t2e5['STATUS']==1])))
sd.plot_lorentz_dist(ax[1,0],t0[t0['STATUS']!=0],fontsize=fontsize,fontweight=fontweight,linestyle='dashed',xlabel=False,ylabel=False)

sd.plot_lorentz_dist(ax[1,1],t5e5[t5e5['STATUS']!=0], label='T = 5e5\n s#={}'.format(len(t5e5[t5e5['STATUS']==1])),ylabel=False)
sd.plot_lorentz_dist(ax[1,1],t0[t0['STATUS']!=0],fontsize=fontsize,fontweight=fontweight,linestyle='dashed',xlabel=False,ylabel=False)

for i in range(2):
	for j in range(2):
		ax[i,j].legend(fontsize=fontsize)
		beautify_plot(ax[i,j])
		ax[i,j].set_ylim(0,550)

fig.tight_layout()


# Load spectrum data for the two different te:
spec_const_te = np.genfromtxt('sim_results/n5000_const_te/ordlor_spectrum.txt',dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float),('gammar',float),('e',float),('delt',float)])
spec_free_te = np.genfromtxt('sim_results/n5000_free_te/ordlor_spectrum.txt',dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float),('gammar',float),('e',float),('delt',float)])

## Figure 2

# fig, ax = plt.subplots(2,1,sharex=True)

# ax[0].scatter(spec_const_te['ta'],spec_const_te['te'],marker='.')
# ax2 = ax[0].twinx()
# ax2.scatter(spec_const_te['ta'],spec_const_te['delt'],marker='.',c='orange')

# ax[0].set_ylabel(r't$_e$ (s)',fontsize=fontsize,fontweight=fontweight)
# ax[0].yaxis.set_major_formatter(StrMethodFormatter('{x:,.1e}')) # 1 decimal places
# beautify_plot(ax[0])
# ax2.set_ylabel(r'$\Delta$t (s)',fontsize=fontsize,fontweight=fontweight)
# ax2.set_ylim(0,4)
# beautify_plot(ax2)



# ax[1].scatter(spec_const_te['ta'],spec_const_te['e']/spec_const_te['e'][0],marker='.')
# ax2 = ax[1].twinx()
# ax2.scatter(spec_const_te['ta'],spec_const_te['gammar']/100,marker='.',c='orange')


# ax[1].set_ylabel('e',fontsize=fontsize,fontweight=fontweight)
# beautify_plot(ax[1])
# ax2.set_ylabel(r'$\Gamma_r$/100',fontsize=fontsize,fontweight=fontweight)
# ax2.set_ylim(0,4.2)
# beautify_plot(ax2)

# ax[1].set_xlabel(r't$_a$ (s)',fontsize=fontsize,fontweight=fontweight)

# plt.subplots_adjust(wspace=0, hspace=0)

## Figure 3

# fig, ax = plt.subplots(2,2)

# ax[0,0].scatter(spec_const_te['ta'],spec_const_te['asyn'],marker='.')
# ax[0,0].scatter(spec_free_te['ta'],spec_free_te['asyn'],marker='.',c='orange')
# ax[0,0].set_ylabel(r'$\alpha_{syn}$',fontsize=fontsize,fontweight=fontweight)
# ax[0,0].set_ylim(-0.1,1.1)
# beautify_plot(ax[0,0])


# ax[0,1].scatter(spec_const_te['ta'],spec_const_te['gammae']/1e4,marker='.')
# ax[0,1].scatter(spec_free_te['ta'],spec_free_te['gammae']/1e4,marker='.',c='orange')
# ax[0,1].set_ylabel(r'$\Gamma_e/10^4$',fontsize=fontsize,fontweight=fontweight)
# ax[0,1].set_ylim(-0.1,4)
# ax[0,1].yaxis.set_label_position("right")
# ax[0,1].yaxis.tick_right()
# beautify_plot(ax[0,1])


# ax[1,0].scatter(spec_const_te['ta'],spec_const_te['Beq'],marker='.')
# ax[1,0].scatter(spec_free_te['ta'],spec_free_te['Beq'],marker='.',c='orange')
# ax[1,0].set_ylabel(r'B$_{eq}$',fontsize=fontsize,fontweight=fontweight)
# ax[1,0].set_yscale('log')
# ax[1,0].set_ylim(0.1,2e6)
# ax[1,0].set_xlabel(r't$_{a}$ (s)',fontsize=fontsize,fontweight=fontweight)
# beautify_plot(ax[1,0])


# ax[1,1].scatter(spec_const_te['ta'],spec_const_te['Esyn'],marker='.')
# ax[1,1].scatter(spec_free_te['ta'],spec_free_te['Esyn'],marker='.',c='orange')
# ax[1,1].set_ylabel(r'E$_{syn}$ (keV)',fontsize=fontsize,fontweight=fontweight)
# ax[1,1].set_yscale('log')
# ax[1,1].set_ylim(1)
# ax[1,1].yaxis.set_label_position("right")
# ax[1,1].yaxis.tick_right()
# ax[1,1].set_xlabel(r't$_{a}$ (s)',fontsize=fontsize,fontweight=fontweight)
# beautify_plot(ax[1,1])


# plt.subplots_adjust(wspace=0, hspace=0)



















plt.show()