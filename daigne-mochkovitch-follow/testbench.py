import numpy as np
import matplotlib.pyplot as plt

import promptjet as pj
import shelldists as sd
import spectrum as sp

import cosmologicalconstants as cc

import timeit



jet_inst = pj.PromptJet(numshells=500,dte=0.2)
jet_inst.jet_evolution(tb=0)

# Load shell distribution data 
# t0=np.genfromtxt('sim_results/t0_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
# t3e4=np.genfromtxt('sim_results/t3e4_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
# t2e5=np.genfromtxt('sim_results/t2e5_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
# t5e5=np.genfromtxt('sim_results/t5e5_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
# tfin=np.genfromtxt('sim_results/ordlor_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])

# Load spectrum data 
spec = np.genfromtxt('sim_results/ordlor_spectrum.txt',dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float),('gammar',float),('e',float),('delt',float)])


# Plot shell distribution
# fig = plt.figure()
# ax = fig.gca()
# sd.plot_lorentz_dist(ax,t0[t0['STATUS']!=0], label='t0')
# sd.plot_lorentz_dist(ax,t3e4[t3e4['STATUS']!=0], label='t3e4, s#={}'.format(len(t3e4[t3e4['STATUS']==1])))
# sd.plot_lorentz_dist(ax,t2e5[t2e5['STATUS']!=0], label='t2e5, s#={}'.format(len(t2e5[t2e5['STATUS']==1])))
# sd.plot_lorentz_dist(ax,t5e5[t5e5['STATUS']!=0], label='t5e5, s#={}'.format(len(t5e5[t5e5['STATUS']==1])))
# sd.plot_lorentz_dist(ax,tfin[tfin['STATUS']!=0], label='tfin')

# plt.legend()

# Plot spectrum
# fig = plt.figure()
# ax = fig.gca()
# sp.plot_spectrum(ax,nuFnu=True,spectrum=spec)
# fig.tight_layout()


# Plot jet emission parameters 
fig = plt.figure()
ax = fig.gca()
ax.scatter(spec['ta'],spec['e']/spec['e'][0])
# ax.scatter(spec['ta'],spec['gammar']/100)
# ax.scatter(spec['ta'],spec['asyn'])
# ax.scatter(spec['ta'],spec['gammae']/1e4)
# ax.scatter(spec['ta'],spec['Beq'],marker='.')
# ax.scatter(spec['ta'],spec['Esyn']/1000)
# ax.set_yscale('log')

plt.show()

