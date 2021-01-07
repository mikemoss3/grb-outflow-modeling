
import numpy as np
import matplotlib.pyplot as plt

import promptjet as pj
import shelldists as sd
import spectrum as sp

import timeit



# jet_inst = pj.PromptJet(numshells=5000,dte=0.002)
# jet_inst.jet_evolution()

# Load shell distribution data 
t0=np.genfromtxt('sim_results/t0_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
t3e4=np.genfromtxt('sim_results/t3e4_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
t2e5=np.genfromtxt('sim_results/t2e5_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
t5e5=np.genfromtxt('sim_results/t5e5_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])
# tfin=np.genfromtxt('sim_results/ordlor_shells.txt',dtype=[('RADIUS',float),('GAMMA',float),('MASS',float),('TE',float),('STATUS',float)])

# Load spectrum data 
spec = np.genfromtxt('sim_results/ordlor_spectrum.txt',dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float)])


# Plot shell distribution
# fig = plt.figure()
# ax = fig.gca()
# sd.plot_lorentz_dist(ax,t0[t0['STATUS']!=0], label='t0')
# sd.plot_lorentz_dist(ax,t3e4[t3e4['STATUS']!=0], label='t3e4, s#={}'.format(len(t3e4[t3e4['STATUS']==1])))
# sd.plot_lorentz_dist(ax,t2e5[t2e5['STATUS']!=0], label='t2e5, s#={}'.format(len(t2e5[t2e5['STATUS']==1])))
# sd.plot_lorentz_dist(ax,t5e5[t5e5['STATUS']!=0], label='t5e5, s#={}'.format(len(t5e5[t5e5['STATUS']==1])))
# # sd.plot_lorentz_dist(ax,tfin[tfin['STATUS']!=0], label='tfin')

# plt.legend()


# Plot spectrum
plt.figure()
sp.plot_spectrum(spectrum=spec)


plt.show()