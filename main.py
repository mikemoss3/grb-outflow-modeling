"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-02-15


The main script to run the prompt jet simulation code

"""

import numpy as np
import matplotlib.pyplot as plt
import os

import promptjet as pj
import shelldists as sd

os.system('rm sim_results/*.txt')
os.system('rm sim_results/*.png')

# jet_inst = pj.PromptJet(numshells=5000,dte=0.002,shelldist=sd.oscillatory)
jet_inst = pj.PromptJet(numshells=5000,dte=0.002,shelldist=sd.step)
jet_inst.jet_evolution(tb=0)



# import radiation as rd
# spec_therm = np.genfromtxt('sim_results/ordlor_spectrum_therm.txt',dtype=[('te',float),('ta',float),('delt',float),('T',float),('L',float)])
# spec_synch = np.genfromtxt('sim_results/ordlor_spectrum_synch.txt',dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float),('gammar',float),('e',float),('delt',float)])

# emission_inst = rd.Emission()
# # time,en,spec = emission_inst.make_tot_spec(spec_therm=spec_therm,spec_synch=spec_synch,dt=1)
# # time,en,spec = emission_inst.make_tot_spec(spec_therm=spec_therm,dt=1)
# time,en,spec = emission_inst.make_tot_spec(spec_synch=spec_synch,dt=1)


# plt.plot(en,en**2*spec[1])
# plt.yscale('log')
# plt.xscale('log')
# plt.show()


# a = np.array([4,5,6,7,8,9])
# print(np.argmax(a>5.5))
# print(a[np.argmax(a>5.5)])