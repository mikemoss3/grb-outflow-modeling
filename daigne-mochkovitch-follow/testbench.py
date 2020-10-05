
import numpy as np
import matplotlib.pyplot as plt

from shell import Shell
import lorentzdists as ld
import promptjet as pj

import timeit



# jet_inst = pj.PromptJet(numshells=500,dt=0.2)
# jet_inst.jet_evolution()
# spectrum = jet_inst.spectrum

# Write out spectrum
# np.savetxt('simspect.txt',spectrum.spectrum)

# Load spectrum if already written
spectrum = np.genfromtxt('simspect.txt',dtype=[('te',float),('ta',float),('asyn',float),('Beq',float),('gammae',float),('Esyn',float)])
# Make spectrum bins
num_bins=1000
emin = np.min(spectrum['Esyn'])
emax = np.max(spectrum['Esyn'])
enlogbins = np.logspace(np.log10(emin),np.log10(emax),num_bins)

plt.hist(spectrum['Esyn'],bins=enlogbins)
plt.xscale('log')

# plt.figure()
# spectrum.plot_spectrum()
plt.show()

