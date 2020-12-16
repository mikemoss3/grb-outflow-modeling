
import numpy as np
import matplotlib.pyplot as plt

from shell import Shell
import lorentzdists as ld
import promptjet as pj

import timeit



# Initialize a prompt jet instance
jet_inst = pj.PromptJet(numshells=5000,dt=1)
# Simulate jet dynamics
# jet_inst.jet_evolution()



# Load Lorentz distribution and Spectrum from files
ld_fn = 'sim_results/jet_shells_t1e5.txt'
spec_fn = 'sim_results/spectrum_t1e5.txt'
# Make data structure
lor_dist = ld.load_lorentz_dist(ld_fn)
jet_inst.spectrum.load_spectrum(spec_fn)
# Plot data
plt.figure()
ld.plot_lorentz_dist(lor_dist)
plt.figure()
jet_inst.spectrum.plot_spectrum()

plt.show()
