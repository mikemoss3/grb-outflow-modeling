
import numpy as np
import matplotlib.pyplot as plt

from shell import Shell
import lorentzdists as ld
import promptjet as pj

import timeit



jet_inst = pj.PromptJet(numshells=200,dt=0.1)
jet_inst.jet_evolution()
plt.figure()
jet_inst.spectrum.plot_spectrum()

plt.show()
