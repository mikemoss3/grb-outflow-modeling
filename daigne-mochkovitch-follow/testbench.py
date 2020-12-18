
import numpy as np
import matplotlib.pyplot as plt

import promptjet as pj
import shelldists as sd

import timeit



jet_inst = pj.PromptJet(numshells=500,dte=0.1)

# print(jet_inst.jet_shells)

jet_inst.jet_evolution()


# x = np.ndarray(shape=4,dtype=[('first',float),('second',float),('third',bool)])
# x['first']=np.array([1,2,3,4])
# x['second']=np.array([1,2,3,1])
# x['third']=np.array([True,True,True,False])
# print(pj.is_sorted(x['first']))
# print(pj.is_sorted(x['second']))
# print(pj.is_sorted(x[x['third']==True]['second']))


# plt.figure()
# jet_inst.spectrum.plot_spectrum()

