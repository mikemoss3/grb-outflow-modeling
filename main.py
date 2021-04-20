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

# jet_inst = pj.PromptJet(numshells=500,dte=0.002,shelldist=sd.oscillatory)
jet_inst = pj.PromptJet(numshells=5000,dte=0.002,shelldist=sd.step,E_dot=10e52)
jet_inst.jet_evolution(tb=0)
