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

# Remove previous simulation results and figures 
# os.system('rm sim_results/*.txt')
# os.system('rm sim_results/*.png')

# Initialize a jet with specific physical characteristics
jet_inst = pj.PromptJet(numshells=5000,dte=0.002,shelldist=sd.oscillatory)
# jet_inst = pj.PromptJet(numshells=5000,dte=0.002,shelldist=sd.step)

# Simulate the jet evolution
jet_inst.jet_evolution()
