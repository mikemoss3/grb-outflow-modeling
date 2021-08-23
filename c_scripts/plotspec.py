import numpy as np
import matplotlib.pyplot as plt

source_data = np.genfromtxt("testsourcespec.txt",dtype=[("ENERG",float),("RATE",float)])
obs_data = np.genfromtxt("testobsspec.txt",dtype=[("ENERG",float),("RATE",float)])

fig = plt.figure()
ax = fig.gca()

ax.scatter(source_data['ENERG'],source_data['RATE'])
ax.scatter(obs_data['ENERG'],obs_data['RATE'])

ax.set_yscale('log')
ax.set_xscale('log')

plt.show()