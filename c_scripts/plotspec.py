import numpy as np
import matplotlib.pyplot as plt


## Spectrum 

spec_tot = np.genfromtxt("source_spectrum_tot.txt",dtype=[("ENERG",float),("RATE",float)])
spec_therm = np.genfromtxt("source_spectrum_therm.txt",dtype=[("ENERG",float),("RATE",float)])
spec_synch = np.genfromtxt("source_spectrum_synch.txt",dtype=[("ENERG",float),("RATE",float)])
# obs_data = np.genfromtxt("testobsspec.txt",dtype=[("ENERG",float),("RATE",float)])

fig = plt.figure()
ax = fig.gca()

ax.scatter(spec_tot['ENERG'],np.power(spec_tot['ENERG'],2)*spec_tot['RATE'],label='Tot')
ax.scatter(spec_therm['ENERG'],np.power(spec_synch['ENERG'],2)*spec_therm['RATE'],label='Therm')
ax.scatter(spec_synch['ENERG'],np.power(spec_therm['ENERG'],2)*spec_synch['RATE'],label='Synch')

ax.set_yscale('log')
ax.set_xscale('log')

# ax.set_xlim(15,350)
ax.legend()


## Light Curve
light_curve_tot = np.genfromtxt("test_tot.txt",dtype=[("TIME",float),("RATE",float)])
light_curve_therm = np.genfromtxt("test_therm.txt",dtype=[("TIME",float),("RATE",float)])
light_curve_synch = np.genfromtxt("test_synch.txt",dtype=[("TIME",float),("RATE",float)])

fig = plt.figure()
ax = fig.gca()

# ax.scatter(source_data['ENERG'],source_data['RATE'])
ax.scatter(light_curve_therm['TIME'],light_curve_therm['RATE'],label='Therm')
ax.scatter(light_curve_synch['TIME'],light_curve_synch['RATE'],label='Synch')
ax.scatter(light_curve_tot['TIME'],light_curve_tot['RATE'],label='Tot')

ax.legend()






plt.show()
