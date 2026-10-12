import numpy as np


mp = 1.6726*np.power(10., -24.) # proton mass, in g
me = 9.1093*np.power(10., -28.) # electron mass, in g
qe = 4.8032*np.power(10., -10.) # charge of the electron, in cm^3/2 g^1/2 s^-1

c_cm = 2.99792*np.power(10., 10.) # speed of light, in cm/s

sigma_T = 6.6524*np.power(10., -25.) # cm^2, Thomson cross-section 
kappa_T = 0.2 # cm^2 / g, Thomson opacity 

kb = 8.62*1*np.power(10., -5.) # ev/K, Boltzmann constant
kb_kev = 8.617*1*np.power(10., -8. ) # keV/K, Boltzmann constant

sigma_sb = 5.67*np.power(10, -5.) # erg cm^-2 K^-4 s^-1, Stephan-Boltzmann constant 

h_planck=6.6*1*np.power(10., -27.) # cm^2 g s^-1, Planck constant 
h_planck_kev = 4.136 * np.power(10., -18.) # keV Hz^-1


a = 7.566 * np.power(10., -15.) # erg cm^-3 K^-4 = g cm^-1 s^-2 K^-4, Radiation constant

omega_m = 0.3 # matter density of the universe
omega_lam = 0.7 # dark energy density of the universe
H0 = 67.4*np.power(10., 5.) # Hubbles Constant cm/s/Mpc

# Conversion factors: 
kev_to_erg = 1.6022*np.power(10., -9.)
ev_to_erg = 1.6022*np.power(10., -12.)
erg_to_kev = 6.242*np.power(10., 8.)
erg_to_ev = 6.242*np.power(10., 11.)

# Bolometric energy bounds
E_bol_min = 1e-12 # keV
E_bol_max = 1e10 # keV