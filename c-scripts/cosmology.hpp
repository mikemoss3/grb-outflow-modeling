/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-20

Header file for cosmology.cpp
This file also stores all useful constants and values used by other scripts in the project

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef COSMO_H
#define COSMO_H

#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;


const double mp = 1.6726*pow(10.,-24.); // proton mass, in g
const double me = 9.1093*pow(10.,-28.); // electron mass, in g

const double c_cm = 2.99792*pow(10.,10.); // speed of light, in cm/s

const double sigma_T = 6.6524*pow(10.,-25.); // cm^2, Thomson cross-section 
const float kappa_T = 0.2; // cm^2 / g, Thomson opacity 

const double kb = 8.62*1*pow(10.,-5. ); // ev/K, Boltzmann constant
const double kb_kev = 8.617*1*pow(10.,-8. ); // keV/K, Boltzmann constant

const double sigma_sb = 5.67*pow(10,-5.); // erg cm^-2 K^-4 s^-1, Stephan-Boltzmann constant 

const double h=6.6*1*pow(10.,-27.); // cm^2 g s^-1, Planck constant 

const double a = 7.566 * pow(10.,-15.); // erg cm^-3 K^-4, Radiation constant

const float omega_m = 0.3; // matter density of the universe
const float omega_lam = 0.7; // dark energy density of the universe
const float H0 = 67.4*pow(10.,5.); // Hubbles Constant cm/s/Mpc

// Conversion factors: 
const double kev_to_erg = 1.6022*pow(10.,-9.);
const double ev_to_erg = 1.6022*pow(10.,-12.);
const double erg_to_kev = 6.242*pow(10.,8.);
const double erg_to_ev = 6.242*pow(10.,11.);

// Bolometric energy bounds
const float E_bol_min = 0.01; // keV
const float E_bol_max = 1e5; // keV

// Luminosity distance function
double lum_dist(float z);


#endif 
