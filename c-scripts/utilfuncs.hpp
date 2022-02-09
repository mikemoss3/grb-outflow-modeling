/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-20

Header file for utilfuncs.cpp

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries
#include "cosmology.hpp"
#include "Response.hpp"
#include "Spectrum.hpp"

using namespace std;

// Method to calculate velocity (cm/s) from Lorentz factor (gamma)
double vel(float gamma);
// Method to calculate the beta factor from Lorentz factor (gamma), where beta == v/c
double beta(float gamma);

// Convolve a spectrum with an instrument response matrix
int ConvolveSpectra(Spectrum * & folded_spectrum, const Spectrum & unfoleded_spectrum, const Response & instrument_response);
// Convolve a spectrum with an instrument response matrix, changes the bins of the unfolded spectrum to be compatible with the instrument response matrix 
int ConvolveSpectra(Spectrum * & folded_spectrum, const Spectrum & unfoleded_spectrum, const Response & instrument_response, bool forcebinning);

// Useful empirical spectral functions
// Power Law
double PL(float energy, double * param_list);
// Broken Power Law
double BPL(float energy, double * param_list);
// Band spectrum function form (Band et. al., 1993)
double Band(float energy, double * param_list);
// Broadened Blackbody
double BB(float energy, double * param_list);
// Synchrotron spectrum
double Synchrotron(float energy, double * param_list);

#endif 
