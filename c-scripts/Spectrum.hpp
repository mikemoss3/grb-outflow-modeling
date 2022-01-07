/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-11

Header file for the Spectrum.hpp class

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef SPECCLASS_H
#define SPECCLASS_H

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <random>
#include "fitsio.h"

// Import Custom Libraries

using namespace std;

// Declare Spectrum class
class Spectrum
{
public:
	// Spectrum constructor
	Spectrum(float energ_min=50, float energ_max=300, int num_energ_bins=80);
	// Explicitly set Spectrum object variables
	Spectrum(std::vector<double> spectrum_rate, std::vector<double> spectrum_unc, std::vector<float> energ_lo, std::vector<float> energ_mid, std::vector<float> energ_hi);
	// Spectrum Copy function
	Spectrum(const Spectrum& tmp_spec);

	// Spectrum member variables
	float energ_min; // Minimum energy of the spectrum used to create the light curve 
	float energ_max; // Maximum energy of the spectrum used to create the light curve 
	int num_energ_bins; // Number of energy bins within the energy range
	
	double spectrum_sum; // Records the summation of the spectrum 
	std::vector<double> spectrum_rate; // Records the spectrum rate in each energy bin
	std::vector<double> spectrum_unc; // Records the uncertainty for the spectrum rate in each energy bin
	std::vector<float> energ_lo; // Energy bin lower limits
	std::vector<float> energ_mid; // Center of each energy bin
	std::vector<float> energ_hi; // Energy bin upper limits

	// Spectrum member functions

	// Make energy axis of the spectrum 
	void make_energ_axes(float energ_min, float energ_max, float num_energ_bins, bool logscale);
	void make_energ_axes(bool logscale);
	void set_energ_axes_length(int num_energ_bins);
	// Set the spectrum to all zeros
	void ZeroSpectrum();
	// Add an uncertainty column
	void add_unc(std::vector<double> unc_vector); // Supply an uncertainty vector, must be same length as spectrum rate vector
	void add_unc(float unc_factor); // Supply an uncertainty factor, uncertainty = rate * factor
	// Add fluctuations to the spectrum
	void add_fluctuations();

	// Get the energy axis designating the middle of each energy bin
	// std::vector<float> get_energ_mid();
	
	// Load from text file
	void LoadFromTXT(std::string in_file_name);
	// Load from FITS file
	int LoadFromFITS(std::string in_file_name);
	// Write spectrum to text file
	void WriteToTXT(std::string out_file_name);
	// Write spectrum to FITS file
	int WriteToFITS(std::string out_file_name);

// protected:
// private:

};

#endif 
