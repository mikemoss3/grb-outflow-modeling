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
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries
// #include "spectrum_funcs.hpp"

using namespace std;

// Declare Spectrum class
class Spectrum
{
public:
	// Spectrum constructor
	Spectrum(float E_min=15., float E_max=200., int num_E_bins=80.);

	// Spectrum member functions
	void make_ENERG_arrs(bool logscale = true);
	void add_component(float *comp_rate);
	// void add_fluctuations();
	// void write_to_FITS();

private:
	float _E_min;
	float _E_max;
	int _num_E_bins;

	float spectrum_sum;
	std::vector<float> spectrum;
	std::vector<float> ENERG_LO;
	std::vector<float> ENERG_MID;
	std::vector<float> ENERG_HI;

	void zero_spectrum();

};

#endif 
