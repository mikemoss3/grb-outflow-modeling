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

using namespace std;

// Declare Spectrum class
class Spectrum
{
public:
	// Spectrum constructor
	Spectrum();

	float E_min;
	float E_max;
	int num_E_bins;
	float z;

	// Spectrum member functions
	void make_ENERG_arrs(bool logscale = true);
	void add_component(float *comp_rate);
	void zero_spectrum();
	std::vector<float> get_ENERG_MID();
	// void set_redshift(float z);
	// void add_fluctuations();
	// void write_to_FITS();
	// void add_to_spec_sum(float val);
	// void add_to_spec_dE(float val, int index);

	double spectrum_sum;
	std::vector<double> spectrum_dE;
	std::vector<float> ENERG_LO;
	std::vector<float> ENERG_MID;
	std::vector<float> ENERG_HI;

// private:

};

#endif 
