/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-31

Header file for FitStats.cpp class

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef FITSTATSCLASS_H
#define FITSTATSCLASS_H

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <numeric>

// Import Custom Libraries
#include "Spectrum.hpp"

using namespace std;

// Declare light curve class
class FitStats
{
public:
	// FitStats constructor
	FitStats();

	// FitStats member variables
	int dof;
	double fit_stat_val;
	std::vector<double> fit_stat_terms;

	// FitStats member methods
	// Calculate fit statistic 
	double calc_fit_stat(const Spectrum & observed_spectrum, const Spectrum & folded_spectrum, bool reset = true);
	// Set the fit statistic type
	void set_fit_stat_type(std::string fit_type);
	// Get fit stat type
	int get_fit_stat_type();
	// Set the fit statistic value
	void set_fit_stat_val(double val);

	// Reset fit stat
	void reset_fit_stat();

private: 
	int fit_stat_type=1; // The dictionary/definitions of each integer value is defined within the set_fit_type()

};

#endif 