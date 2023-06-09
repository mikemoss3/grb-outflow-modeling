/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-12-01

Header file for TTEs.cpp class

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef TTECLASS_H
#define TTECLASS_H

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries
#include "Spectrum.hpp"
#include "LightCurve.hpp"

using namespace std;

// Declare light curve class
class TTEs
{
public:
	// TTE constructor
	TTEs();

	// Member variables
	std::vector<double> TTE_TIME; // Vector containing the time stamps for all TTE
	std::vector<float> TTE_PHA; // Vector containing the PHA for all TTE
	std::vector<float> TTE_CHANS; // Vector containing the instrument channels numbers
	std::vector<float> TTE_ELO; // Vector containing the lower energy bound for the instrument channel at the same index 
	std::vector<float> TTE_EHI; // Vector containing the upper energy bound for the instrument channel at the same index

	// Member functions
	
	// Load GRB spectrum and light curve from TTE data designated by the given file name
	int LoadTTEData(std::string file_name);
	// Resize the events vectors (i.e., the TTE TIME and PHA vectors)
	void resize_events(int num_evts);
	// Resize ebounds vectors (i.e., the TTE CHANS, ELO, and EHI vectors)
	void resize_ebounds(int num_chans);

	// Make spectrum from TTE data
	void make_spectrum(Spectrum  * & p_spectrum);
	// Make spectrum from TTE data within a certain time interval and energy range
	void make_spectrum(Spectrum  * & p_spectrum, float TSTART, float TEND, float EMIN, float EMAX);
	
	// Slices vectors between the new start and stop points
	vector<float> slicing(vector<float>& arr, int new_start, int new_end);
	vector<double> slicing(vector<double>& arr, int new_start, int new_end);
	
	// Make light curve from TTE data
	void make_light_curve(LightCurve * & p_light_curve, float dt);
	// Make light curve from TTE data within a certain time interval and energy range
	void make_light_curve(LightCurve  * & p_light_curve, float TSTART, float TEND, float EMIN, float EMAX, float dt = 0.1);


// private:
};

#endif 