/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-11

Header file for light_curve.cpp class

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef LIGHTCURVECLASS_H
#define LIGHTCURVECLASS_H

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries
#include "Spectrum.hpp"
#include "Response.hpp"

using namespace std;

// Declare light curve class
class LightCurve
{
public:
	// LightCurve constructor
	LightCurve(float energ_min=50, float energ_max=300, int num_energ_bins=80, float tmin=0, float tmax=20, float dt=0.1, float z=1);

	// LightCurve member variables 
	float energ_min;
	float energ_max;
	int num_energ_bins;
	float tmin;
	float tmax;
	float dt;
	float z;
	Response *p_instrument_response;

	std::vector<double> lc_rate;
	std::vector<float> lc_time;

	// LightCurve member functions
	// Make time axis
	void make_time_axis(float tmin, float tmax, float dt);
	void make_time_axis(float dt);
	// Set the light curve rates to zero
	void ZeroLightCurve();

	// Set response instrument matrix of the observation
	void set_instrument_response(Response* instrument_response);

	// Add thermal component light curve
	void AddThermalLightCurve(std::string file_name);
	// Add synchrotron component light curve
	void AddSynchLightCurve(std::string file_name);

	// Write light curve to text files
	void WriteToTXT(std::string out_file_name);
	// Write light curve to FITS file
	void WriteToFITS(std::string out_file_name);

private:
	int num_time_bins; 

};

#endif 