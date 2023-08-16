/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-11

Header file for LightCurve.cpp class

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef LIGHTCURVECLASS_H
#define LIGHTCURVECLASS_H

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries

using namespace std;

// Declare light curve class
class LightCurve
{
public:
	// LightCurve constructor
	LightCurve(float tmin=0, float tmax=20, float dt=0.1, bool logscale = false);
	// Explicitly provide the member variables
	LightCurve(std::vector<double> lc_time, std::vector<double> lc_rate);
	// Copy Constructor
	LightCurve(const LightCurve& tmp_light_curve);

	// LightCurve member variables 
	float tmin;
	float tmax;
	float dt;
	int num_time_bins; 

	std::vector<double> lc_rate;
	std::vector<double> lc_time;

	// LightCurve member functions
	// Make time axis
	void make_time_axis(float tmin, float tmax, float dt, bool logscale = false);
	void make_time_axis(bool logscale = false);
	// Set the lengths of the time axes
	void set_time_axis_length(int num_time_bins);
	// Set the light curve rates to zero
	void ZeroLightCurve();

	// Write light curve to text files
	void WriteToTXT(std::string out_file_name);
	// Write light curve to FITS file
	void WriteToFITS(std::string out_file_name);

// private:

};

#endif 