/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-20

This script can be used to create the light curve for given emission data.
This script calls sub-scripts to calculate the spectrum of the source for each time bin, 
the count rate is found from the summation of the spectrum within the specified energy range.

The light curve calculation is done in light_curve_funcs.cpp and the spectrum calculations are completed in spectrum_funcs.cpp

*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "fitsio.h"

// Import Custom Libraries
#include "cosmology.hpp"
#include "LightCurve.hpp"

using namespace std;

// LightCurve constructors
LightCurve::LightCurve(float tmin, float tmax, float dt)
{
	// Assign input variables to class member variables
	this->tmin = tmin;
	this->tmax = tmax;
	this->dt = dt;
	num_time_bins = (tmax-tmin)/dt;

	// Set time axis length
	set_time_axis_length(num_time_bins);
	// Make time axis
	make_time_axis();
	// Zero light curve
	ZeroLightCurve();

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

LightCurve::LightCurve(std::vector<double> lc_time, std::vector<double> lc_rate)
{
	// Define the member variables of this light curve object from the input vectors
	this->tmin = lc_time.at(0);
	this->tmax = lc_time.back();
	this->dt = lc_time.at(1) - lc_time.at(0);
	this->num_time_bins = lc_time.size(); 

	// Assign the member vectors of this light curve object to the input vectors.
	this->lc_time = lc_time;
	this->lc_rate = lc_rate;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Copy constructor 
LightCurve::LightCurve(const LightCurve& tmp_light_curve)
{
	LightCurve(tmp_light_curve.tmin, tmp_light_curve.tmax, tmp_light_curve.dt);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make time axis
void LightCurve::make_time_axis(float tmin, float tmax, float dt)
{
	this->tmin = tmin;
	this->tmax = tmax;
	this->dt = dt;
	num_time_bins = (tmax-tmin)/dt;

	set_time_axis_length(num_time_bins);
	make_time_axis();
	// std::cout << "Light curve has been set back to zeros." << "\n";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void LightCurve::make_time_axis()
{
	for(int i=0; i<num_time_bins; ++i)
	{
		lc_time.at(i) = tmin+(dt*i);
	}
	// Zero out light curve
	ZeroLightCurve();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void LightCurve::set_time_axis_length(int num_time_bins)
{
	lc_time.resize(num_time_bins);
	lc_rate.resize(num_time_bins);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void LightCurve::ZeroLightCurve()
{
	for(int i=0; i<num_time_bins; ++i)
	{
		lc_rate.at(i) = 0.;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Write light curve to text files
void LightCurve::WriteToTXT(std::string out_file_name)
{
	// Write the light curve to a text file
	ofstream light_curve_file; // Construct file 
	light_curve_file.open(out_file_name); // Open text file with this name
	
	// Set precision
	light_curve_file << std::fixed;

	int i=0;
	// For each time bin, write the time and count rate to the file.
	while ( i < num_time_bins)
	{
		light_curve_file << lc_time.at(i);
		light_curve_file << " ";
		light_curve_file << lc_rate.at(i);
		light_curve_file << "\n";		
		++i;
	}
	light_curve_file.close(); // Close file
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Write light curve to FITS file
void LightCurve::WriteToFITS(std::string out_file_name)
{
	// Maybe next time
}