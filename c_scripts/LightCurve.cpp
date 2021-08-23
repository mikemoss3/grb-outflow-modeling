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
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries
#include "cosmology.hpp"
#include "Spectrum.hpp"
#include "LightCurve.hpp"

using namespace std;

// LightCurve constructor
LightCurve::LightCurve(float energ_min, float energ_max, int num_energ_bins, float tmin, float tmax, float dt, float z)
{
	this->energ_min = energ_min;
	this->energ_max = energ_max;
	this->num_energ_bins = num_energ_bins;
	this->tmin = tmin;
	this->tmax = tmax;
	this->dt = dt;
	this->z = z;
	make_time_axis(tmin, tmax, dt);
}
// Make time axis
void LightCurve::make_time_axis(float tmin, float tmax, float dt)
{
	this->tmin = tmin;
	this->tmax = tmax;
	this->dt = dt;
	num_time_bins = (tmax-tmin)/dt;

	for(int i=0; i<num_time_bins; ++i)
	{
		lc_time.push_back(tmin+(dt*i));
	}

	std::cout << "Light curve has been set back to zeros." << "\n";
	ZeroLightCurve();
}
void LightCurve::make_time_axis(float dt)
{
	make_time_axis(tmin, tmax, dt);
}
void LightCurve::ZeroLightCurve()
{
	for(int i=0; i<num_time_bins; ++i)
	{
		lc_rate.push_back(0.);
	}
}

// Add thermal component light curve
void LightCurve::AddThermalLightCurve(std::string file_name)
{	
	// Make spectrum class
	Spectrum spectrum = Spectrum(energ_min, energ_max, num_energ_bins, tmin, tmax, z);

	// Load emission data
	spectrum.ReadThermalEmissionData(file_name);

	// For each time bin, calculate the photon rate.
	for(int i=0; i<num_time_bins-1;++i)
	{
		spectrum.spectrum_sum = 0.; // Reset summation
		// Find the spectrum sum for each emission event which occurs between (light_curve_time[i], light_curve_time[i+1]) 
		spectrum.MakeThermalSpec(lc_time.at(i), lc_time.at(i+1) );
		
		// Ensure that the light curve rate is set to zero before adding values to it.
		lc_rate.at(i) = spectrum.spectrum_sum;

		// Convert units from erg / s to keV / s
		lc_rate.at(i) *= erg_to_kev;
		// Normalize by the time bin size 
		lc_rate.at(i) /= (lc_time.at(1)-lc_time.at(0));
		// Apply distance corrections
		lc_rate.at(i) /= 4 * M_PI * pow(lum_dist(z),2);
	}
}

// Add synchrotron component light curve
void LightCurve::AddSynchLightCurve(std::string file_name)
{
	// Make spectrum class
	Spectrum spectrum = Spectrum(energ_min, energ_max, num_energ_bins, tmin, tmax, z);
	
	// Load emission data 
	spectrum.ReadSynchEmissionData(file_name);
	
	// For each time bin, calculate the photon rate.
	for(int i=0; i<num_time_bins-1;++i)
	{
		spectrum.spectrum_sum = 0.; // Reset summation
		// Find the spectrum sum for each emission event which occurs between (light_curve_time[i], light_curve_time[i+1]) 		
		spectrum.MakeSynchSpec(lc_time.at(i), lc_time.at(i+1));

		// Ensure that the light curve rate is set to zero before adding values to it.
		lc_rate.at(i) = spectrum.spectrum_sum;

		// Convert units from erg / s to keV / s
		lc_rate.at(i) *= erg_to_kev;
		// Normalize by the time bin size 
		lc_rate.at(i) /= (lc_time.at(1)-lc_time.at(0));
		// Apply distance corrections
		lc_rate.at(i) /= 4 * M_PI * pow(lum_dist(z),2);
	}
}


// Write light curve to text files
void LightCurve::WriteToTXT(std::string out_file_name)
{
	// Write the light curve to a text file
	ofstream light_curve_file; // Construct file 
	light_curve_file.open(out_file_name); // Open text file with this name
	int i=0;
	// For each time bin, write the time and count rate to the file.
	while ( i < num_time_bins)
	{
		light_curve_file << lc_time.at(i);
		light_curve_file << " ";
		light_curve_file << lc_rate.at(i);
		light_curve_file << "\n";		
		i++;
	}
	light_curve_file.close(); // Close file
}
// Write light curve to FITS file
void LightCurve::WriteToFITS(std::string out_file_name)
{

}