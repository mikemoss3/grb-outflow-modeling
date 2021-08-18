/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-20

This file defines all the functions necessary to read in emission data from a given data file name in order to calculate a light curve for the emission.
This script then calls methods from spectrum_funcs.cpp to calculate the spectrum of the source for each time bin, 
the count rate is found from the summation of the spectrum within the specified energy range.

*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries
#include "spectrum_funcs.hpp"
#include "cosmology.hpp"
// Reference Self-header
#include "light_curve_funcs.hpp"

// Initialize necessary arrays and call function to calculate light curve for thermal emission
void make_thermal_LC(double light_curve_rate[], double light_curve_time[], char * filename, int time_length, float z, float Emin, float Emax)
{
	// Find number of lines in emission data file. This is the number of emission events.
	int num_lines = 0;
	ifstream file_spec_data(filename); // Open the data file
	string line_spec_data; // Initialize a line
	// For each line increase the num_lines iterator
	while ( getline( file_spec_data, line_spec_data) ){num_lines++;}
	// Close files and free memory 
	file_spec_data.close(); 

	// Find the number of energy bins to use when creating the spectrum.
	int num_en_bins = 20*log10(Emax/Emin);

	// Initialize arrays to store emission data 
	float te[num_lines];
	float ta[num_lines];
	float delt[num_lines];
	float Temp[num_lines];
	double Flux[num_lines];
	float Rphot[num_lines];

	// Read in emission data and store parameters
	read_in_thermal_emission_data(filename, te, ta, delt, Temp, Flux, Rphot);
	
	// Make spectrum class
	Spectrum spectrum;
	spectrum.E_min = Emin;
	spectrum.E_max = Emax;
	spectrum.num_E_bins = num_en_bins;
	spectrum.z = z;
	// Fill in the energy vector
	spectrum.make_ENERG_arrs(true);
	// Zero out the spectrum (this also serves to set the length of the vector)
	spectrum.zero_spectrum();
	spectrum.spectrum_sum=0;

	// For each time bin, calculate the photon rate.
	for(int i=0; i<time_length;i++)
	{
		spectrum.spectrum_sum = 0.; // Reset summation
		// Find the spectrum sum for each emission event which occurs between (light_curve_time[i], light_curve_time[i+1]) 
		make_thermal_spec( &spectrum, ta, Temp, Flux, delt, num_lines, light_curve_time[i],light_curve_time[i+1], false);
		
		// Ensure that the light curve rate is set to zero before adding values to it.
		light_curve_rate[i] = spectrum.spectrum_sum;

		// Convert units from erg / s to keV / s
		light_curve_rate[i] *= erg_to_kev;
		// Normalize by the time bin size 
		light_curve_rate[i] /= (light_curve_time[1]-light_curve_time[0]);
		// Apply distance corrections
		light_curve_rate[i] /= 	4 * M_PI * pow(lum_dist(z),2);

	}
}

// Initialize necessary arrays and call function to calculate light curve for synchrotron emission
void make_synch_LC(double light_curve_rate[], double light_curve_time[], char * filename, int time_length, float z, float Emin, float Emax)
{
	// Find number of lines in emission data file. This is the number of emission events.
	int num_lines = 0;
	ifstream file_spec_data(filename); // Open the data file
	string line_spec_data; // Initialize a line
	// For each line increase the num_lines iterator
	while ( getline( file_spec_data, line_spec_data) ){num_lines++;}
	// Close files and free memory 
	file_spec_data.close(); 

	// Find the number of energy bins to use when creating the spectrum.
	int num_en_bins = 20*log10(Emax/Emin);

	// Initialize arrays to store emission data 
	float te[num_lines];
	float ta[num_lines];
	float asyn[num_lines];
	float Beq[num_lines];
	float gammae[num_lines];
	float Esyn[num_lines];
	float gammar[num_lines];
	double e_diss[num_lines];
	float delt[num_lines];
	float tau[num_lines];
	float relvel[num_lines];

	// Read in emission data and store parameters
	read_in_synch_emission_data(filename, te, ta, asyn, Beq, gammae, Esyn, gammar, e_diss, delt, tau, relvel);

	// Make spectrum class
	Spectrum spectrum;
	spectrum.E_min = Emin;
	spectrum.E_max = Emax;
	spectrum.num_E_bins = num_en_bins;
	spectrum.z = z;
	// Fill in the energy vector
	spectrum.make_ENERG_arrs(true);
	// Zero out the spectrum (this also serves to set the length of the vector)
	spectrum.zero_spectrum();
	spectrum.spectrum_sum=0;

	// For each time bin, calculate the photon rate.
	for(int i=0; i<time_length;i++)
	{
		spectrum.spectrum_sum = 0.; // Reset summation
		
		// Find the spectrum sum for each emission event which occurs between (light_curve_time[i], light_curve_time[i+1]) 
		make_synch_spec(&spectrum, ta, Esyn, e_diss, delt, tau, relvel, num_lines, light_curve_time[i],light_curve_time[i+1], false);

		// Ensure that the light curve rate is set to zero before adding values to it.
		light_curve_rate[i] = spectrum.spectrum_sum;
		
		// Convert units from erg / s to keV / s
		light_curve_rate[i] *= erg_to_kev;
		// Normalize by the time bin
		light_curve_rate[i] /= (light_curve_time[1]-light_curve_time[0]);
		// Apply distance correction
		light_curve_rate[i] /= 	4 * M_PI * pow(lum_dist(z),2) ;

	}
}
