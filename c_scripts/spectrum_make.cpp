/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-20

This script can be used to create the spectrum for given emission data.
This script calls methods from spectrum_funcs.cpp to calculate the spectrum of the source for each energy bin.

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
#include "Spectrum.hpp"
#include "cosmology.hpp"

using namespace std;


int main(int argc, char *argv[])
{
	// argv[1] = spectrum type, i.e., thermal, synchrotron, Compton
	// argv[2] = file name containing emission data
	// argv[3] , argv[4] = Tmin, Tmax. Specifies the time interval over which the spectrum will be calculated
	// argv[5] , argv[6] = Emin, Emax. Specifies the energy interval over which the spectrum will be calculated
	// argv[7] = Number of energy bins to use in the spectrum creation
	// argv[8] = z, the redshift of the source
	// argv[9] = nuFnu (true/false), indicates if the photon spectrum (false) or the energy density (true) should be calculated.
	// argv[10] = sum_only (true/false), indicates if only the total emission should be found. (Used to make light curves and save computation time.)

	// Check if all arguments were given. If more or less were given, print to screen a description of the inputs  
	if (argc != 11)
	{
		printf(" All arguments must be specified: \n");
		printf("argv[1] 		= spectrum type, i.e., thermal, synchrotron, Compton \n");
		printf("argv[2] 		= file name containing emission data \n");
		printf("argv[3] , argv[4] 	= Tmin, Tmax. Specifies the observer frame time interval over which the spectrum will be calculated \n");
		printf("argv[5] , argv[6] 	= Emin, Emax. Specifies the observer frame energy interval over which the spectrum will be calculated \n");
		printf("argv[7] 		= Number of energy bins to use in the spectrum creation \n");
		printf("argv[8] 		= z, the redshift of the source \n");
		printf("argv[9] 		= nuFnu (true/false), indicates if the photon spectrum (false) or the energy density (true) should be calculated.\n");
		printf("argv[10] 		= sum_only (true/false), indicates if only the total emission should be found. (Used to make light curves and save computation time.) \n");
		return 1;
	}

	// Assign the input arguments
	char *spec_type = argv[1]; // Spectrum type
	char *filename = argv[2]; // Emission data file name

	float Tmin = atof(argv[3]); // Minimum time of the light curve 
	float Tmax = atof(argv[4]); // Maximum time of the light curve 
	
	float Emin = atof(argv[5]); // Minimum energy of the spectrum used to create the light curve 
	float Emax = atof(argv[6]); // Maximum energy of the spectrum used to create the light curve 
	int num_en_bins = atoi(argv[7]); // Number of energy bins within the energy range

	float z = atof(argv[8]); // Redshift to the source

	std::stringstream ss_nuFnu(argv[9]); // Initialize a string stream to read in a boolean
	std::stringstream ss_sum_only(argv[10]); // Initialize a string stream to read in a boolean
	bool nuFnu; // Initialize a boolean
	bool sum_only; // Initialize a boolean 
	ss_nuFnu >> std::boolalpha >> nuFnu; // Indicates whether a photon spectrum or a energy density spectrum should be calculated
	ss_sum_only >> std::boolalpha >> sum_only; // Indicates if only the spectrum sum should be calculated. 
	// setting sum_only = true is used to save computation time and can be used when interested in creating a light curve 

	// Find number of lines in emission data file. This is the number of emission events.
	int num_lines = 0;
	ifstream file_spec_data(filename); // Open the data file
	string line_spec_data; // Initialize a line
	// For each line increase the num_lines iterator
	while ( getline( file_spec_data, line_spec_data) ){num_lines++;}
	// Close files and free memory 
	file_spec_data.close(); 


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
		
	// If the spectrum type is thermal: 
	if ( strcmp(spec_type,"thermal") == 0 )
	{
		// Initialize arrays to store emission event data 
		float te[num_lines];
		float ta[num_lines];
		float delt[num_lines];
		float Temp[num_lines];
		double Flux[num_lines];
		float Rphot[num_lines];

		// Read in emission event data and store parameters
		read_in_thermal_emission_data(filename, te, ta, delt, Temp, Flux, Rphot);

		// Use the emission event data to create a spectrum within the defined energy range
		make_thermal_spec(&spectrum, ta, Temp, Flux, delt, num_lines, Tmin, Tmax, nuFnu);
	}
	// Else, if the spectrum type is synchrotron
	else if ( strcmp(spec_type,"synchrotron") == 0 )
	{
		// Initialize arrays to store emission event data 
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

		// Read in emission event data and store parameters
		read_in_synch_emission_data(filename, te, ta, asyn, Beq, gammae, Esyn, gammar, e_diss, delt, tau, relvel);

		// Use the emission event data to create a spectrum within the defined energy range
		make_synch_spec(&spectrum, ta, Esyn, e_diss, delt, tau, relvel, num_lines, Tmin, Tmax, nuFnu);
	}

	// Write the spectrum sum out to a file.
	// Open a file and write the total sum to file
	ofstream spec_sum_file;
	spec_sum_file.open("./sim_results/spectrum_total.txt");
	spec_sum_file << spectrum.spectrum_sum;
	spec_sum_file.close(); // Close file
	
	// If the spectrum rate per energy bin was calculated, print it to file as well
	if(sum_only == false)
	{
		ofstream spec_dE_file;
		spec_dE_file.open("./sim_results/spectrum_points.txt");
		int i=0;
		// For each energy bin, write the energy bin value and the spectrum rate to file.
		while ( i < num_en_bins)
		{
			spec_dE_file << spectrum.ENERG_MID[i];
			spec_dE_file << " ";
			spec_dE_file << spectrum.spectrum_dE[i];
			spec_dE_file << "\n";		
			i++;
		}
		spec_dE_file.close(); // Close file
	}

	return 0;
}
