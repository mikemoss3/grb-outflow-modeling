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
#include "spectrum_funcs.hpp"
#include "light_curve_funcs.hpp"

using namespace std;

int main(int argc, char *argv[])
{
	// argv[1] = spectrum type, i.e., thermal, synchrotron, Compton
	// argv[2] = file name containing emission data
	// argv[3] , argv[4] = Tmin, Tmax. Specifies the time interval over which the spectrum will be calculated
	// argv[5] , argv[6] = Emin, Emax. Specifies the energy interval over which the spectrum will be calculated
	// argv[8] = z, the redshift of the source

	// Check if all arguments were given. If more or less were given, print to screen a description of the inputs  
	if (argc != 9)
	{
		printf("All arguments must be specified: \n");
		printf("argv[1] 		= spectrum type, i.e., thermal, synchrotron, Compton \n");
		printf("argv[2] 		= file name containing emission data \n");
		printf("argv[3] , argv[4] 	= Tmin, Tmax. Specifies the observer frame time interval over which the spectrum will be calculated \n");
		printf("argv[5] 		= dt, Time resolution of the light curve. \n");
		printf("argv[6] , argv[7] 	= Emin, Emax. Specifies the observer frame energy interval over which the spectrum will be calculated \n");
		printf("argv[8] 		= z, the redshift of the source \n");
		return 1;
	}

	// Assign the input arguments
	char *spec_type = argv[1]; // Spectrum type
	char *filename = argv[2]; // Emission data file name

	float Tmin = atof(argv[3]); // Minimum time of the light curve 
	float Tmax = atof(argv[4]); // Maximum time of the light curve 
	float dt = atof(argv[5]); // Time resolution of the light curve 
	
	float Emin = atof(argv[6]); // Minimum energy of the spectrum used to create the light curve 
	float Emax = atof(argv[7]); // Maximum energy of the spectrum used to create the light curve 

	float z = atof(argv[8]); // Redshift to the source

	// Find the number of time bins
	int time_length = (Tmax-Tmin)/dt; 
	// Create a time-axis for the light curve
	double light_curve_time[time_length]; // Initialize the axis
	for(int i=0; i<=time_length; i++)
	{
		light_curve_time[i] = Tmin+(i*dt);
	}

	// Initialize light curve count rate array
	double light_curve_rate[time_length];
	// Depending on the spectrum type, send to a different light curve making function
	// These two functions are defined in light_curve_funcs.cpp
	if ( strcmp(spec_type,"thermal")==0 )
	{
		make_thermal_LC(light_curve_rate,light_curve_time,filename,time_length,z,Emin,Emax);
	}
	else if ( strcmp(spec_type,"synchrotron")==0 )
	{
		make_synch_LC(light_curve_rate,light_curve_time,filename,time_length,z,Emin,Emax);

	}

	// Write the light curve to a text file
	ofstream light_curve_file; // Construct file 
	light_curve_file.open("./sim_results/light_curve_points.txt"); // Open text file with this name
	int i=0;
	// For each time bin, write the time and count rate to the file.
	while ( i < time_length)
	{
		light_curve_file << light_curve_time[i];
		light_curve_file << " ";
		light_curve_file << light_curve_rate[i];
		light_curve_file << "\n";		
		i++;
	}
	light_curve_file.close(); // Close file
	
	return 0;
}
