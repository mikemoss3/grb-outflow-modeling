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
#include "Spectrum.hpp"
#include "cosmology.hpp"

using namespace std;


// Define Spectrum class and member functions 

// Spectrum constructor
Spectrum::Spectrum(float E_min, float E_max, int num_E_bins)
{
	// Set class variables
	_E_min = E_min;
	_E_max = E_max;
	_num_E_bins = num_E_bins;

	// Fill in the energy vector
	make_ENERG_arrs(true);
	// Zero out the spectrum (this also serves to set the length of the vector)
	zero_spectrum();
}
 
// Spectrum member functions
// Fill in energy vectors 
void Spectrum::make_ENERG_arrs(bool logscale)
{

	if(logscale == true)
	{
		float logdE = log( _E_max/_E_min ) / _num_E_bins;
		for(int i=0; i<=_num_E_bins;i++)
		{
			ENERG_LO.push_back( _E_min + exp(logdE*i) );
			ENERG_HI.push_back( _E_min + exp(logdE*(i+1)) );
			ENERG_MID.push_back( (ENERG_LO[i] + ENERG_HI[i]) /2 );
		}
	}
	else
	{
		float dE = (_E_max - _E_min) / _num_E_bins;
		for(int i=0; i<_num_E_bins;i++)
		{
			ENERG_LO.push_back( _E_min + (dE*i) );
			ENERG_HI.push_back( _E_min + (dE*(i+1)) );
			ENERG_MID.push_back( (ENERG_LO[i] + ENERG_HI[i]) /2 );
		}		
	}
}
// Zero out the spectrum (this also serves to set the length of the vector)
void Spectrum::zero_spectrum()
{
	for(int i=0; i<=_num_E_bins; i++)
	{
		spectrum.push_back(0);
	}
}
void Spectrum::add_component(float *comp_rate)
{
	for(int i=0; i<=_num_E_bins;i++)
	{
		spectrum[i] += comp_rate[i];
		spectrum_sum+= comp_rate[i];
	}
}

// void Spectrum::add_fluctuations()
// {

// }

// void Spectrum::write_to_FITS()
// {

// }