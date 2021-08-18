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
Spectrum::Spectrum()
{

}
 
// Spectrum member functions
// Fill in energy vectors 
void Spectrum::make_ENERG_arrs(bool logscale)
{

	if(logscale == true)
	{
		float logdE = log( E_max/E_min ) / num_E_bins;
		for(int i=0; i< num_E_bins; i++)
		{
			ENERG_LO.push_back( E_min + exp(logdE*i) );
			ENERG_HI.push_back( E_min + exp(logdE*(i+1)) );
			ENERG_MID.push_back( (ENERG_LO.at(i) + ENERG_HI.at(i)) /2 );
		}
	}
	else
	{
		float dE = (E_max - E_min) / num_E_bins;
		for(int i=0; i<num_E_bins;i++)
		{
			ENERG_LO.push_back( E_min + (dE*i) );
			ENERG_HI.push_back( E_min + (dE*(i+1)) );
			ENERG_MID.push_back( (ENERG_LO.at(i) + ENERG_HI.at(i)) /2 );
		}		
	}
}
// Zero out the spectrum (this also serves to set the length of the vector)
void Spectrum::zero_spectrum()
{
	for(int i=0; i<num_E_bins; i++)
	{
		spectrum_dE.push_back(0.);
	}
}

std::vector<float> Spectrum::get_ENERG_MID()
{
	return ENERG_MID;
}

// void Spectrum::set_redshift(float z)
// {
// 	_z = z;
// }

// void Spectrum::add_fluctuations()
// {

// }

// void Spectrum::write_to_FITS()
// {

// }

// void Spectrum::add_component(float *comp_rate)
// {
// 	for(int i=0; i<=num_E_bins;i++)
// 	{
// 		spectrum_dE.at(i) += comp_rate.at(i);
// 		spectrum_sum+= comp_rate.at(i);
// 	}
// }

// void Spectrum::add_to_spec_sum(float val)
// {
// 	spectrum_sum += val;
// }
// void Spectrum::add_to_spec_dE(float val, int index)
// {
// 	spectrum_dE[index] += val;
// }