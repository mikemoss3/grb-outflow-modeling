/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-29

This file defines all the functions necessary to create a custom response function 

*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <random>
#include "fitsio.h"

using namespace std;

// Include custom libraries and headers
#include "Response.hpp"
#include "response_funcs.hpp"
#include "Spectrum.hpp"
#include "spectrum_funcs.hpp"


// Define response function shape
// Fold a spectrum with a response matrix
int make_folded_spectrum(Spectrum *folded_spectrum, Response response, Spectrum source_spectrum)
{
	// Correctly assign values for folded spectrum object
	(*folded_spectrum).E_min = response.ECHAN_LO.at(0);
	(*folded_spectrum).E_max = response.ECHAN_HI.back();
	(*folded_spectrum).num_E_bins = response.num_chans;
	(*folded_spectrum).z = source_spectrum.z;

	// Fill in the energy vector
	(*folded_spectrum).ENERG_LO = response.ECHAN_LO;
	(*folded_spectrum).ENERG_MID = response.ECHAN_MID;
	(*folded_spectrum).ENERG_HI = response.ECHAN_HI;
	// Zero out the spectrum (this also serves to set the length of the vector)
	(*folded_spectrum).zero_spectrum();
	(*folded_spectrum).spectrum_sum=0;
	(*folded_spectrum).num_E_bins = response.num_chans;

	// Check if the number of source photon energy bins is the same as the number of
	// response photon energy bins 
	if(source_spectrum.num_E_bins != response.num_phot_bins)
	{
		std::cout << "Source spectrum and Response Matrix do not have same number of photon bins.\n";
		return 1;
	}

	/* Convolve source spectrum with instrument response matrix */
	float tmp_sum; // Keeps track of matrix multiplication sum 
	// For each instrument energy channel
	for( int i=0; i < response.num_chans; i++)
	{
		tmp_sum = 0; // Set sum to zero
		// For for each photon energy bin
		for( int j=0 ; j < source_spectrum.num_E_bins; j++)
		{
			// tmp_sum += response.MATRIX[j][i] * source_spectrum.spectrum_dE.at(j);
			tmp_sum += response.MATRIX.at(j).at(i) * source_spectrum.spectrum_dE.at(j);
		}
		
		// Write the summation to the corresponding energy channel of the folded spectrum 
		(*folded_spectrum).spectrum_dE.at(i) = tmp_sum;
	}

	return 0;
}