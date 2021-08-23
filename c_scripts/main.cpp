/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-20

Main function to create and use Spectra, Light Curves, and Response functions

*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries
#include "Response.hpp"
#include "LightCurve.hpp"
#include "Spectrum.hpp"
#include "cosmology.hpp"

using namespace std;

int ConvolveSpectrum(Spectrum *folded_spectrum, const Spectrum & source_spectrum, const Response & response);


int main(int argc, char const *argv[])
{
	Response response = Response();
	std::string resp_file_name = "/Users/mjmoss/Research/instr-rsp-investigate/example-instr-files/sw00883832000b_1speak.rsp";
	response.LoadRespFromFile(resp_file_name);

	Spectrum spectrum = Spectrum(response.phot_energ_min, response.phot_energ_max, response.num_phot_bins, 0, 20, 1);
	std::string file_name = "../sim_results/ordlor_spectrum_synch.txt";
	float tmin = 0;
	float tmax = 20;
	spectrum.AddSynchComp(file_name, tmin, tmax);

	Spectrum folded_spectrum = Spectrum(response.chan_energ_min, response.chan_energ_max, response.num_chans, 0, 20, spectrum.z);
	
	ConvolveSpectrum(&folded_spectrum, spectrum, response);
	
	/*
	spectrum.WriteToTXT("testsourcespec.txt");
	unfolded_spectrum.WriteToTXT("testobsspec.txt");
	*/
	// LightCurve lightcurve = LightCurve();
	// lightcurve.AddSynchLightCurve("../sim_results/ordlor_spectrum_synch.txt");
	// lightcurve.WriteToTXT("test_lc.txt");

	return 0;
}

// Fold a spectrum with a response matrix
int ConvolveSpectrum(Spectrum * folded_spectrum, const Spectrum & source_spectrum, const Response & response)
{
	// Fill in the energy vector
	(*folded_spectrum).energ_lo = response.chan_energ_lo;
	(*folded_spectrum).energ_mid = response.chan_energ_mid;
	(*folded_spectrum).energ_hi = response.chan_energ_hi;

	// Check if the number of source photon energy bins is the same as the number of
	// response photon energy bins 
	if(source_spectrum.num_energ_bins != response.num_phot_bins)
	{
		std::cout << "Source spectrum and Response Matrix do not have same number of photon bins.\n";
		return 1;
	}

	/* Convolve source spectrum with instrument response matrix */
	double tmp_sum; // Keeps track of matrix multiplication sum 
	// For each instrument energy channel
	for( int i=0; i < response.num_chans; i++)
	{
		tmp_sum = 0; // Set sum to zero
		// For for each photon energy bin
		for( int j=0 ; j < source_spectrum.num_energ_bins; j++)
		{
			tmp_sum += response.prob_matrix.at(j).at(i) * source_spectrum.spectrum_rate.at(j);
		}
		
		// Write the summation to the corresponding energy channel of the folded spectrum 
		(*folded_spectrum).spectrum_rate.at(i) = tmp_sum;
	}

	return 0;
}