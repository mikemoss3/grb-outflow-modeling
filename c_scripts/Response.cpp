/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-29

This file is used to generate custom response matrices 

*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

// Include custom libraries and headers
#include "response_funcs.hpp"
#include "spectrum_funcs.hpp"

int main(int argc, char const *argv[])
{
	// argv[1] , argv[2] = Emin, Emax. Specifies the energy interval over which the spectrum will be calculated
	
	if (argc != 11)
	{
		printf(" All arguments must be specified: \n");
		printf("argv[1] , argv[2] 	= Emin, Emax. Specifies the observer frame energy interval over which the spectrum will be calculated \n");
		return 1;
	}

	// Make Source and Observed spectrum energy array
	float Emin = atof(argv[1]);
	float Emax = atof(argv[2]);
	int num_en_bins = log10(Emax/Emin)*20;
	float E_src[num_en_bins];
	float E_obs[num_en_bins];
	make_en_axis(E_src,Emin,Emax,num_en_bins);
	make_en_axis(E_obs,Emin,Emax,num_en_bins);
	

	// Make rsp array
	float rsp[num_en_bins][num_en_bins];
	make_rsp_gauss(rsp,E_obs,E_src,num_en_bins);


	// Fill it up

	// Print it out in a FITS file (or text file maybe?)

	return 0;

}

