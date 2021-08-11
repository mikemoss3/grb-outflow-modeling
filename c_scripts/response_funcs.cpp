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

using namespace std;

// Include custom libraries and headers
#include "response_funcs.hpp"
#include "spectrum_funcs.hpp"


// Define response function shape
// Identity matrix, any source spectrum will be perfectly observed
void make_rsp_identity(float rsp[], float E_src[], float E_obs[], int num_en_bins)
{
	for(int i=0; i<=num_en_bins; i++)
	{
		for(int j=0; j<=num_en_bins; j++)
		{
			if(i==j){rsp[i][j] = 1;}
		}

	}

}
// Decrease as 1/DeltaE^alpha from E_true
void make_rsp_overDeltaE(float rsp[], float E_src[], float E_obs[], int num_en_bins, float alpha)
{
	float tmp_sum;
	for(int i=0; i<=num_en_bins; i++)
	{
		tmp_sum = 0;
		for(int j=0; j<=num_en_bins; j++)
		{
			rsp[i][j] = 1/(1+ pow(abs(E_src[j]-E_obs),2) );
			tmp_sum += 1/(1+ pow(abs(E_src[j]-E_obs),2) );
		}
		// Normalize column
		for(int j=0; j<=num_en_bins; j++)
		{
			rsp[i][j] /= tmp_sum;
		}

	}
}

// Define the Gaussian PDF 
double gauss(float x, float mu, float sigma)
{
	return (1/sigma/sqrt(2*M_PI)) * exp(-pow(x-mu,2)/2/pow(sigma,2));
}

// Decrease as a Gaussian from E_true
void make_rsp_gauss(float **rsp, float E_src[], float E_obs[], int num_en_bins)
{
	for(int i=0; i<=num_en_bins; i++)
	{
		for(int j=0; j<=num_en_bins; j++)
		{
			rsp[i][j] = gauss(E_src[j],E_obs[i],sqrt(E_obs[i]));
		}
	}
}