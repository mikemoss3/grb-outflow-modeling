/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-29

This header file initializes all the functions necessary to create a custom response function 

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef RESPONSEFUNCS_H
#define RESPONSEFUNCS_H

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries
#include "spectrum_funcs.hpp"

using namespace std;

// Define response function shape
// Identity matrix, any source spectrum will be perfectly observed
void make_rsp_identity(float rsp[], float E_src[], float E_obs [], int num_en_bins);
// Decrease as 1/DeltaE^alpha from E_true
void make_rsp_overDeltaE(float rsp[], float E_src[], float E_obs [], int num_en_bins, float alpha=2);
// Define the Gaussian PDF 
double gauss(float x, float mu, float sigma);
// Decrease as a Gaussian from E_true
void make_rsp_gauss(float **rsp, float E_src[], float E_obs [], int num_en_bins);


#endif 
