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
#include "fitsio.h"

// Import Custom Libraries
#include "Response.hpp"
#include "Spectrum.hpp"
#include "spectrum_funcs.hpp"

using namespace std;

// Define response function shape
// Fold a spectrum with a response matrix
int make_folded_spectrum(Spectrum *folded_spectrum, Response response, Spectrum source_spectrum);


#endif 
