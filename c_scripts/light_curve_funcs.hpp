/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-20

Header file for light_curve_funcs.cpp

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef LIGHTCURVEFUNCS_H
#define LIGHTCURVEFUNCS_H

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

// Declare functions
// Initialize necessary arrays and call function to calculate light curve for thermal emission
void make_thermal_LC(double light_curve_rate[], double light_curve_time[], char * filename, int time_length, float z = 1, float Emin = 0.1, float Emax = 10000);
// Initialize necessary arrays and call function to calculate light curve for synchrotron emission
void make_synch_LC(double light_curve_rate[], double light_curve_time[], char * filename, int time_length, float z = 1, float Emin = 0.1, float Emax = 10000);

#endif 
