/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2020-09-02

This file defines useful function for cosmological calculations
*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Reference self-header
#include "cosmology.hpp"

// Luminosity distance function
double lum_dist(float z)
{

    float lum_dist_int=0.;
    float num_steps = 1000.;
    float dz = z/num_steps;
    float z_curr=0.;
    for(int i=0;i<num_steps;i++)
    {
    	z_curr = i*dz;
    	lum_dist_int += dz * ((1.+z)*c_cm/H0 ) * 1./sqrt( ((omega_m * pow(1.+z_curr,3.) )+omega_lam) );
    }

    // return in centimeters
    return lum_dist_int*3.086*pow(10.,24.);
}
