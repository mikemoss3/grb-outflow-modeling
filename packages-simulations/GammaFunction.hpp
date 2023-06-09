/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2022-02-28

Defines the necessary scripts to calculate the incomplete gamma function. 
This script follows the procedures outlined in Numerical Recipes by Press, Teukolsky, Vetterling, and Flannery.

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef GAMMAFUNC_H
#define GAMMAFUNC_H

#include <math.h>
#include <cmath>
#include <iostream>

using namespace std;


// Returns the value of the gamma function, Gamma(xx)
float gamma_func(float xx);

// Returns the complement of the incomplete gamma function, Gamma(a,x)
float gamma_inc(float a, float x);

// Returns the incomplete gamma function Gamma(a,x) evaluated by its series representation
// Also returns Gamma(a) as gamf.
void gser(float *gamser, float a, float x, float *gamf);

// Returns the incomplete gamma function Gamma(a,x), evaluated by is complement, using the Lentz's method
// Also returns Gamma(a) as gamf.
void gcf(float *gammacf, float a, float x, float *gamf);



#endif 
