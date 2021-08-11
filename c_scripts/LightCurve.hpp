/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-11

Header file for light_curve.cpp class

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef LIGHTCURVECLASS_H
#define LIGHTCURVECLASS_H

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

// Declare light curve class
class Light_Curve
{
private:
	int secret;
public:
	Light_Curve();
	
	void func();
};

#endif 