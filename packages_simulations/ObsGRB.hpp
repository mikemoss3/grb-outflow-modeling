/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-11

Header file for GRB.cpp class

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef GRBCLASS_H
#define GRBCLASS_H

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries
#include "Spectrum.hpp"
#include "LightCurve.hpp"
#include "TTEs.hpp"

using namespace std;

// Declare light curve class
class ObsGRB
{
public:
	// GRB constructor
	ObsGRB();

	// GRB member variables 
	Spectrum * p_obs_spectrum; // Pointer to the source spectrum
	LightCurve * p_obs_light_curve; // Pointer to the observer light curve
	TTEs * p_TTEs; // Pointer to the Time Tagged Event object 

	// GRB member function
	
	// Load GRB spectrum and light curve from GBM TTE data designated by the given file name
	int LoadTTEData(std::string file_name);
	
	// Make spectrum from the TTE data
	void MakeTTESpectrum();
	// Make spectrum from TTE data within a certain time interval and energy range
	void MakeTTESpectrum(float TSTART, float TEND, float EMIN, float EMAX);

	// Make light curve from the TTE data
	void MakeTTELightCurve(float dt = 0.1);
	// Make light curve from TTE data within a certain time interval and energy range
	void MakeTTELightCurve(float TSTART, float TEND, float EMIN, float EMAX, float dt = 0.1);

	/*
	// Load spectrum from FITS file specified by file_name
	void LoadSpectrumFromFITS(std::string file_name);
	// Load spectrum from text file specified by file_name
	void LoadSpectrumFromTXT(std::string file_name);

	// Load light curve from FITS file specified by file_name
	void LoadLightCurveFromFITS(std::string file_name);
	// Load light curve from text file specified by file_name
	void LoadLightCurveFromTXT(std::string file_name);
	*/

// private:
};

#endif 