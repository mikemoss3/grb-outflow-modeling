/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-30

ObsGRB class which holds and interfaces with LightCurve and Spectrum which will hold observation data.
*/

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
#include "ObsGRB.hpp"

using namespace std;

// ObsGRB constructor
ObsGRB::ObsGRB()
{
	// Assign pointers to NULL that point to the spectrum and light curve objects associated with this GRB.
	p_obs_spectrum = NULL;
	p_obs_light_curve = NULL;
	p_TTEs = NULL;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load GRB spectrum and light curve from TTE data designated by the given file name
int ObsGRB::LoadTTEData(std::string file_name)
{
	/*
	This method is used to load a spectrum and light curve from a FITS file defined by file_name.

	The given FITS file should be of the standard form of the TTE data reported by the Fermi/GBM instrument.
	*/

	// Initialize TTE object
	p_TTEs = new TTEs();

	// Fill TTE object by loading TTE data (using given file name)
	int status = 0; // Status marker for loading FITS data
	status = (*p_TTEs).LoadTTEData(file_name);


	// Make Spectrum
	// (*p_TTEs).make_spectrum(p_obs_spectrum);
	MakeTTESpectrum();

	// Make Light Curve
	// (*p_TTEs).make_light_curve(p_obs_light_curve, 0.1);
	MakeTTELightCurve(0.1);

	return status;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make spectrum for the TTE data using complete time interval and energy range
void ObsGRB::MakeTTESpectrum()
{
	/*
	Make spectrum for the TTE data for events within the entire time range and energy range.
	
	*/

	(*p_TTEs).make_spectrum(p_obs_spectrum);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make spectrum for the TTE data
void ObsGRB::MakeTTESpectrum(float TSTART, float TEND, float EMIN, float EMAX)
{
	/*
	Make spectrum for the TTE data for events within the given time range and energy range.
	
	Attributes:
	TSTART, TEND = time interval bounds, all events that occur within this time interval will be used to make the spectrum
	EMIN, EMAX = energy band bounds, all events that occur within this energy band will be used to make the spectrum
	*/

	// Make Spectrum
	(*p_TTEs).make_spectrum(p_obs_spectrum,TSTART,TEND,EMIN,EMAX);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make Light Curve for the TTE data using complete time interval and energy range
void ObsGRB::MakeTTELightCurve(float dt)
{
	/*
	Make Light Curve for the TTE data for events within the entire time range and energy range.
	
	*/

	(*p_TTEs).make_light_curve(p_obs_light_curve, dt);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make Light Curve for the TTE data
void ObsGRB::MakeTTELightCurve(float TSTART, float TEND, float EMIN, float EMAX, float dt)
{
	/*
	Make Light Curve for the TTE data for events within the given time range and energy range.
	
	Attributes:
	TSTART, TEND = time interval bounds, all events that occur within this time interval will be used to make the spectrum
	EMIN, EMAX = energy band bounds, all events that occur within this energy band will be used to make the spectrum
	*/

	// Make Spectrum
	(*p_TTEs).make_light_curve(p_obs_light_curve, TSTART, TEND, EMIN, EMAX, dt);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
// Load spectrum from FITS file specified by file_name
void ObsGRB::LoadSpectrumFromFITS(std::string file_name)
{
}
// Load spectrum from text file specified by file_name
void ObsGRB::LoadSpectrumFromTXT(std::string file_name)
{

}

// Load light curve from FITS file specified by file_name
void ObsGRB::LoadLightCurveFromFITS(std::string file_name)
{

}
// Load light curve from text file specified by file_name
void ObsGRB::LoadLightCurveFromTXT(std::string file_name)
{
*/
