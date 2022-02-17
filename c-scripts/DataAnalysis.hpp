/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-31

Header file for DataAnalysis.cpp class

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef DATAANALYSISCLASS_H
#define DATAANALYSISCLASS_H

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <list>

// Import Custom Libraries
#include "cosmology.hpp"
#include "utilfuncs.hpp"
#include "Spectrum.hpp"
#include "LightCurve.hpp"
#include "ModelParams.hpp"
#include "SynthGRB.hpp"
#include "Response.hpp"
#include "FitStats.hpp"
#include "ObsGRB.hpp"

using namespace std;

// Declare light curve class
class DataAnalysis
{
public:
	// DataAnalysis constructor
	DataAnalysis();

	// DataAnalysis member variables
	FitStats * p_curr_fit_stats; // A pointer to the current fit statistic object 
	ModelParams * p_curr_model_params; // A pointer to the current model parameter 
	std::vector<double> best_fit_params; // List of best fit parameters, will be filled after a fit is performed

	// Make an empty function that will be set to the desired spectral function to use during fitting 
	// The float argument will be the energy to evaluate the spectrum at, 
	// the pointer points to a list of doubles/floats which are the input parameters to the function
	typedef double (*functiontype)(float, double * ); 
	std::list<functiontype> spec_funcs;
	
	// Initialize the number of parameters to be given, and the parameter array
	std::vector<std::vector<std::vector<double>>> model_param_bounds;

	std::vector<Spectrum*> spec_list; // Vector containing pointers to all loaded spectra
	std::vector<Response*> resp_list; // Vector containing pointers to all loaded response matrices

	// DataAnalysis member methods	
	// Set the fit statistic type
	void set_fit_stat_type(std::string fit_type_string);
	// Set the fitting method 
	void set_fit_method(std::string fit_type_string);
	
	// Set fitting function 
	// void set_fit_func( double (*functiontype)(float, double * ) );
	void add_fit_func(std::string fit_func);
	// Set boundary for a specific parameter
	void set_param_bound(int model_num, int param_num, double * param_bounds);

	// Load an observed spectrum + instrument response combination
	void LoadSpecAndResp(Spectrum * observed_spectrum, Response * instrument_response);
	// Remove observed spectrum and instrument response combination
	void RemoveSpecAndResp(int index);

	// Fitting methods using empirical functions
	// Fit the currently loaded data with an empirical model (e.g., a power law, broken power law, or Band function)
	void FitSpectrum();
	// Make parameter space grid
	void make_param_space_grid(std::vector<std::vector<double>> &param_combo_list);
	// Brute Force fitting method 
	void BruteForce();
	// Brute Fore fitting, but using multiple cores (uses multi processing)
	void BruteForce(int cores);

	// Calculate spectrum from the loaded model components
	void make_model_spectrum(Spectrum * p_model_spectrum, std::vector<double> parameter_list);

private:
	int fit_method = 1;
	bool flag_set_param_init = false;
	bool flag_set_param_space = false;
};

#endif 
