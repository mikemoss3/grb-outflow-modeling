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

// Import Custom Libraries
#include "cosmology.hpp"
#include "utilfuncs.hpp"
#include "Spectrum.hpp"
#include "LightCurve.hpp"
#include "ModelParams.hpp"
#include "SynthGRB.hpp"
#include "SynthGRBLibrary.hpp"
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

	// Make an empty function that will be set to the desired spectral function to use during fitting 
	// The float argument will be the energy to evaluate the spectrum at, 
	// the pointer points to a list of doubles/floats which are the input parameters to the function
	typedef double (*functiontype)(float, double * ); 
	functiontype spec_func = NULL;
	// Initialize the number of parameters to be given, and the parameter array
	// int num_elems = 0;
	float * param_arr;

	std::vector<Spectrum*> spec_list; // Vector containing pointers to all loaded spectra
	std::vector<Response*> resp_list; // Vector containing pointers to all loaded response matrices

	// DataAnalysis member methods	
	// Set the fit statistic type
	void set_fit_stat_type(std::string fit_type_string);
	// Set the fitting method 
	void set_fit_method(std::string fit_type_string);
	
	// Set fitting function 
	// void set_fit_func( double (*functiontype)(float, double * ) );
	void set_fit_func(std::string fit_func);

	// Set the initial parameter set
	void set_init_params(ModelParams * p_input_init_params);
	void set_init_params(float tw, float dte, float eps_e, float eps_b, float zeta, double E_dot_iso, float theta, float r_open, float eps_th, float sigma, float p, std::string LorentzDist, std::string ShellDistParamsFile);
	// Set parameter space 
	void set_param_space(
		std::vector<float> eps_e_vec, 
		std::vector<float> eps_b_vec, 
		std::vector<float> zeta_vec, 
		std::vector<double> E_dot_iso_vec, 
		std::vector<float> theta_vec, 
		std::vector<float> r_open_vec, 
		std::vector<float> eps_th_vec, 
		std::vector<float> sigma_vec,
		std::vector<float> p_vec);
	void set_param_space(float * param_space);


	// Load an observed spectrum + instrument response combination
	void LoadSpecAndResp(Spectrum * observed_spectrum, Response * instrument_response);
	// Remove observed spectrum and instrument response combination
	void RemoveSpecAndResp(int index);

	// Fitting methods using empirical functions
	// Fit the currently loaded data with an empirical model (e.g., a power law, broken power law, or Band function)
	void FitSpectrum();
	// Brute Force fitting method 
	void BruteForce();
	// Brute Fore fitting, but using multiple cores (uses multi processing)
	void BruteForce(int cores);

private:
	int fit_method = 1;
	bool flag_set_param_init = false;
	bool flag_set_param_space = false;
};

#endif 
