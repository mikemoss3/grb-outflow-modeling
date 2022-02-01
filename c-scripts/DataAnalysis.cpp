/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-30

DataAnalysis class which contains all the necessary methods to perform data analysis and fitting between observed GRB data and synthetic grb data
*/

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
#include "FitStats.hpp"
#include "ObsGRB.hpp"
#include "DataAnalysis.hpp"

using namespace std;

// DataAnalysis constructor
DataAnalysis::DataAnalysis()
{
	/*
	Constructor for the DataAnalysis class, initializes a pointer for the current fit statistics 
	and a pointer for the model parameters associated with the current fit statistics.
	*/
	
	p_curr_fit_stats = new FitStats(); // A pointer to the current fit statistic object 
	p_curr_model_params = NULL; // A pointer to the current model parameter

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void DataAnalysis::set_init_params(ModelParams * p_input_init_params)
{
	/*
	Set the initial parameters for the fitting algorithm to begin at.
	The input parameter is a ModelParams object.
	*/
	this->flag_set_param_init = true;
	this->p_curr_model_params = p_input_init_params;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void DataAnalysis::set_init_params(float tw, float dte, float eps_e, float eps_b, float zeta, double E_dot_iso, float theta, float r_open, float eps_th, float sigma, float p, std::string LorentzDist, std::string ShellDistParamsFile)
{
	/*
	Set the initial parameters for the fitting algorithm to begin at.
	Set each parameter separately, the input values will be used to construct a ModelParams object.
	*/
	this->flag_set_param_init = true;
	this->p_curr_model_params = new ModelParams(tw, dte, eps_e, eps_b, zeta, E_dot_iso, theta, r_open, eps_th, sigma, p, LorentzDist, ShellDistParamsFile);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void DataAnalysis::set_param_space(
		std::vector<float> eps_e_vec, 
		std::vector<float> eps_b_vec, 
		std::vector<float> zeta_vec, 
		std::vector<double> E_dot_iso_vec, 
		std::vector<float> theta_vec, 
		std::vector<float> r_open_vec, 
		std::vector<float> eps_th_vec, 
		std::vector<float> sigma_vec,
		std::vector<float> p_vec)
{
	/*
	Set the parameter space for the fitting algorithm to explore. The parameter space is specified by creating a list of value the fitting algorithm should use.
	This is essentially creating a grid of parameter combinations.
	*/
	this->flag_set_param_space = true;
	(*p_curr_model_params).set_param_space( 
	eps_e_vec, 
	eps_b_vec, 
	zeta_vec, 
	E_dot_iso_vec, 
	theta_vec, 
	r_open_vec, 
	eps_th_vec, 
	sigma_vec,
	p_vec);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Set the fit statistic type
void DataAnalysis::set_fit_stat_type(std::string fit_type_string)
{
	/*
	Set the fit statistic to use during fitting.
	Default is chi-squared

	fit_stat = 1 : chi-squared
	fit_stat = 2 : reduced chi-squared 
	fit_stat = 3 : cstat
	*/
	(*p_curr_fit_stats).set_fit_stat_type(fit_type_string);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Private method to set the fitting method 
void DataAnalysis::set_fit_method(std::string fit_type_string)
{
	/*
	Set the fit method to use for fitting.
	The default fitting method is brute force over the designated parameter space. 
	fit_method == 1 : Brute Force
	fit_method == 2 : Monte Carlo
	fit_method == 3 : Genetic 
	fit_method == 4 : Annealing
	*/
	if(fit_type_string == "brute")
	{
		this->fit_method = 1;
	}
	else if(fit_type_string == "MC")
	{
		this->fit_method = 2;
	}
	else if(fit_type_string == "genetic")
	{
		this->fit_method = 3;
	}
	else if(fit_type_string == "annealing")
	{
		this->fit_method = 4;
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Set the function to use to fit the data
// void DataAnalysis::set_fit_func(double (*func)(float, double *))
void DataAnalysis::set_fit_func(std::string fit_func)
{
	if(fit_type_string == "brute")
	{
		this->spec_func = PL;
	}
	else if(fit_type_string == "MC")
	{
		this->spec_func = BPL;
	}
	else if(fit_type_string == "genetic")
	{
		this->spec_func = Band;
	}
	else if(fit_type_string == "annealing")
	{
		this->spec_func = BB;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load an observed spectrum + instrument response combination
void DataAnalysis::LoadSpecAndResp(Spectrum * p_observed_spectrum, Response * p_instrument_response)
{
	/*
	Load an observed spectrum + instrument response combination. Both will be added (by reference) to the list of
	previously loaded spectra and response matrices, if any.  
	*/

	// Add spectrum to spectrum list
	spec_list.push_back(p_observed_spectrum);
	// Add response to response list
	resp_list.push_back(p_instrument_response);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Remove observed spectrum and instrument response combination
void DataAnalysis::RemoveSpecAndResp(int index)
{
	/*
	Removes the observed spectrum + instrument response combination stored at index "index" in the list of loaded spectra+response combinations.
	
	This will not delete the spectrum and response from memory, but simply remove the references from the list.
	*/
	spec_list.erase(spec_list.begin()+index);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Fit the currently loaded data with an empirical model (e.g., a power law, broken power law, or Band function) 
void DataAnalysis::FitSpectrum()
{
	/*
	Fit an observed spectrum with an empirical model (e.g., a power law, broken power law, or Band function) 
	The empirical spectra are folded with the supplied instrument response matrix.

	The default fitting method is brute force over the designated parameter space. 
	fit_method == 1 : Brute Force
	fit_method == 2 : Monte Carlo
	fit_method == 3 : Genetic 
	fit_method == 4 : Annealing
	*/

	// Reset currently best fit statistics
	(*p_curr_fit_stats).reset_fit_stat();

	// Call particular fitting method:
	if(fit_method == 1)
	{
		// BruteForce(observed_spectrum, instrument_response);
		BruteForce();
	}
	/*
	if(fit_method == 2)
	{
		MonteCarlo(observed_spectrum, instrument_response);
	}
	if(fit_method == 3)
	{
		Genetic(observed_spectrum, instrument_response);
	}
	if(fit_method == 4)
	{
		Annealing(observed_spectrum, instrument_response);
	}
	*/

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Brute Force fitting method 
void DataAnalysis::BruteForce()
{

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Brute Fore fitting, but using multiple cores (uses multi processing)
void DataAnalysis::BruteForce(int cores)
{

}