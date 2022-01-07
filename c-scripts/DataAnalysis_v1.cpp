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
// void DataAnalysis::BruteForce(int cores)
// {

// }

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Fit the currently loaded data with a Synthetic GRB
void DataAnalysis::FitSpectrum_SynthGRB()
{
	/*
	Fit an observed spectrum with synthetic spectra generated from the model of Hascoet, Daigne, and Mochkovitch 2013 (e.g., a thermal + a non-thermal component).
	The synthetic spectra are folded with the supplied instrument response matrix.

	The default fitting method is brute force over the designated parameter space. 
	fit_method == 1 : Brute Force
	fit_method == 2 : Monte Carlo
	fit_method == 3 : Genetic 
	fit_method == 4 : Annealing
	*/

	// Reset currently best fit statistics
	(*p_curr_fit_stats).reset_fit_stat();

	if(flag_set_param_init == false)
	{
		std::cout << "Initial parameters have not been set.\n";
		return;
	}

	if(flag_set_param_space == false)
	{
		std::cout << "Parameter space has not been set.\n";
		return;
	}

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
void DataAnalysis::BruteForce_SynthGRB()
{
	/*
	Using the brute force method, fit the observed spectra with synthetic spectra generated from the model of Hascoet, Daigne, and Mochkovitch 2013 (e.g., a thermal + a non-thermal component).
	The synthetic spectra are folded with the supplied instrument response matrix.
	*/

	if(flag_set_param_init == false)
	{
		std::cout << "Initial parameters have not been set.\n";
		return;
	}
	if(flag_set_param_space == false)
	{
		std::cout << "Parameter space has not been set.\n";
		return;
	}

	// Initialize a SynthGRB object to calculate a spectrum for
	SynthGRB * p_template_grb = new SynthGRB(p_curr_model_params); // Make object
	(*p_template_grb).SimulateJetDynamics(); // Simulate jet dynamics
	
	// Pointer that points to a Spectrum object that will hold the convolved model spectrum
	Spectrum * p_folded_spectrum = NULL;

	////////////////////////////////////////////
	// Calculate fit statistic of initial parameters for each ObsSpec+Resp combination

	for(size_t i=0; i < resp_list.size(); ++i)
	{

		// Using current (initial) parameters, simulate jet dynamics and calculate the synthetic spectrum													
		(*p_template_grb).make_source_spectrum((*resp_list.at(i)).phot_energ_min, (*resp_list.at(i)).phot_energ_max, (*resp_list.at(i)).num_phot_bins);

		// Convolve synthetic spectrum with response matrix
		ConvolveSpectra(p_folded_spectrum, (*(*p_template_grb).p_source_spectrum), (*resp_list.at(i)),true);

		// Find fit statistic between observed data and convolved synthetic spectrum
		(*p_curr_fit_stats).calc_fit_stat((*spec_list.at(i)),(*p_folded_spectrum),false);
	}	

	// Delete
	delete p_template_grb;
	delete p_folded_spectrum;

	////////////////////////////////////////////
	// Explore every parameter combination
	
	// Until the new fitting statistic differs from the old fit statistic by 0.1%, keep refitting.
	for(size_t i1=0; i1< (*p_curr_model_params).eps_e_vec.size(); ++i1)
	{
		for(size_t i3=0; i3< (*p_curr_model_params).eps_b_vec.size(); ++i3)
		{
			for(size_t i4=0; i4< (*p_curr_model_params).zeta_vec.size(); ++i4)
			{
				for(size_t i6=0; i6< (*p_curr_model_params).E_dot_iso_vec.size(); ++i6)
				{
					for(size_t i7=0; i7< (*p_curr_model_params).theta_vec.size(); ++i7)
					{
						for(size_t i8=0; i8< (*p_curr_model_params).r_open_vec.size(); ++i8)
						{
							for(size_t i10=0; i10< (*p_curr_model_params).eps_th_vec.size(); ++i10)
							{
								for(size_t i11=0; i11< (*p_curr_model_params).sigma_vec.size(); ++i11)
								{

									// Make a model parameters object from the current combinations of parameters 
									ModelParams * p_new_model_params = new ModelParams((*p_curr_model_params).tw, (*p_curr_model_params).dte,(*p_curr_model_params).eps_e_vec.at(i1), 
										(*p_curr_model_params).eps_b_vec.at(i3), (*p_curr_model_params).zeta_vec.at(i4), (*p_curr_model_params).E_dot_iso_vec.at(i6), 
										(*p_curr_model_params).theta_vec.at(i7), (*p_curr_model_params).r_open_vec.at(i8), (*p_curr_model_params).eps_th_vec.at(i10), 
										(*p_curr_model_params).sigma_vec.at(i11), (*p_curr_model_params).LorentzDist, (*p_curr_model_params).ShellDistParamsFile);

									// Using current parameters, simulate jet dynamics and calculate the synthetic spectrum
									SynthGRB * p_template_grb = new SynthGRB(p_new_model_params);
									(*p_template_grb).SimulateJetDynamics();

									// Make pointer to convolved model spectrum
									Spectrum * p_folded_spectrum = NULL;
									
									// Find fit statistic between observed data and convolved synthetic spectrum
									FitStats * p_new_fit_stats = new FitStats();

									for(size_t i=0; i < resp_list.size(); ++i)
									{
										// Using current (initial) parameters, simulate jet dynamics and calculate the synthetic spectrum													
										(*p_template_grb).make_source_spectrum((*resp_list.at(i)).phot_energ_min, (*resp_list.at(i)).phot_energ_max, (*resp_list.at(i)).num_phot_bins);

										// Convolve synthetic spectrum with response matrix
										ConvolveSpectra(p_folded_spectrum, (*(*p_template_grb).p_source_spectrum), (*resp_list.at(i)),true);

										// Find fit statistic between observed data and convolved synthetic spectrum
										(*p_new_fit_stats).calc_fit_stat((*spec_list.at(i)),(*p_folded_spectrum),false);
									}	

									// The new fit statistic is lower than the current one:
									if((*p_new_fit_stats).fit_stat_val < (*p_curr_fit_stats).fit_stat_val)
									{
										// Update the current best fit parameters to the best ones
										(*p_curr_model_params).copy(p_new_model_params);

										// Update the current fit statistic to the new one
										p_curr_fit_stats = p_new_fit_stats;
									}

									// Remove temp objects
									delete p_template_grb;
									delete p_folded_spectrum;
									delete p_new_model_params;

								}
							}
						}
					}
				}
			}
		}
	}		
	
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////


// void DataAnalysis::BruteForce_SynthGRB(int cores)
// {
	/*
	Using the brute force method and multiprocessing, fit the observed spectra with synthetic spectra generated from the model of Hascoet, Daigne, and Mochkovitch 2013 (e.g., a thermal + a non-thermal component).
	The synthetic spectra are folded with the supplied instrument response matrix.
	*/

// 	// Check if the parameter space has been specified 
// 	if(flag_set_param_space == false)
// 	{
// 		std::cout << "Parameter space has not been set.\n";
// 		return;
// 	}

// 	// Make a list of all the parameter combinations 
	

// 	/* Initialize necessary objects*/
// 	// Make temporary fit stat object to store intermediate calculations
// 	FitStats * p_curr_fit_stats = new FitStats();
// 	// Initialize a SynthGRB object to calc a spectrum for
// 	SynthGRB * p_template_grb = new SynthGRB(p_curr_model_params);
// 	Spectrum * p_folded_spectrum = new Spectrum(instrument_response.chan_energ_min, instrument_response.chan_energ_max, instrument_response.num_chans);


// 	////////////////////////////////////////////
// 	/* Calculate fit statistic of initial parameters */ 

// 	(*p_template_grb).set_model_params(p_curr_model_params);

// 	// Using current (initial) parameters, simulate jet dynamics and calculate the synthetic spectrum													
// 	(*p_template_grb).SimulateJetDynamics();
// 	(*p_template_grb).make_source_spectrum(instrument_response.phot_energ_min, instrument_response.phot_energ_max, instrument_response.num_phot_bins);

// 	// Convolve synthetic spectrum with response matrix
// 	ConvolveSpectra(p_folded_spectrum, (*(*p_template_grb).p_source_spectrum), instrument_response,true);

// 	// Find fit statistic between observed data and convolved synthetic spectrum
// 	(*p_curr_fit_stats).calc_fit_stat(observed_spectrum,(*p_folded_spectrum));


// 	////////////////////////////////////////////
// 	/* Explore every parameter combination */


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

void DataAnalysis::set_init_params(float tw, float dte, float eps_e, float eps_b, float zeta, double E_dot_iso, float theta, float r_open, float eps_th, float sigma, std::string LorentzDist, std::string ShellDistParamsFile)
{
	/*
	Set the initial parameters for the fitting algorithm to begin at.
	Set each parameter separately, the input values will be used to construct a ModelParams object.
	*/
	this->flag_set_param_init = true;
	this->p_curr_model_params = new ModelParams(tw, dte, eps_e, eps_b, zeta, E_dot_iso, theta, r_open, eps_th, sigma, LorentzDist, ShellDistParamsFile);
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
		std::vector<float> sigma_vec)
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
	sigma_vec);
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
void DataAnalysis::set_fit_func(string func_name)
{

}
