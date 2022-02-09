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
void DataAnalysis::add_fit_func(std::string fit_func)
{
	// Add another model row
	model_param_bounds.resize(model_param_bounds.size()+1);

	// Add model
	if(fit_func == "PL")
	{
		// Add a Power Law to the function list
		this->spec_funcs.push_back(PL);
		// Set default parameter bounds: alpha, norm
		this->model_param_bounds.at( model_param_bounds.size()-1 ) = {{-2.,2.,5},{1e40,1e55,50}};
		this->best_fit_params.resize( best_fit_params.size() + 2 );

	}
	else if(fit_func == "BPL")
	{
		// Add a Broken Power Law to the function list
		this->spec_funcs.push_back(BPL);
		// Set default parameter bounds: ebreak, alpha, beta, norm
		this->model_param_bounds.at( model_param_bounds.size()-1 ) = {{1.,1e3,5},{-2.,2.,5},{-1.,-4.,5},{1e40,1e55,50}};
		this->best_fit_params.resize( best_fit_params.size() + 4 );
	}
	else if(fit_func == "Band")
	{
		// Add a Band function to the function list
		this->spec_funcs.push_back(Band);
		// Set default parameter bounds: ebreak, alpha, beta, norm
		this->model_param_bounds.at( model_param_bounds.size()-1 ) = {{1.,1e5,5},{-2.,2.,5},{-1.,-4.,5},{1e35,1e55,50}};
		this->best_fit_params.resize( best_fit_params.size() + 4 );
	}
	else if(fit_func == "BB")
	{
		// Add a Black Body to the function list
		this->spec_funcs.push_back(BB);
		// Set default parameter bounds: temperature, alpha, norm
		this->model_param_bounds.at( model_param_bounds.size()-1 ) = {{0.1,1e3,5},{-2.,2.,5},{-1.,-4.,5},{1e35,1e55,50}};
		this->best_fit_params.resize( best_fit_params.size() + 4 );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Set boundary for a specific parameter
void DataAnalysis::set_param_bound(int model_num, int param_num, double * model_param_bounds)
{
	this->model_param_bounds.at(model_num).at(param_num).at(0) = model_param_bounds[0];
	this->model_param_bounds.at(model_num).at(param_num).at(1) = model_param_bounds[1];
	this->model_param_bounds.at(model_num).at(param_num).at(2) = model_param_bounds[2];
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

// Make parameter space grid
void DataAnalysis::make_param_space_grid(std::vector<std::vector<double>> &param_combo_list)
{
	// Number of models
	size_t num_models = model_param_bounds.size();
	// Number of parameters
	size_t num_params = 0;
	// Number of total parameters
	size_t tot_num_params = 0;
	// Number of parameter combinations
	size_t num_combos=1;
	// Temporary index tracker
	size_t tmp_ind = 0;


	// Find the number of parameters and combinations
	for(size_t i=0; i < num_models; ++i)
	{
		num_params = model_param_bounds.at(i).size();
		tot_num_params += num_params; 
		for(size_t j=0; j < num_params; ++j)
		{
			num_combos *= model_param_bounds.at(i).at(j).at(2); // Multiply by the number of steps for this parameter
		}

	}

	// Create parameter grid
	std::vector<std::vector<double>> param_space_grid;
	param_space_grid.resize(tot_num_params);
	
	// Resize the parameter combination list
	param_combo_list.resize(num_combos);

	// For each model, make a parameter space for each parameter
	for(size_t i=0; i < num_models; ++i)
	{
		for(size_t j = 0; j < model_param_bounds.at(i).size(); ++j)
		{	
			std::vector<double> tmp_arr;

			// If the parameter space is reasonably small, explore the parameter space linearly.
			if( ( model_param_bounds.at(i).at(j).at(1) - model_param_bounds.at(i).at(j).at(0) ) < 100)
			{
				double step_size = ( model_param_bounds.at(i).at(j).at(1) - model_param_bounds.at(i).at(j).at(0) ) / model_param_bounds.at(i).at(j).at(2);

				for(int k = 0; k <= model_param_bounds.at(i).at(j).at(2); ++k)
				{
					tmp_arr.push_back(model_param_bounds.at(i).at(j).at(0) + (k*step_size));
				}	
			}
			else
			{
				// If the parameter space is reasonably large, explore the parameter space logarithmically
				double step_size = ( log(model_param_bounds.at(i).at(j).at(1)) - log(model_param_bounds.at(i).at(j).at(0)) ) / model_param_bounds.at(i).at(j).at(2);

				for(int k = 0; k <= model_param_bounds.at(i).at(j).at(2); ++k)
				{
					tmp_arr.push_back( model_param_bounds.at(i).at(j).at(0)*exp(k*step_size));
				}	

			}

			param_space_grid.at(tmp_ind) = tmp_arr;
			tmp_ind+=1;
		}
	}

	// To keep track of next element in each of the parameter
	unsigned int* indices = new unsigned int[tot_num_params];
	// Initialize each to zero
	for (size_t i = 0; i < tot_num_params; ++i)
	{
		indices[i] = 0;
	}

	// Find all parameter combinations
	for(size_t k=0; k<num_combos; ++k)
	{
		// Current combination
		for (size_t i = 0; i < tot_num_params; ++i)
		{
			param_combo_list.at(k).push_back( param_space_grid.at(i).at(indices[i]) );
		}
	
		// Find the "rightmost" parameter that has more elements left to use after this current element 
		int next = tot_num_params - 1;
		while (next >= 0 && (indices[next] + 1 >= param_space_grid.at(next).size()))
		{
			--next;
		}

		// If found move to next element in that array
		++indices[next];

		// For all arrays to the right of this array current index again points to first element
		for (size_t i = next + 1; i < tot_num_params; ++i)
		{
			indices[i] = 0;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Brute Force fitting method 
void DataAnalysis::BruteForce()
{
	// Make list of all parameter combinations
	std::vector<std::vector<double>> param_combo_list;
	make_param_space_grid(param_combo_list);

	// Initialize a fit statistics
	FitStats * p_curr_fit_stats = new FitStats();
	FitStats * p_best_fit_stats = new FitStats();

	// Which parameter combo is best?
	int ind_best_param_combo;

	// Make spectrum objects to use, one for the unfolded and folded spectrum
	Spectrum * p_model_spectrum_unf = new Spectrum();
	Spectrum * p_model_spectrum_fold = new Spectrum();

	// Use the first parameter combination as the initial starting point to compare to
	// For each loaded data set and response matrix
	for(size_t i=0; i < spec_list.size(); ++i)
	{
		// Set the energy axis of the unfolded spectrum to the same as the current data/response
		(*p_model_spectrum_unf).make_energ_axes( (*resp_list.at(i)).phot_energ_min, (*resp_list.at(i)).phot_energ_min, (*resp_list.at(i)).num_phot_bins, true);

		// Calculate spectrum from the models using current parameter combination
		make_model_spectrum(p_model_spectrum_unf, param_combo_list.at(i));

		// Fold the spectrum with the loaded response matrix 
		ConvolveSpectra(p_model_spectrum_fold, (*p_model_spectrum_unf), (*resp_list.at(i)) );

		// Calculate fit statistic between the folded model spectrum and the data
		// The false argument indicates that the fit_statistic value will not reset
		(*p_curr_fit_stats).calc_fit_stat( (*spec_list.at(i)), (*p_model_spectrum_fold), false);
	}
	// Record the current fit parameters
	(*p_best_fit_stats).set_fit_stat_val( (*p_curr_fit_stats).fit_stat_val ); 
	ind_best_param_combo=0;


	// For each parameter combination in the list
	for(size_t i=1; i < param_combo_list.size(); ++i)
	{
		// For each loaded data set and response matrix
		for(size_t i=0; i < spec_list.size(); ++i)
		{
			// Set the energy axis of the unfolded spectrum to the same as the current data/response
			(*p_model_spectrum_unf).make_energ_axes( (*resp_list.at(i)).phot_energ_min, (*resp_list.at(i)).phot_energ_min, (*resp_list.at(i)).num_phot_bins, true);

			// Calculate spectrum from the models using current parameter combination
			make_model_spectrum(p_model_spectrum_unf, param_combo_list.at(0));

			// Fold the spectrum with the loaded response matrix 
			ConvolveSpectra(p_model_spectrum_fold, (*p_model_spectrum_unf), (*resp_list.at(i)) );

			// Calculate fit statistic between the folded model spectrum and the data
			(*p_curr_fit_stats).calc_fit_stat( (*spec_list.at(i)), (*p_model_spectrum_fold), false);
		}

		// If the new fit statistic is better than the previous one, update the best fit parameters
		if( (*p_curr_fit_stats).fit_stat_val < (*p_best_fit_stats).fit_stat_val)
		{
			(*p_best_fit_stats).set_fit_stat_val( (*p_curr_fit_stats).fit_stat_val ); 
			ind_best_param_combo=i;
		}

		// Resent fit_stat
		(*p_curr_fit_stats).reset_fit_stat();
	}

	this->best_fit_params = param_combo_list.at(ind_best_param_combo);

	for(size_t i = 0; i < param_combo_list.at(0).size(); ++i)
	{
		cout << param_combo_list.at(ind_best_param_combo).at(i) << " ";
	}
	cout << endl;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Brute Fore fitting, but using multiple cores (uses multi processing)
void DataAnalysis::BruteForce(int cores)
{

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Calculate spectrum from the loaded model components
void DataAnalysis::make_model_spectrum(Spectrum * p_model_spectrum, std::vector<double> parameter_list)
{
	// Initialize temp values
	size_t num_params = 0;
	int index_param_list = 0;

	// Reset function list position
	auto it = spec_funcs.begin();

	// For each model
	for( size_t i=0; i < spec_funcs.size(); ++i )
	{	
		// How many parameters does this function have?
		num_params = model_param_bounds.at(i).size();
		// Make list of the parameters for just this component
		double sub_param_list[num_params];
		for(size_t j = 0; j<num_params; ++j)
		{
			sub_param_list[j] = parameter_list[index_param_list+j];
		}

		// For each energy bin
		for(int j=0; j<(*p_model_spectrum).num_energ_bins; ++j)
		{
			// Evaluate the function at the current energy using the parameters in the parameter list
			(*p_model_spectrum).spectrum_rate.at(j) += (*it)( (*p_model_spectrum).energ_mid.at(j), sub_param_list );
			// Add the rate the spectrum total
			(*p_model_spectrum).spectrum_sum += (*p_model_spectrum).spectrum_rate.at(j); 
		}

		// Advance to next function
		std::advance(it,1);
		// Advance to the next set of parameters
		index_param_list+=num_params;
	}
}