/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-30

FitStats class which contains all the necessary methods to interface with a set of fit parameters and statistics
*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <numeric>

// Import Custom Libraries
#include "FitStats.hpp"
#include "Spectrum.hpp"

using namespace std;

// FitStats constructor
FitStats::FitStats()
{
	reset_fit_stat();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Calculate fit statistic 
double FitStats::calc_fit_stat(const Spectrum & observed_spectrum, const Spectrum & folded_spectrum, bool reset)
{
	/*
	Calculate the fit statistic value between the two supplied spectra.
	It assumed the first supplied spectrum is the observed spectrum and the second supplied spectrum is the model spectrum 
	convolved with ax instrument response matrix. This is important, because the uncertainties of the observe spectrum are used in the calculation
	of the fit statistic.

	*/
	if(reset == true)
	{
		reset_fit_stat();
	}

	this->dof += folded_spectrum.num_energ_bins; // Set degrees of freedom to the size of the observation

	// Calc chi-squared
	if(fit_stat_type == 1)
	{
		double tot_chi_val = 0.;
		for(int i=0; i < observed_spectrum.num_energ_bins; ++i)
		{
			// Calculate each chi^2 term, assuming that the error is the square root of the count rate. 
			if(observed_spectrum.spectrum_rate.at(i) > 0.)
			{
				fit_stat_terms.push_back( pow(observed_spectrum.spectrum_rate.at(i) - folded_spectrum.spectrum_rate.at(i),2.)/observed_spectrum.spectrum_unc.at(i) );
			}
			else if (folded_spectrum.spectrum_rate.at(i) > 0.)
			{
				fit_stat_terms.push_back( pow(observed_spectrum.spectrum_rate.at(i) - folded_spectrum.spectrum_rate.at(i),2.)/folded_spectrum.spectrum_unc.at(i) );
			}
			else
			{
				fit_stat_terms.push_back( 0. );
			}
			tot_chi_val += fit_stat_terms.back();
		}
		// Set the chi^2 as the sum of the individual terms 
		set_fit_stat_val(tot_chi_val);
		return tot_chi_val;
	}
	
	// Calc reduced chi-squared
	else if(fit_stat_type == 2)
	{
		float tot_chi_red = 0;
		for(int i=0; i < observed_spectrum.num_energ_bins; ++i)
		{
			// Calculate each chi^2 term, assuming that the error is the square root of the count rate. 
			if(observed_spectrum.spectrum_rate.at(i) > 0.)
			{
				fit_stat_terms.push_back( pow(observed_spectrum.spectrum_rate.at(i) - folded_spectrum.spectrum_rate.at(i),2.)/observed_spectrum.spectrum_unc.at(i)/dof );
			}
			else if (folded_spectrum.spectrum_rate.at(i) > 0.)
			{
				fit_stat_terms.push_back( pow(observed_spectrum.spectrum_rate.at(i) - folded_spectrum.spectrum_rate.at(i),2.)/folded_spectrum.spectrum_unc.at(i)/dof );
			}
			else
			{
				fit_stat_terms.push_back( 0. );
			}
		}
		// Set the chi^2 as the sum of the individual terms 
		tot_chi_red += fit_stat_terms.back();
		set_fit_stat_val(tot_chi_red);
		return tot_chi_red;
	}
	
	// Calc c-stat
	else if(fit_stat_type == 3)
	{
		set_fit_stat_val(69.);
		return 69.;
	}
	else
	{
		std::cout << "Please choose a usable fit_stat_type.\n";
		return 1.;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Set the fit statistic type
void FitStats::set_fit_stat_type(std::string fit_type_string)
{
	/*
	Set the fit statistic to use during fitting.
	Default is chi-squared

	fit_stat = 1 : chi-squared
	fit_stat = 2 : reduced chi-squared 
	fit_stat = 3 : cstat
	*/
	if(fit_type_string == "chisq")
	{
		this->fit_stat_type = 1;
	}
	else if(fit_type_string == "redchisq")
	{
		this->fit_stat_type = 2;
	}
	else if(fit_type_string == "cstat")
	{
		this->fit_stat_type = 3;
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Get fit stat type
int FitStats::get_fit_stat_type()
{
	return this->fit_stat_type;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Private method to set the fit statistic value 
void FitStats::set_fit_stat_val(double val)
{
	/*
	Set the value of the fit statistic.
	*/
	this->fit_stat_val += val;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Reset fit statistic value and the number of degree of freedom
void FitStats::reset_fit_stat()
{
	/*
	Reset fit statistic value and the number of degree of freedom
	*/
	this->dof=0;
	set_fit_stat_val(std::numeric_limits<float>::max());
	fit_stat_terms.resize(0);
}