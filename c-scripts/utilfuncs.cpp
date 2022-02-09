/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-12-04

This file defines useful utility functions 
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
#include "Response.hpp"
#include "Spectrum.hpp"

// Reference self-header
#include "utilfuncs.hpp"


double vel(float gamma)
{
    /*Method to calculate velocity (cm/s) from Lorentz factor (gamma) */

    return c_cm * sqrt(1.-(1./pow(gamma,2.) ));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

double beta(float gamma)
{
    /*Method to calculate the beta factor from Lorentz factor (gamma), where beta == v/c */

    // Approximation 
    // return 1.-1./(2.*pow(gamma,2.));
    // Exact
    return sqrt(1. - 1./pow(gamma,2.));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int ConvolveSpectra(Spectrum * & p_folded_spectrum, const Spectrum & unfolded_spectrum, const Response & instrument_response)
{
    /*
    Convolves an unfolded spectrum with an instrument response matrix. 
    */ 

    // Define the folded spectrum to match the instrument energy channels.
    p_folded_spectrum = new Spectrum(instrument_response.chan_energ_min, instrument_response.chan_energ_max, instrument_response.num_chans);
    // Fill in the energy vector    
    (*p_folded_spectrum).energ_lo = instrument_response.chan_energ_lo;
    (*p_folded_spectrum).energ_mid = instrument_response.chan_energ_mid;
    (*p_folded_spectrum).energ_hi = instrument_response.chan_energ_hi;

    // Check if the number of source photon energy bins is the same as the number of
    // response photon energy bins 
    if(unfolded_spectrum.num_energ_bins != instrument_response.num_phot_bins)
    {
        std::cout << "Source spectrum and Response Matrix do not have same number of photon bins.\n";
        return 1;
    }

    // std::cout << instrument_response.num_phot_bins << "\n";
    // std::cout << unfolded_spectrum.num_energ_bins << "\n";
    // std::cout << instrument_response.num_chans << "\n";
    // std::cout << (*folded_spectrum).num_energ_bins << "\n";

    /* Convolve source spectrum with instrument response matrix */
    double tmp_col_sum = 0; // Keeps track of matrix multiplication sum 
    // For each instrument energy channel
    for( int i=0; i < instrument_response.num_chans; ++i)
    {
        tmp_col_sum = 0; // Reset column sum to zero

        // For for each photon energy bin
        // for( int j=0 ; j < unfolded_spectrum.num_energ_bins; ++j)
        // For Fermi Fits files, the final energy bin is not formatted as a matrix, only a single value, which fucks things up. So I just subtract one here.
        for( int j=0 ; j < unfolded_spectrum.num_energ_bins-1; ++j)
        {
            tmp_col_sum += instrument_response.prob_matrix.at(j).at(i) * unfolded_spectrum.spectrum_rate.at(j);
        }
        
        // Write the summation to the corresponding energy channel of the folded spectrum 
        (*p_folded_spectrum).spectrum_rate.at(i) = tmp_col_sum;
        // Increase spectrum sum 
        (*p_folded_spectrum).spectrum_sum += tmp_col_sum;
    }

    return 0;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int ConvolveSpectra(Spectrum * & folded_spectrum, const Spectrum & unfolded_spectrum, const Response & instrument_response, bool forcebinning)
{
    /*
    Convolves an unfolded spectrum with an instrument response matrix. 

    If forcebinning == true, the unfolded spectrum is re-binned to match the instrument response matrix detector bins.
    */ 

    if(forcebinning == true)
    {
        // Bin the unfolded source spectrum to match the bins of the instrument response 

        // Create new spectrum object to store the binned spectrum
        Spectrum binned_unfolded_spectrum = Spectrum(instrument_response.phot_energ_min, instrument_response.phot_energ_max, instrument_response.num_phot_bins);

        // Initialize variable to keep track of contribution
        double contribution = 0.;

        // Binning algorithm:
        if(unfolded_spectrum.num_energ_bins == binned_unfolded_spectrum.num_energ_bins )
        {
            // If the two have an equal number of bins, call ConvoleSpectrum
            return ConvolveSpectra(folded_spectrum, unfolded_spectrum, instrument_response);
        }
        // Depends on if their are more source energy bins or response matrix photon bins 
        else if( unfolded_spectrum.num_energ_bins > binned_unfolded_spectrum.num_energ_bins )
        {
            // For each energy bin in the source spectrum, find the bin they would be included within, sum their rate
            for(int i=0; i < binned_unfolded_spectrum.num_energ_bins; ++i)
            {
                // Reset the contribution variable
                contribution = 0.;

                // For each original source energy bin within the new binned energy bin, add its contribution
                for(int j=0; j < unfolded_spectrum.num_energ_bins; ++j)
                {
                    if( (unfolded_spectrum.energ_mid.at(j) > binned_unfolded_spectrum.energ_lo.at(i)) & (unfolded_spectrum.energ_mid.at(j) < binned_unfolded_spectrum.energ_hi.at(i) ) )
                    {
                        contribution +=  unfolded_spectrum.spectrum_rate.at(j);
                    }
                } 

                binned_unfolded_spectrum.spectrum_rate.at(i) += contribution;
            }
        }
        else // binned_unfolded_spectrum.num_energ_bins > unfolded_spectrum.num_energ_bins
        {
            // Initialize a number of bins tracker, this is used to normalize after
            int num_bins_in = 0;

            // For each source energy bin, spread it out over the bins
            for(int i=0; i < unfolded_spectrum.num_energ_bins; ++i )
            {
                // Find the number of new bins that comprise the larger source bin
                for(int j=0; j < binned_unfolded_spectrum.num_energ_bins; ++j)
                {
                    if( (binned_unfolded_spectrum.energ_mid.at(j) > unfolded_spectrum.energ_lo.at(i)) & ( binned_unfolded_spectrum.energ_mid.at(j) < unfolded_spectrum.energ_hi.at(i)) )
                    {
                        num_bins_in+=1;
                    }
                }
                // Redistribute the rate in the larger source bin over the new bins comprise it
                for(int j=0; j < binned_unfolded_spectrum.num_energ_bins; ++j)
                {
                    if( (binned_unfolded_spectrum.energ_mid.at(j) > unfolded_spectrum.energ_lo.at(i)) & ( binned_unfolded_spectrum.energ_mid.at(j) < unfolded_spectrum.energ_hi.at(i)) )
                    {
                        binned_unfolded_spectrum.spectrum_rate.at(j) = unfolded_spectrum.spectrum_rate.at(i)/num_bins_in;
                    }
                }
            }
        }

        // Now run ConvolveSpectra
        return ConvolveSpectra(folded_spectrum, binned_unfolded_spectrum, instrument_response);
    }
    else
    {
        return ConvolveSpectra(folded_spectrum, unfolded_spectrum, instrument_response);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Power Law 
double PL(float energy, double * param_list)
{
    /*
    param_list[0] = alpha, power law slope
    param_list[1] = norm, normalization of the spectral function
    */
    float alpha = param_list[0]; 
    double norm = param_list[1];

    return norm*pow(energy/100.,alpha);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Broken Power Law 
double BPL(float energy, double * param_list)
{
    /*
    Compute the broken power law spectrum at a particular energy
    
    param_list[0] = e0, the peak energy
    param_list[1] = alpha, the low energy power law slope
    param_list[2] = beta, the high energy power law slope
    param_list[3] = norm, normalization of the spectral function

    */

    double e0 = param_list[0];
    float alpha = param_list[1];
    float beta = param_list[2];
    double norm = param_list[3];

    // If the energy is below the peak energy
    if (energy < e0)
    {
        return norm*pow(energy/e0, alpha);
    }
    // If the energy is above the peak energy
    else
    {
        return norm*pow(energy/e0, beta);
        // return pow(e0/e0, alpha-beta)*pow(energy/e0, beta);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Band spectrum function form (Band et. al., 1993)
// Defaults alpha = -1, beta = -2.5
double Band(float energy, double * param_list)
{
    /*
    Compute the Band spectrum at a particular energy
    
    param_list[0] = e0, the break energy
    param_list[1] = alpha, the low energy power law slope
    param_list[2] = beta, the high energy power law slope
    param_list[3] = norm, normalization of the spectral function

    */
    
    double e0 = param_list[0];
    float alpha = param_list[1];
    float beta = param_list[2];
    double norm = param_list[3];

    // If the energy is below the break energy
    if (energy <= (alpha-beta)*e0)
    {
        return norm*pow(energy/100.,alpha) * exp(-energy/e0);
    }
    // If the energy is above the peak energy
    else
    {
        return norm*pow((alpha-beta)*e0/100.,alpha-beta) * exp(beta-alpha) * pow(energy/100.,beta);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Broadened Blackbody
double BB(float energy, double * param_list)
{
    /*
    param_list[0] = temp, temperature of the black body (in keV)
    param_list[1] = alpha, low energy power law slope of the black body
    param_list[2] = norm, normalization of the spectral function
    */

    double temp = param_list[0]*kb_kev;
    float alpha = param_list[1];
    double norm = param_list[2];

    return norm*pow(energy/temp,1.+alpha)/(exp(energy/temp) - 1.);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Synchrotron spectrum function
double Synchrotron(float energy, double * param_list)
{
    /*
    Compute the synchrotron spectrum at a particular energy 
    
    param_list[0] = nu_c, Hz, the frequency associated with the critical Lorentz factor of the electron population
    param_list[1] = nu_m, Hz, the frequency associated with the minimum Lorentz factor of the electron population
    param_list[2] = p, the power law index of the electron population distribution
    param_list[3] = norm, normalization of the spectral function

    */
    
    double nu_c = param_list[0];
    double nu_m = param_list[1];
    float p = param_list[2];
    double norm = param_list[3];


    // Decide if the electron population is in the slow- of fast- cooling regimes
    if(nu_m < nu_c)
    {
        // Slow-cooling

        // If the energy is below h*nu_m,
        if(energy < h_planck*nu_m)
        {
            return pow(energy,-2./3.);
        }
        else if(energy < h_planck*nu_c)
        {
            // If the energy is between h*nu_m and h*nu_c
            return pow(energy,-(p+1.)/2.);
        }
        else
        {
            // If the energy is above h*nu_c
            return pow(energy,-(p+2.)/2.);
        }   
    }
    else
    {
        // Fast-cooling

        // If the energy is below h*nu_c,
        if(energy < h_planck*nu_m)
        {
            return pow(energy,-2./3.);
        }
        else if(energy < h_planck*nu_c)
        {
            // If the energy is between h*nu_c and h*nu_m
            return pow(energy,-3./2.);
        }
        else
        {
            // If the energy is above h*nu_m
            return pow(energy,-(p+2.)/2.);
        }
    }
}