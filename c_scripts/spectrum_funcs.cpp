/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-20

This file defines all the functions necessary to read in emission data from a given data file name in order to calculate a spectrum for the emission.

*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries
#include "cosmology.hpp"
#include "Spectrum.hpp"
// Reference Self-header
#include "spectrum_funcs.hpp"

using namespace std;

// Calls function to calculate the thermal spectrum rate for each energy bin
void make_thermal_spec(Spectrum * spectrum, float ta[], float Temp[], double Flux[], float delt[], int num_lines, float Tmin, float Tmax, bool nuFnu)
{
	float z = (*spectrum).z;
	// For each emission event, calculate the emitted spectrum
	for(int i=0;i<num_lines;i++)
	{		
		// We only want to take the emission that occurs between the specified Tmin and Tmax. 
		// The emission occurs between (ta+delt), if any of it overlaps with Tmin and Tmax, calculate its contribution.
		if ( ta[i]*(1+z) <= Tmax and (ta[i]+delt[i])*(1+z) >= Tmin )
		{
			// Call function to calculate thermal spectrum rate
			calc_thermal_spec(spectrum, Temp[i],Flux[i], nuFnu);
			}
		// Else, if this emission occurred outside of the time interval, don't add it to the spectrum.
	}
}
// Read in emission data for thermal emission
void read_in_thermal_emission_data(std::string file_name, float te[], float ta[], float delt[], float Temp[], double Flux[], float Rphot[])
{
	// Load in the file
	ifstream file_spec_data(file_name);
	string line_spec_data;
	int iter=0;
	while ( getline( file_spec_data, line_spec_data) ) 
	{
		std::istringstream iss(line_spec_data);
		string col1_val;
		string col2_val;
		string col3_val;
		string col4_val;
		string col5_val;
		string col6_val;
		if ( iss >> col1_val >> col2_val >> col3_val >> col4_val >> col5_val >> col6_val )
		{
			te[iter] = stof(col1_val);
			ta[iter] = stof(col2_val);
			delt[iter] = stof(col3_val);
			Temp[iter] = stof(col4_val);
			Flux[iter] = stod(col5_val);
			Rphot[iter] = stof(col6_val);
		}
		iter++;
	}

	// Close files and free memory 
	file_spec_data.close();
}
// Calculate the thermal spectrum from the given temperature and flux of the emission
void calc_thermal_spec(Spectrum * spectrum, float temp, double flux, bool nuFnu)
{
	// Calculate the normalization
	float norm=0; // Set normalization to zero
	float Emin_bol = 0.01; // Minimum energy of the Bolometric band
	float Emax_bol = 1e5; // Maximum energy of the Bolometric band
	int norm_num_bin=40*log(Emax_bol/Emin_bol); // Number of energy bins to use to calculate normalization
	float norm_energy_axis[norm_num_bin+1]; // Initialize energy axis for the normalization, add one extra in order to use Left-Riemann-Sum
	make_en_axis(norm_energy_axis,Emin_bol,Emax_bol,norm_num_bin+1); // Make energy axis

	float en_curr; // Current energy to evaluate the addition to the normalization
	// For each energy bin along the normalization energy axis, calculate the addition to the normalization and add it. 
	for(int i=0;i < norm_num_bin; i++)
	{
		en_curr = norm_energy_axis[i]; // Set the current energy
		// Calculate the contribution according to a Left-Riemann-Sum
		norm += (norm_energy_axis[i+1]-en_curr) * en_curr * thermal_spectrum(en_curr,temp);
	}

	// Calculate the rate of the thermal spectrum
	double tmp_val=0; // Temporary value to store spectrum contributions
	// For each energy bin along the energy axis
	for(int i=0; i < (*spectrum).num_E_bins; i++)
	{
		en_curr = (*spectrum).ENERG_MID[i]; // Set the current energy
		
		// Check whether it was specified to calculate the energy density spectrum (nuFnu==true) or the photon spectrum (false)
		// The current temp and energy bin define the count rate, the normalization found above is applied.
		// This must still be multiplied by the flux of the source.
		if(nuFnu == true)
			{tmp_val = flux * pow(en_curr,2) * thermal_spectrum(en_curr, temp) / norm;}
		else
			{tmp_val = flux * thermal_spectrum(en_curr,temp) / norm;}
		
	
		// Add the contribution to the total spectrum according to a Left-Riemann-Sum
		(*spectrum).spectrum_sum += ((*spectrum).ENERG_HI[i] - (*spectrum).ENERG_LO[i]) * tmp_val;
		// Check if the spectrum rate per energy bin is requested and store it if so. 
		(*spectrum).spectrum_dE[i] += tmp_val;

	}
}
// Thermal spectrum function form based on a modified Planck function. 
// The modification comes in the form of setting the default for alpha = 0.4, which accounts for the observation of multiple black bodies simultaneously
double thermal_spectrum(float energy, float temp, float alpha)
{
	return pow(energy/(kb_kev*temp),1+alpha)/(exp(energy/(kb_kev*temp))-1);
}

// Calls function to calculate the synchrotron spectrum rate for each energy bin
void make_synch_spec(Spectrum * spectrum, float ta[], float Esyn[], double e_diss[], float delt[], float tau[], float relvel[], int num_lines, float Tmin, float Tmax, bool nuFnu)
{
	float z = (*spectrum).z;
	// For each emission event, calculate the emitted spectrum
	for(int i=0;i<num_lines;i++)
	{
		// We only want to take the emission that occurs between the specified Tmin and Tmax. 
		// The emission occurs between (ta+delt), if any of it overlaps with Tmin and Tmax, calculate its contribution.
		if ( ta[i]*(1+z) <= Tmax and (ta[i]+delt[i])*(1+z) >= Tmin )
		{
			// The emission will only be observable if the relativistic velocity is great than the local sound speed v_s/c = 0.1
			// And if the wind is transparent to the radiation
			if ( relvel[i] > 0.1 and tau[i] < 1)
			{
				// Call function to calculate thermal spectrum rate
				calc_synch_spec(spectrum,Esyn[i],e_diss[i],delt[i],nuFnu);	
			}
			
		}
		// Else, if this emission occurred outside of the time interval, don't add it to the spectrum. 
	}
}
// Read in emission data for synchrotron emission
void read_in_synch_emission_data(std::string file_name, float te[], float ta[], float asyn[], float Beq[], float gammae[], float Esyn[], float gammar[], double e_diss[],float delt[], float tau[], float relvel[])
{
	// Load in the file
	ifstream file_spec_data(file_name);
	string line_spec_data;
	int iter=0;
	while ( getline( file_spec_data, line_spec_data) ) 
	{
		std::istringstream iss(line_spec_data);
		string col1_val;
		string col2_val;
		string col3_val;
		string col4_val;
		string col5_val;
		string col6_val;
		string col7_val;
		string col8_val;
		string col9_val;
		string col10_val;
		string col11_val;
		if ( iss >> col1_val >> col2_val >> col3_val >> col4_val >> col5_val >> col6_val >> col7_val >> col8_val >> col9_val >> col10_val >> col11_val)
		{
			te[iter] = stof(col1_val);
			ta[iter] = stof(col2_val);
			asyn[iter] = stof(col3_val);
			Beq[iter] = stof(col4_val);
			gammae[iter] = stof(col5_val);
			Esyn[iter] = stof(col6_val);
			gammar[iter] = stof(col7_val);
			e_diss[iter] = stod(col8_val);
			delt[iter] = stof(col9_val);
			tau[iter] = stof(col10_val);
			relvel[iter] = stod(col11_val);
		}
		iter++;
	}

	// Close files and free memory 
	file_spec_data.close();
}
// Calculate the synchrotron spectrum from the synchrotron energy and flux of the emission
void calc_synch_spec(Spectrum * spectrum, float Esyn, double e_diss, float delt, bool nuFnu)
{

	// Calculate the normalization
	float norm=0; // Set normalization to zero
	float Emin_bol = 0.01; // Minimum energy of the Bolometric band
	float Emax_bol = 1e5; // Maximum energy of the Bolometric band
	int norm_num_bin=20*log(Emax_bol/Emin_bol); // Number of energy bins to use to calculate normalization
	float norm_energy_axis[norm_num_bin+1]; // Initialize energy axis for the normalization, add one extra in order to use Left-Riemann-Sum
	make_en_axis(norm_energy_axis,Emin_bol,Emax_bol,norm_num_bin+1); // Make energy axis

	float en_curr; // Current energy to evaluate the addition to the normalization
	// For each energy bin along the normalization energy axis, calculate the addition to the normalization and add it. 
	for(int i=0;i < norm_num_bin; i++)
	{
		en_curr = norm_energy_axis[i]; // Set the current energy
		// Calculate the contribution according to a Left-Riemann-Sum
		norm += (norm_energy_axis[i+1]-en_curr) * en_curr * synch_spectrum(en_curr,Esyn);
	}

	// Calculate the rate of the thermal spectrum
	double tmp_val=0; // Temporary value to store spectrum contributions
	// For each energy bin along the energy axis
	for(int i=0; i < (*spectrum).num_E_bins; i++)
	{
		float en_curr = (*spectrum).ENERG_MID[i]; // Set the current energy

		// Check whether it was specified to calculate the energy density spectrum (nuFnu==true) or the photon spectrum (false)
		// The current temp and energy bin define the count rate, the normalization found above is applied.
		// This must still be multiplied by the energy dissipated during the emission event.
		// The energy dissipated can be turned into Flux by dividing the energy dissipated by the observed emission duration (delt).
		if(nuFnu == true)
			{tmp_val = e_diss * pow(en_curr,2)* synch_spectrum(en_curr,Esyn) / norm / delt;}
		else
			{tmp_val = e_diss * synch_spectrum(en_curr,Esyn) / norm / delt;}

		// Add the contribution to the total spectrum according to a Left-Riemann-Sum
		(*spectrum).spectrum_sum += ((*spectrum).ENERG_HI[i] - (*spectrum).ENERG_LO[i]) * tmp_val;

		// Check if the spectrum rate per energy bin is requested and store it if so. 
		(*spectrum).spectrum_dE[i] += tmp_val;
	}

}
// Synchrotron spectrum function form,
// Defaults alpha = -1, beta = -2.5
double synch_spectrum(float energy, float Esyn, float alpha, float beta)
{
	return (1/Esyn) * band(energy,Esyn,alpha,beta);
}

// Band spectrum function form (Band et. al., 1993)
// Defaults alpha = -1, beta = -2.5
double band(float energy, float e0, float alpha, float beta)
{
	// Compute the power law spectrum at a particular energy
	// If the energy is below the peak energy
	if (energy <= (alpha-beta)*e0)
	{
		return pow(energy/100,alpha) * exp(-energy/e0);
	}
	// If the energy is above the peak energy
	else
	{
		return pow((alpha-beta)*e0/100,alpha-beta) * exp(beta-alpha) * pow(energy/100,beta);
	}
}

// Make the energy axis of a spectrum, bounded by Emin and Emax with a number of bins = num_en_bins
int make_en_axis(float energy_axis[], float Emin, float Emax, int num_en_bins)
{
	// Move to log space to define the energy interval with equally spaced points (in log space)
	float log_Emin = log10(Emin);
	float log_Emax = log10(Emax);
	float log_dE = (log_Emax - log_Emin) / num_en_bins;

	// For each bin, calculate the energy axis
	for(int i=0;i<num_en_bins;i++)
	{
		energy_axis[i] = pow( 10, log_Emin + (i*log_dE));
	}

	return 0;

}