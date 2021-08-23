/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-20

This script can be used to create the spectrum for given emission data.
This script calls methods from spectrum_funcs.cpp to calculate the spectrum of the source for each energy bin.

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
#include "cosmology.hpp"

using namespace std;


// Define Spectrum class and member functions 

// Spectrum constructors
Spectrum::Spectrum(float energ_min, float energ_max, float num_energ_bins, float tmin, float tmax, float z)
{
	this->tmin = tmin; // Minimum time of the light curve 
	this->tmax = tmax; // Maximum time of the light curve 
	
	this->energ_min = energ_min; // Minimum energy of the spectrum used to create the light curve 
	this->energ_max = energ_max; // Maximum energy of the spectrum used to create the light curve 
	this->num_energ_bins = num_energ_bins; // Number of energy bins within the energy range

	this->z = z; // Redshift to the source

	// Fill in the energy vector
	make_energ_axes(energ_min, energ_max, num_energ_bins, true);
}
 
Spectrum::Spectrum(const Spectrum& tmp_spec)
{
	Spectrum(tmp_spec.energ_min, tmp_spec.energ_max, tmp_spec.num_energ_bins, tmp_spec.tmin, tmp_spec.tmax, tmp_spec.z);
}

// Spectrum member functions
// Fill in energy vectors 
void Spectrum::make_energ_axes(float energ_min, float energ_max, float num_energ_bins, bool logscale)
{

	this->energ_min=energ_min;
	this->energ_max=energ_max;
	this->num_energ_bins=num_energ_bins;


	if(logscale == true)
	{
		float logdE = log( energ_max - energ_min ) / num_energ_bins;
		for(int i=0; i< num_energ_bins; ++i)
		{
			energ_lo.push_back( energ_min + exp(logdE*i) );
			energ_hi.push_back( energ_min + exp(logdE*(i+1)) );
			energ_mid.push_back(  energ_min + ( (exp(logdE*i) + exp(logdE*(i+1)) ) /2) );
		}
	}
	else
	{
		float dE = (energ_max - energ_min) / num_energ_bins;
		for(int i=0; i<num_energ_bins; ++i)
		{
			energ_lo.push_back( energ_min + (dE*i) );
			energ_hi.push_back( energ_min + (dE*(i+1)) );
			energ_mid.push_back( energ_min + ( (2*i+1)/2) );
		}		
	}
	std::cout << "Spectrum has been set back to zeros." << "\n";
	ZeroSpectrum();
	spectrum_sum = 0;
}
void Spectrum::make_energ_axes(bool logscale)
{
	make_energ_axes(energ_min, energ_max, num_energ_bins, logscale);
}
// Zero out the spectrum (this also serves to set the length of the vector)
void Spectrum::ZeroSpectrum()
{
	for(int i=0; i < num_energ_bins; i++)
	{
		spectrum_rate.push_back(0.);
	}
	spectrum_sum = 0.;
}

// Add thermal component spectrum to total spectrum
void Spectrum::AddThermComp(std::string file_name, float tmin, float tmax)
{

	// Read in emission event data and store parameters
	ReadThermalEmissionData(file_name);

	// Use the emission event data to create a spectrum within the defined time range
	MakeThermalSpec(tmin, tmax);
}
// Read in emission data for thermal emission
void Spectrum::ReadThermalEmissionData(std::string file_name)
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
			te.push_back( stof(col1_val) );
			ta.push_back( stof(col2_val) );
			delt.push_back( stof(col3_val) );
			temp.push_back( stof(col4_val) );
			flux.push_back( stod(col5_val) );
			rphot.push_back( stof(col6_val) );
		}
		++iter;
	}

	num_emission_events = iter;

	// Close files and free memory 
	file_spec_data.close();
}
// Calls function to calculate the thermal spectrum rate for each energy bin
void Spectrum::MakeThermalSpec(float t_min, float t_max)
{
	// For each emission event, calculate the emitted spectrum
	for(int i=0;i<num_emission_events;++i)
	{		
		// We only want to take the emission that occurs between the specified tmin and tmax. 
		// The emission occurs between (ta+delt), if any of it overlaps with tmin and tmax, calculate its contribution.
		if ( ta.at(i)*(1+z) <= tmax and (ta.at(i)+delt.at(i))*(1+z) >= tmin )
		{
			// Call function to calculate thermal spectrum rate
			CalcThermalContribution(temp.at(i),flux.at(i));
			}
		// Else, if this emission occurred outside of the time interval, don't add it to the spectrum.
	}
}

// Calculate the thermal spectrum from the given temperature and flux of the emission
void Spectrum::CalcThermalContribution(float temp, double flux)
{
	// Calculate the normalization
	float norm=0; // Set normalization to zero
	float e_min_bol = 0.01; // Minimum energy of the Bolometric band
	float e_max_bol = 1e5; // Maximum energy of the Bolometric band
	int norm_num_bin=40*log(e_max_bol/e_min_bol); // Number of energy bins to use to calculate normalization
	std::vector<float> norm_energy_axis(norm_num_bin+1); // Make vector to store energy axis (initialized to zero's) 
	make_en_axis(norm_energy_axis,e_min_bol,e_max_bol,norm_num_bin+1); // Fill in energy values

	float en_curr; // Current energy to evaluate the addition to the normalization
	// For each energy bin along the normalization energy axis, calculate the addition to the normalization and add it. 
	for(int i=0;i < norm_num_bin; ++i)
	{
		en_curr = norm_energy_axis.at(i); // Set the current energy
		// Calculate the contribution according to a Left-Riemann-Sum
		norm += (norm_energy_axis.at(i+1)-en_curr) * en_curr * ThermalSpec(en_curr,temp);
	}

	// Calculate the rate of the thermal spectrum
	double tmp_val=0; // Temporary value to store spectrum contributions
	// For each energy bin along the energy axis
	for(int i=0; i < num_energ_bins; ++i)
	{
		en_curr = energ_mid.at(i); // Set the current energy
		
		// The current temp and energy bin define the count rate, the normalization found above is applied.
		// This must still be multiplied by the flux of the source.
		tmp_val = flux * ThermalSpec(en_curr,temp) / norm;

		// Add the contribution to the total spectrum according to a Left-Riemann-Sum
		spectrum_sum += (energ_hi.at(i) - energ_lo.at(i)) * tmp_val;
		// Check if the spectrum rate per energy bin is requested and store it if so. 
		spectrum_rate[i] += tmp_val;

	}
}
// Thermal spectrum function form based on a modified Planck function. 
double Spectrum::ThermalSpec(float energy, float temp, float alpha)
{
	return pow(energy/(kb_kev*temp),1+alpha)/(exp(energy/(kb_kev*temp))-1);
}

// Add Synchrotron component spectrum to total spectrum
void Spectrum::AddSynchComp(std::string file_name, float tmin, float tmax)
{

	// Read in emission event data and store parameters
	ReadSynchEmissionData(file_name);

	// Use the emission event data to create a spectrum within the defined time range
	MakeSynchSpec(tmin, tmax);
}
// Read in emission data for synchrotron emission
void Spectrum::ReadSynchEmissionData(std::string file_name)
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
			te.push_back( stof(col1_val) );
			ta.push_back( stof(col2_val) );
			asyn.push_back( stof(col3_val) );
			beq.push_back( stof(col4_val) );
			gammae.push_back( stof(col5_val) );
			esyn.push_back( stof(col6_val) );
			gammar.push_back( stof(col7_val) );
			e_diss.push_back( stod(col8_val) );
			delt.push_back( stof(col9_val) );
			tau.push_back( stof(col10_val) );
			relvel.push_back( stod(col11_val) );
		}
		++iter;
	}

	num_emission_events = iter;

	// Close files and free memory 
	file_spec_data.close();
}
// Calls function to calculate the synchrotron spectrum rate for each energy bin
void Spectrum::MakeSynchSpec(float tmin, float tmax)
{
	// For each emission event, calculate the emitted spectrum
	for(int i=0;i<num_emission_events;++i)
	{
		// We only want to take the emission that occurs between the specified Tmin and Tmax. 
		// The emission occurs between (ta+delt), if any of it overlaps with Tmin and Tmax, calculate its contribution.
		if ( ta.at(i)*(1+z) <= tmax and (ta.at(i)+delt.at(i))*(1+z) >= tmin )
		{
			// The emission will only be observable if the relativistic velocity is great than the local sound speed v_s/c = 0.1
			// And if the wind is transparent to the radiation
			if ( relvel.at(i) > 0.1 and tau.at(i) < 1)
			{
				// Call function to calculate thermal spectrum rate
				CalcSynchContribution(esyn.at(i),e_diss.at(i),delt.at(i));	
			}
			
		}
		// Else, if this emission occurred outside of the time interval, don't add it to the spectrum. 
	}
}
// Calculate the synchrotron spectrum from the synchrotron energy and flux of the emission
void Spectrum::CalcSynchContribution(float esyn, double e_diss, float delt)
{
	// Calculate the normalization
	float norm=0; // Set normalization to zero
	float e_min_bol = 0.01; // Minimum energy of the Bolometric band
	float e_max_bol = 1e5; // Maximum energy of the Bolometric band
	int norm_num_bin=20*log(e_max_bol/e_min_bol); // Number of energy bins to use to calculate normalization
	std::vector<float> norm_energy_axis(norm_num_bin+1); // Make vector to store energy axis (initialized to zero's) 
	make_en_axis(norm_energy_axis,e_min_bol,e_max_bol,norm_num_bin+1); // Make energy axis

	float en_curr; // Current energy to evaluate the addition to the normalization
	// For each energy bin along the normalization energy axis, calculate the addition to the normalization and add it. 
	for(int i=0;i < norm_num_bin; ++i)
	{
		en_curr = norm_energy_axis.at(i); // Set the current energy
		// Calculate the contribution according to a Left-Riemann-Sum
		norm += (norm_energy_axis.at(i+1)-en_curr) * en_curr * SynchSpec(en_curr, esyn);
	}

	// Calculate the rate of the thermal spectrum
	double tmp_val=0; // Temporary value to store spectrum contributions
	// For each energy bin along the energy axis
	for(int i=0; i < num_energ_bins; ++i)
	{
		float en_curr = energ_mid.at(i); // Set the current energy

		// The current temp and energy bin define the count rate, the normalization found above is applied.
		// This must still be multiplied by the energy dissipated during the emission event.
		// The energy dissipated can be turned into Flux by dividing the energy dissipated by the observed emission duration (delt).
		tmp_val = e_diss * SynchSpec(en_curr, esyn) / norm / delt;

		// Add the contribution to the total spectrum according to a Left-Riemann-Sum
		spectrum_sum += (energ_hi.at(i) - energ_lo.at(i)) * tmp_val;

		// Check if the spectrum rate per energy bin is requested and store it if so. 
		spectrum_rate.at(i) += tmp_val;
	}

}
// Synchrotron spectrum function form,
// Defaults alpha = -1, beta = -2.5
double Spectrum::SynchSpec(float energy, float esyn, float alpha, float beta)
{
	return (1/esyn) * Band(energy, esyn, alpha, beta);
}

// Band spectrum function form (Band et. al., 1993)
// Defaults alpha = -1, beta = -2.5
double Spectrum::Band(float energy, float e0, float alpha, float beta)
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

std::vector<float> Spectrum::get_energ_mid()
{
	return energ_mid;
}

// Make the energy axis of a spectrum, bounded by Emin and Emax with a number of bins = num_en_bins
int Spectrum::make_en_axis(std::vector<float> & energy_axis, float emin, float emax, int num_en_bins)
{
	// Move to log space to define the energy interval with equally spaced points (in log space)
	float log_emin = log10(emin);
	float log_emax = log10(emax);
	float log_de = (log_emax - log_emin) / num_en_bins;

	// For each bin, calculate the energy axis
	for(int i=0;i<num_en_bins;++i)
	{
		energy_axis.at(i) = pow( 10, log_emin + (i*log_de));
	}

	return 0;
}

void Spectrum::WriteToTXT(std::string out_file_name)
{
	// Write the light curve to a text file
	ofstream spec_fule; // Construct file 
	spec_fule.open(out_file_name); // Open text file with this name
	int i=0;
	// For each time bin, write the time and count rate to the file.
	while ( i < num_energ_bins)
	{
		spec_fule << energ_mid.at(i);
		spec_fule << " ";
		spec_fule << spectrum_rate.at(i);
		spec_fule << "\n";		
		i++;
	}
	spec_fule.close(); // Close file
}

void Spectrum::WriteToFITS(std::string out_file_name)
{
}





