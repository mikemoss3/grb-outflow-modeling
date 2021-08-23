/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-11

Header file for the Spectrum.hpp class

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef SPECCLASS_H
#define SPECCLASS_H

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries

using namespace std;

// Declare Spectrum class
class Spectrum
{
public:
	// Spectrum constructor
	Spectrum(float energ_min=50, float energ_max=300, float num_energ_bins=80, float tmin=0, float tmax=20, float z=1);

	// Spectrum Copy function
	Spectrum(const Spectrum& tmp_spec);

	// Spectrum member variables
	float tmin; // Minimum time of the light curve 
	float tmax; // Maximum time of the light curve 
	float energ_min; // Minimum energy of the spectrum used to create the light curve 
	float energ_max; // Maximum energy of the spectrum used to create the light curve 
	int num_energ_bins; // Number of energy bins within the energy range
	float z; // Redshift of the source 
	
	double spectrum_sum; // Records the summation of the spectrum 
	std::vector<double> spectrum_rate; // Records the spectrum rate in each energy bin
	std::vector<float> energ_lo; // Energy bin lower limits
	std::vector<float> energ_mid; // Center of each energy bin
	std::vector<float> energ_hi; // Energy bin upper limits

	// Spectrum member functions
	// Make energy axis of the spectrum 
	void make_energ_axes(float energ_min, float energ_max, float num_energ_bins, bool logscale);
	void make_energ_axes(bool logscale);
	// Set the spectrum to all zeros
	void ZeroSpectrum();

	// Add thermal component spectrum to total spectrum
	void AddThermComp(std::string file_name, float tmin, float tmax);
	// Calls function to calculate the thermal spectrum rate for each energy bin
	void ReadThermalEmissionData(std::string file_name);
	// Read in emission data for thermal emission
	void MakeThermalSpec(float tmin, float tmax);
	// Calculate the thermal spectrum from the given temperature and flux of the emission
	void CalcThermalContribution(float temp, double flux);
	// Thermal spectrum function form based on a modified Planck function. 
	double ThermalSpec(float energy, float temp, float alpha = 0.4);

	// Add Synchrotron component spectrum to total spectrum
	void AddSynchComp(std::string file_name, float tmin, float tmax);
	// Read in emission data for synchrotron emission
	void ReadSynchEmissionData(std::string file_name);
	// Calls function to calculate the synchrotron spectrum rate for each energy bin
	void MakeSynchSpec(float tmin, float tmax);
	// Calculate the synchrotron spectrum from the synchrotron energy and flux of the emission
	void CalcSynchContribution(float Esyn, double e_diss, float delt);
	// Synchrotron spectrum function form,
	double SynchSpec(float energy, float Esyn, float alpha = -1, float beta = -2.5);

	// Band spectrum function form (Band et. al., 1993)
	double Band(float energy, float e0, float alpha = -1, float beta = -2.5);

	// void AddFluctuations();

	std::vector<float> get_energ_mid();
	// Make the energy axis of a spectrum, bounded by Emin and Emax with a number of bins = num_en_bins
	int make_en_axis(std::vector<float> & energy_axis, float Emin, float Emax, int num_en_bins);
	
	// Write spectrum to text file
	void WriteToTXT(std::string out_file_name);
	// Write spectrum to FITS file
	void WriteToFITS(std::string out_file_name);

// protected:
private:
	int num_emission_events;

	// Initialize arrays to store thermal emission event data 
	std::vector<float> te;
	std::vector<float> ta;
	std::vector<float> delt;
	std::vector<float> temp;
	std::vector<double> flux;
	std::vector<float> rphot;

	// Initialize arrays to store Synchrotron emission event data 
	// std::vector<float> te;
	// std::vector<float> ta;
	std::vector<float> asyn;
	std::vector<float> beq;
	std::vector<float> gammae;
	std::vector<float> esyn;
	std::vector<float> gammar;
	std::vector<double> e_diss;
	// std::vector<float> delt;
	std::vector<float> tau;
	std::vector<float> relvel;

};

#endif 
