/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-20

Header file for light_curve_funcs.cpp

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef SPECTRUMFUNCS_H
#define SPECTRUMFUNCS_H

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

// Declaring function
// Calls function to calculate the thermal spectrum rate for each energy bin
void read_in_thermal_emission_data(std::string file_name, float te[], float ta[], float delt[], float Temp[], double Flux[], float Rphot[]);
// Read in emission data for thermal emission
void make_thermal_spec(double spectrum_sum[], double spectrum_dE[], float ta[], float Temp[], double Flux[], float delt[], int num_lines, float Tmin, float Tmax, float z = 1, float Emin = 0.1, float Emax = 10000, int num_en_bins = 50, bool sum_only = false , bool nuFnu = true);
// Calculate the thermal spectrum from the given temperature and flux of the emission
void calc_thermal_spec(double &spectrum_sum, double spectrum_dE[], float temp, double flux, float z = 1, float Emin = 0.1, float Emax = 10000, int num_en_bins = 50, bool sum_only = false , bool nuFnu = true);
// Thermal spectrum function form based on a modified Planck function. 
double thermal_spectrum(float energy, float temp, float alpha = 0.4);

// Calls function to calculate the synchrotron spectrum rate for each energy bin
void read_in_synch_emission_data(std::string file_name, float te[], float ta[], float asyn[], float Beq[], float gammae[], float Esyn[], float gammar[], double e_diss[],float delt[], float tau[], float relvel[]);
// Read in emission data for synchrotron emission
void make_synch_spec(double spectrum_sum[], double spectrum_dE[], float ta[], float Esyn[], double e_diss[], float delt[], float tau[], float relvel[], int num_lines, float Tmin, float Tmax, float z = 1, float Emin = 0.1, float Emax = 10000, int num_en_bins = 50, bool sum_only = false , bool nuFnu = true);
// Calculate the synchrotron spectrum from the synchrotron energy and flux of the emission
void calc_synch_spec(double &spectrum_sum, double spectrum_dE[], float Esyn, double e_diss, float delt, float z = 1, float Emin = 0.1, float Emax = 10000, int num_en_bins = 50, bool sum_only = false , bool nuFnu = true);
// Synchrotron spectrum function form,
double synch_spectrum(float energy, float Esyn, float alpha = -1, float beta = -2.5);

// Band spectrum function form (Band et. al., 1993)
double band(float energy, float e0, float alpha = -1, float beta = -2.5);
// Make the energy axis of a spectrum, bounded by Emin and Emax with a number of bins = num_en_bins
int make_en_axis(float energy_axis[], float Emin, float Emax, int num_en_bins);

#endif 