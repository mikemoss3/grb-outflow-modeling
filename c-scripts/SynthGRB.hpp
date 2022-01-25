/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-30

Header file for GRB.cpp class

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef SYNTHGRBCLASS_H
#define SYNTHGRBCLASS_H

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "fitsio.h"

// Import Custom Libraries
#include "cosmology.hpp"
#include "utilfuncs.hpp"
#include "ShellDist.hpp"
#include "ModelParams.hpp"
#include "Spectrum.hpp"
#include "LightCurve.hpp"

using namespace std;

// Declare light curve class
class SynthGRB
{
public:
	// SynthGRB constructor
	SynthGRB();
	SynthGRB(float tw, float dte, float eps_e, float eps_b, float zeta, double E_dot_iso, float theta, float r_open, float eps_th, float sigma, float p, std::string LorentzDist, std::string ShellDistParamsFile);
	SynthGRB(ModelParams * input_model_params);

	// GRB member variables 
	// float energ_min; // Lower energy bound of spectrum
	// float energ_max; // Upper energy bound of spectrum 
	// int num_energ_bins; // Number of energy bins in spectrum
	// float tmin; // Lower time bound of light curve
	// float tmax; // Upper time bound of light curve
	// float dt; // Time resolution of light curve 
	// float z; // Redshift of the source
	
	// jet parameters:
	// float tw = 10; // sec, Duration of the wind
	// float dte = 0.002; // sec, Time between successive shell launches
	// float eps_e; //  Fraction of dissipated energy that goes into the electrons 
	// float eps_b; // Fraction of dissipated energy that goes into the magnetic field 
	// float zeta; // Fraction of electrons which are accelerated 
	// double E_dot_iso; // erg/s, Injected isotropic energy rate
	// float theta; // radians, Half-opening angle of the jet
	// float r_open; // cm, Opening radius of the jet
	// float eps_th; // Fraction of energy in the outflow in the form of thermal energy 
	// float sigma; // Magnetization of the outflow 
	// float p; // Power law index of the electron population 
	// std::string LorentzDist; // Distribution of the jet shells
	// std::string ShellDistParamsFile; // File that contains the parameters to create the distribution of jet shells
	ModelParams * p_model_params; // Stores all model parameters	
	Spectrum * p_source_spectrum; // Pointer to the source spectrum
	LightCurve * p_source_light_curve; // Pointer to the source light curve
	ShellDist * p_jet_shells; // Pointer to the jet shells 

	ShellDist jet_shells; // Shell distribution

	double E_dot_kin; // Fraction of energy in Kinetic
	double E_dot_therm; // Fraction of energy in thermal
	
	float gamma_bar; // Average Lorentz factor of all jet shells 
	float m_bar; // Average mass of each shell
	double m_tot; // Average mass of each shell

	// Initialize jet based on current jet parameters and shell distribution
	void InitializeJet();
	// Load jet parameters from a text file
	void LoadJetParamsFromTXT(std::string file_name);

	// GRB member function
	// Set model params
	void set_model_params(ModelParams * new_mod_params);
	// Set jet shells
	void set_jet_shells();

	// Set the source spectrum 
	// This class can make the source spectrum using the given emission data, but it can be manually set here.
	void set_source_spectrum(Spectrum *in_source_spectrum);
	// Set the source light curve 
	// This class can make the source light curve using the given emission data, but it can be manually set here.
	void set_source_light_curve(LightCurve *in_source_light_curve);

	// Reset jet simulation variable vectors
	void reset_simulation();
	// Perform simulation of jet dynamics using the loaded jet parameters
	void SimulateJetDynamics();

	// Make the source spectrum using the emission data 
	void make_source_spectrum(float energ_min = 50., float energ_max = 350., int num_energ_bins = 50, float tmin = 0., float tmax = 30.);
	// Calls function to calculate the thermal spectrum rate for each energy bin
	void MakeThermalSpec(Spectrum * therm_spectrum, float tmin, float tmax);
	// Calculate the thermal spectrum from the given temperature and flux of the emission
	void CalcThermalContribution(Spectrum * therm_spectrum, float temp, double flux);
	// Thermal spectrum function form based on a modified Planck function. 
	double ThermalSpec(float energy, float temp, float alpha=0.4);

	// Calls function to calculate the internal shock spectrum rate for each energy bin
	void MakeISSpec(Spectrum * intsh_spectrum, float tmin, float tmax, float alpha = -0.7, float beta = -3.5);

	// Calls function to calculate the external shock spectrum rate for each energy bin
	void MakeExtShockSpec(Spectrum * extsh_spectrum, float tmin, float tmax);
	// Calls function to calculate the forward shock spectrum rate for each energy bin
	void MakeFSSpec(Spectrum * extsh_spectrum, float tmin, float tmax, float alpha = -0.7, float beta = -3.5);
	// Calls function to calculate the reverse shock spectrum rate for each energy bin
	void MakeRSSpec(Spectrum * extsh_spectrum, float tmin, float tmax, float alpha = -1.5, float beta = -3.5);

	// Calculate the synchrotron spectrum from the synchrotron energy and flux of the emission
	void CalcSynchContribution(Spectrum * synch_spectrum, double Esyn, double e_diss, double delt, float alpha = -0.7, float beta = -3.5);
	// Synchrotron spectrum function form,
	double SynchSpec(float energy, double Esyn, float alpha = -0.7, float beta = -3.5);
	// Make the source light curve using the emission data
	void make_source_light_curve(float energ_min, float energ_max, float Tstart, float Tend, float dt);

	// Load SynthGRB from file
	void LoadFromFile(std::string file_name);

	// Function to write out jet parameters found in the jet simulation
	void write_out_jet_params(std::string dir_path_name);
	// Write spectrum to text file
	void WriteSpectrumToTXT(std::string out_file_name);
	// Write spectrum to FITS file
	void WriteSpectrumToFITS(std::string out_file_name);
	// Write light curve to text file
	void WriteLightCurveToTXT(std::string out_file_name);
	// Write light curve to FITS file
	void WriteLightCurveToFITS(std::string out_file_name);


private:
	// Calculate final Lorentz factor of two colliding shells (after complete redistribution of momenta)
	double _calc_shell_coll_gamma(float s1g, float s2g, float s1m, float s2m);
	// Calculate relative velocity between two shells with Lorentz factors gamma_1 and gamma_2, where gamma_1 > gamma_2
	double _rel_vel(float gamma_1, float gamma_2);
	// Check if shells are ordered
	int _check_if_sorted();
	// Make the energy axis of a spectrum, bounded by Emin and Emax with a number of bins = num_en_bins
	int _make_en_axis(std::vector<float> & energy_axis, float emin, float emax, int num_en_bins);

	// Initialize arrays to store thermal emission event data 
	std::vector<float> te_therm; // sec, Times of emission (in the rest frame of the jet), is also equal to the time a shell crosses the photosphere
	std::vector<float> ta_therm; // sec, Times when the emission arrives at the observer (in the observer frame)
	std::vector<float> delt_therm; // sec, Widths of the emission (in observer frame)
	std::vector<float> T_phot; // K, Temperatures at photosphere
	std::vector<double> L_phot; // erg/s, Luminosities of thermal emission
	std::vector<float> r_phot; // cm, radii of photosphere 

	// Initialize arrays to store Synchrotron emission event data 
	std::vector<double> te_is; // sec, Time of emission (in the rest frame of the jet), is also equal to the time a shell crosses the photosphere
	std::vector<double> ta_is; // sec, Time when the emission arrives at the observer (in the observer frame)
	std::vector<double> delt_is; // sec, Width of the emission (in observer frame)
	std::vector<double> beq_is; // (erg/cm^3)^1/2, Magnetic field density, assuming B = Beq 
	std::vector<float> gamma_e_is; // Lorentz factor of the electron population
	std::vector<double> esyn_is; // erg, Synchrotron energy emitted by accelerated electron
	std::vector<float> gamma_r_is; // Approximation of the combined Lorentz factor of the colliding shells
	std::vector<double> e_diss_is; // erg, Dissipated energy during the collision 
	std::vector<float> asyn_is; // Fraction of the energy in electrons which goes into synchrotron (as opposed to Inverse Compton)
	std::vector<float> tau; // Optical depth at the location of the collision
	std::vector<float> relvel; // Relative velocity between the two shells

	// Initialize arrays to store external shock emission data
	std::vector<float> te_fs; // sec, Time of emission (in the rest frame of the jet)
	std::vector<float> ta_fs; // sec, Time when the emission arrives at the observer (in the observer rest frame)
	std::vector<float> delt_fs; // sec, duration of the emission (in the observer frame)
	std::vector<double> beq_fs; // (erg/cm^3)^1/2, Magnetic field density, assuming B = Beq 
	std::vector<float> gamma_e_fs; // Lorentz factor of the electron population
	std::vector<double> esyn_fs; // erg, Synchrotron energy emitted by the accelerated electrons in the FS
	std::vector<float> gamma_r_fs; // Combined Lorentz factor after collision
	std::vector<double> e_diss_fs; // Energy dissipated in a FS	
	// std::vector<double> eps_fs; // Internal energy dissipated in a collision 
	// std::vector<double> rho_fs; // g cm^-3, Density of the collision region
	
	std::vector<float> te_rs; // sec, Time of emission (in the rest frame of the jet)
	std::vector<float> ta_rs; // sec, Time when the emission arrives at the observer (in the observer rest frame)
	std::vector<float> delt_rs; // sec, duration of the emission (in the observer frame)
	std::vector<double> beq_rs; // (erg/cm^3)^1/2, Magnetic field density, assuming B = Beq 
	std::vector<float> gamma_e_rs; // Lorentz factor of the electron population
	std::vector<double> esyn_rs; // erg, Synchrotron energy emitted by the accelerated electrons in the RS	
	std::vector<float> gamma_r_rs; // Combined Lorentz factor after collision
	std::vector<double> e_diss_rs; // Energy dissipated in a RS
	// std::vector<double> eps_rs; // Internal energy dissipated in a collision 
	// std::vector<double> rho_rs; // g cm^-3, Density of the collision region


};

#endif 