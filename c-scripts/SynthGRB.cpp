/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-30

SynthGRB class which holds and interfaces LightCurve, Spectrum, and Response objects.
*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <sstream>
#include "fitsio.h"

// Import Custom Libraries
#include "cosmology.hpp"
#include "utilfuncs.hpp"
#include "ShellDist.hpp"
#include "ModelParams.hpp"
#include "Spectrum.hpp"
#include "LightCurve.hpp"
#include "SynthGRB.hpp"

using namespace std;

// SynthGRB constructor, uses default arguments
SynthGRB::SynthGRB()
{
	// Assign pointers to NULL
	p_source_spectrum = NULL;
	p_source_light_curve = NULL;
	p_jet_shells = NULL;
	p_model_params = new ModelParams();

	// Initialize jet based on current jet parameters and shell distribution
	InitializeJet();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Overloaded SynthGRB constructor, uses given arguments 
SynthGRB::SynthGRB(float tw, float dte, 
		double E_dot_iso, float theta, float r_open, float eps_th, float sigma,
		float eps_e_int, float eps_b_int, float zeta_int, float p_int,
		float eps_e_ext, float eps_b_ext, float zeta_ext, float p_ext,
		float k_med, double rho_not, 
		std::string LorentzDist, std::string ShellDistParamsFile)
{
	// Assign pointers to NULL
	p_source_spectrum = NULL;
	p_source_light_curve = NULL;
	p_jet_shells = NULL;
	p_model_params = new ModelParams(tw, dte, 
		E_dot_iso, theta, r_open, eps_th, sigma,
		eps_e_int, eps_b_int, zeta_int, p_int,
		eps_e_ext, eps_b_ext, zeta_ext, p_ext,
		k_med, rho_not,
		LorentzDist, ShellDistParamsFile);

	// Initialize jet based on current jet parameters and shell distribution
	InitializeJet();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Overloaded SynthGRB constructor, uses given arguments 
SynthGRB::SynthGRB(ModelParams * input_model_params)
{
	// Assign pointers to NULL
	p_source_spectrum = NULL;
	p_source_light_curve = NULL;
	p_jet_shells = NULL;
	p_model_params = input_model_params;

	// Initialize jet based on current jet parameters and shell distribution
	InitializeJet();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Initialize jet based on current jet parameters and shell distribution
void SynthGRB::InitializeJet()
{
	E_dot_kin = (*p_model_params).E_dot_iso/(1+(*p_model_params).sigma); // Fraction of energy in Kinetic
	E_dot_therm = (*p_model_params).E_dot_iso*(*p_model_params).eps_th; // Fraction of energy in thermal

	// Make shell distribution
	set_jet_shells();

	// Average Lorentz factor of all jet shells  
	gamma_bar = 0.0;
	for(int i=0; i<(*p_model_params).numshells; ++i) 
	{
	     gamma_bar += (*p_jet_shells).shell_gamma.at(i);
	}
	gamma_bar /= (*p_model_params).numshells;
	// Average mass of each shell
	m_bar = (*p_model_params).E_dot_iso*(*p_model_params).tw/gamma_bar/pow(c_cm,2.)/(*p_model_params).numshells;
	m_tot = (*p_model_params).E_dot_iso*(*p_model_params).tw/gamma_bar/pow(c_cm,2.);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load jet parameters from a file
void SynthGRB::LoadJetParamsFromTXT(std::string file_name)
{
	// Array to store params 
	string inputs[19];
	int i = 0;

	// Load in the file
	ifstream file_jet_params(file_name);
	string line_jet_params;
	if ( file_jet_params.is_open() ) 
	{ 
		while(getline(file_jet_params, line_jet_params))
		{
			// If it is a commented line, just pass by
			if(line_jet_params[0] == '#'){ continue; }
			// Else: 
			inputs[i] = line_jet_params;
			++i;
		}	
	
		// Close files and free memory 
		file_jet_params.close();
	}
	else{ std::cout << "Unable to open file."; }

	// Set model parameters 
	// tw
	// dte

	// E_dot_iso
	// theta
	// r_open
	// eps_th
	// sigma
	
	// eps_e_int
	// eps_b_int
	// zeta_int
	// p_int
	
	// eps_e_ext
	// eps_b_ext
	// zeta_ext
	// p_ext
	
	// k_med
	// rho_not
	
	// LorentzDist
	// ShellDistParamsFile	
	p_model_params = new ModelParams(
		stof(inputs[0]), // tw ,sec
		stof(inputs[1]), // dte, sec
		stod(inputs[2]), // E_dot_iso, erg/s
		stof(inputs[3]), // theta, rad
		stof(inputs[4]), // r_open, cm
		stof(inputs[5]), // eps_th
		stod(inputs[6]), // sigma
		stof(inputs[7]), // eps_e_int
		stof(inputs[8]), // eps_b_int
		stof(inputs[9]), // zeta_int
		stof(inputs[10]), // p_int
		stof(inputs[11]), // eps_e_ext
		stof(inputs[12]), // eps_b_ext
		stof(inputs[13]), // zeta_ext
		stof(inputs[14]), // p_ext
		stof(inputs[15]), // k_med
		stod(inputs[16]), // rho_not
		inputs[17], // Lorentz dist.
		inputs[18]); // Lorentz dist. parameter file

	// Initialize jet based on current jet parameters and shell distribution
	InitializeJet();

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Set model parameters
void SynthGRB::set_model_params(ModelParams * new_mod_params)
{
	this->p_model_params = new_mod_params;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Set jet shells 
void SynthGRB::set_jet_shells()
{
	if(p_jet_shells != NULL)
	{
		delete p_jet_shells;
	}

	p_jet_shells = new ShellDist((*p_model_params).numshells, (*p_model_params).E_dot_iso);

	if((*p_model_params).LorentzDist == "step" or (*p_model_params).LorentzDist == "Step" )
	{
		// Load shell distribution parameters from (*p_model_params).ShellDistParamsFile
		if((*p_model_params).ShellDistParamsFile == "Default")
		{
			// Make shell distribution with default arguments
			(*p_jet_shells).step((*p_model_params).dte);
		}
		else
		{
			// Array to store para-ms 
			string inputs[4];
			int i = 0;

			// Load from designated file
			ifstream file_jet_params((*p_model_params).ShellDistParamsFile);
			string line_jet_params;
			if ( file_jet_params.is_open() ) 
			{
				while(getline(file_jet_params, line_jet_params))
				{
					// If it is a commented line, just pass by
					if(line_jet_params[0] == '#'){ continue; }
					// Else: 
					inputs[i] = line_jet_params;
					++i;
				}	

				// Close files and free memory 
				file_jet_params.close();
			}
			else std::cout << "Unable to open file.";
		
			// Make shell distribution with input parameters
			bool fluc;
			istringstream(inputs[3]) >> fluc; 
			(*p_jet_shells).step((*p_model_params).dte, stof(inputs[0]), stof(inputs[1]), stof(inputs[2]), fluc);
		}
		
	}
	if((*p_model_params).LorentzDist == "smoothstep" or (*p_model_params).LorentzDist == "Smoothstep" )
	{
		// Load shell distribution parameters from (*p_model_params).ShellDistParamsFile
		if((*p_model_params).ShellDistParamsFile == "Default")
		{
			// Make shell distribution with default arguments
			(*p_jet_shells).smoothstep((*p_model_params).dte);
		}
		else
		{
			// Array to store params 
			string inputs[4];
			int i = 0;

			// Load from designated file
			ifstream file_jet_params((*p_model_params).ShellDistParamsFile);
			string line_jet_params;
			if ( file_jet_params.is_open() ) 
			{

				while(getline(file_jet_params, line_jet_params))
				{
					// If it is a commented line, just pass by
					if(line_jet_params[0] == '#'){ continue; }
					// Else: 
					inputs[i] = line_jet_params;
					++i;
				}
				
				// Close files and free memory 
				file_jet_params.close();
			}
			else std::cout << "Unable to open file.";

			// Make shell distribution with input parameters
			bool fluc;
			istringstream(inputs[3]) >> fluc; 
			(*p_jet_shells).smoothstep((*p_model_params).dte, stof(inputs[0]), stof(inputs[1]), stof(inputs[2]), fluc);

		}
		
	}
	if((*p_model_params).LorentzDist == "oscillatory" or (*p_model_params).LorentzDist == "Oscillatory" or (*p_model_params).LorentzDist == "osci")
	{
		// Load shell distribution parameters from (*p_model_params).ShellDistParamsFile
		if((*p_model_params).ShellDistParamsFile == "Default")
		{
			// Make shell distribution with default arguments
			(* p_jet_shells).oscillatory((*p_model_params).dte);	
		}
		else
		{
			// Array to store params 
			string inputs[5];
			int i = 0;

			// Load from designated file
			ifstream file_jet_params((*p_model_params).ShellDistParamsFile);
			string line_jet_params;
			if ( file_jet_params.is_open() ) 
			{
				while(getline(file_jet_params, line_jet_params))
				{
					// If it is a commented line, just pass by
					if(line_jet_params[0] == '#'){ continue; }
					// Else: 
					inputs[i] = line_jet_params;
					++i;
				}

				// Close files and free memory 
				file_jet_params.close();
			}
			else std::cout << "Unable to open file.";

			// Make shell distribution with input parameters
			bool fluc;
			istringstream(inputs[4]) >> fluc; 
			(*p_jet_shells).oscillatory((*p_model_params).dte, stof(inputs[0]), stof(inputs[1]), stof(inputs[2]), stof(inputs[3]), fluc);
		}
	}
	if((*p_model_params).LorentzDist == "gauss_inject")
	{
		// Load shell distribution parameters from (*p_model_params).ShellDistParamsFile
		if((*p_model_params).ShellDistParamsFile == "Default")
		{
			// Make shell distribution with default arguments
			(* p_jet_shells).gauss_inject((*p_model_params).dte);	
		}
		else
		{
			// Array to store params 
			string inputs[8];
			std::vector<std::vector<float>> vec_vec(4);
			int i = 0;
			int k = 0;

			// Load from designated file
			ifstream file_jet_params((*p_model_params).ShellDistParamsFile);
			string line_jet_params;
			if ( file_jet_params.is_open() ) 
			{
				while(getline(file_jet_params, line_jet_params))
				{
					// If it is a commented line, just pass by
					if(line_jet_params[0] == '#')
					{ continue; }
					else if(i < 3 )
					{
						inputs[i] = line_jet_params;
						++i;
					}
					else if( i < 7)
					{
						std::vector<float> tmp_v;
						// Build an istream that holds the input string
						std::istringstream iss(line_jet_params);


						// Iterate over the istream, using >> to grab floats and push_back to store them in the vector
						std::copy(std::istream_iterator<float>(iss),std::istream_iterator<float>(), std::back_inserter(tmp_v));
						vec_vec.at(k) = tmp_v;

						++i;
						++k;
					}
					else
					{
						inputs[i] = line_jet_params;
						++i;
					}

				}

				// Close files and free memory 
				file_jet_params.close();
			}
			else std::cout << "Unable to open file.";

			// Make shell distribution with input parameters
			bool fluc;
			istringstream(inputs[7]) >> fluc; 
			(*p_jet_shells).gauss_inject((*p_model_params).dte, stof(inputs[0]), stof(inputs[1]), stof(inputs[2]), vec_vec.at(0), vec_vec.at(1), vec_vec.at(2), vec_vec.at(3), fluc);
		}
	}
	if((*p_model_params).LorentzDist == "square_inject")
	{
		// Load shell distribution parameters from (*p_model_params).ShellDistParamsFile
		if((*p_model_params).ShellDistParamsFile == "Default")
		{
			// Make shell distribution with default arguments
			(* p_jet_shells).square_inject((*p_model_params).dte);	
		}
		else
		{
			// Array to store params 
			string inputs[7];
			std::vector<std::vector<float>> vec_vec(3);
			int i = 0;
			int k = 0;

			// Load from designated file
			ifstream file_jet_params((*p_model_params).ShellDistParamsFile);
			string line_jet_params;
			if ( file_jet_params.is_open() ) 
			{
				while(getline(file_jet_params, line_jet_params))
				{
					// If it is a commented line, just pass by
					if(line_jet_params[0] == '#')
					{ continue; }
					else if(i < 3 )
					{
						inputs[i] = line_jet_params;
						++i;
					}
					else if( i < 6)
					{
						std::vector<float> tmp_v;
						// Build an istream that holds the input string
						std::istringstream iss(line_jet_params);


						// Iterate over the istream, using >> to grab floats and push_back to store them in the vector
						std::copy(std::istream_iterator<float>(iss),std::istream_iterator<float>(), std::back_inserter(tmp_v));
						vec_vec.at(k) = tmp_v;

						++i;
						++k;
					}
					else
					{
						inputs[i] = line_jet_params;
						++i;
					}

				}

				// Close files and free memory 
				file_jet_params.close();
			}
			else std::cout << "Unable to open file.";

			// Make shell distribution with input parameters
			bool fluc;
			istringstream(inputs[6]) >> fluc; 
			(*p_jet_shells).square_inject((*p_model_params).dte, stof(inputs[0]), stof(inputs[1]), stof(inputs[2]), vec_vec.at(0), vec_vec.at(1), vec_vec.at(2), fluc);
		}
	}
	if((*p_model_params).LorentzDist == "fred_inject")
	{
		// Load shell distribution parameters from (*p_model_params).ShellDistParamsFile
		if((*p_model_params).ShellDistParamsFile == "Default")
		{
			// Make shell distribution with default arguments
			(* p_jet_shells).fred_inject((*p_model_params).dte);	
		}
		else
		{
			// Array to store params 
			string inputs[8];
			std::vector<std::vector<float>> vec_vec(4);
			int i = 0;
			int k = 0;

			// Load from designated file
			ifstream file_jet_params((*p_model_params).ShellDistParamsFile);
			string line_jet_params;
			if ( file_jet_params.is_open() ) 
			{
				while(getline(file_jet_params, line_jet_params))
				{
					// If it is a commented line, just pass by
					if(line_jet_params[0] == '#')
					{ continue; }
					else if(i < 3 )
					{
						inputs[i] = line_jet_params;
						++i;
					}
					else if( i < 7)
					{
						std::vector<float> tmp_v;
						// Build an istream that holds the input string
						std::istringstream iss(line_jet_params);


						// Iterate over the istream, using >> to grab floats and push_back to store them in the vector
						std::copy(std::istream_iterator<float>(iss),std::istream_iterator<float>(), std::back_inserter(tmp_v));
						vec_vec.at(k) = tmp_v;

						++i;
						++k;
					}
					else
					{
						inputs[i] = line_jet_params;
						++i;
					}

				}

				// Close files and free memory 
				file_jet_params.close();
			}
			else std::cout << "Unable to open file.";

			// Make shell distribution with input parameters
			bool fluc;
			istringstream(inputs[7]) >> fluc; 			
			(*p_jet_shells).fred_inject((*p_model_params).dte, stof(inputs[0]), stof(inputs[1]), stof(inputs[2]), vec_vec.at(0), vec_vec.at(1), vec_vec.at(2), vec_vec.at(3), fluc);
		}
	}
	if((*p_model_params).LorentzDist == "linear")
	{
		// Load shell distribution parameters from (*p_model_params).ShellDistParamsFile
		if((*p_model_params).ShellDistParamsFile == "Default")
		{
			// Make shell distribution with default arguments
			(* p_jet_shells).oscillatory((*p_model_params).dte);	
		}
		else
		{
			// Array to store params 
			string inputs[2];
			int i = 0;

			// Load from designated file
			ifstream file_jet_params((*p_model_params).ShellDistParamsFile);
			string line_jet_params;
			if ( file_jet_params.is_open() ) 
			{
				while(getline(file_jet_params, line_jet_params))
				{
					// If it is a commented line, just pass by
					if(line_jet_params[0] == '#'){ continue; }
					// Else: 
					inputs[i] = line_jet_params;
					++i;
				}

				// Close files and free memory 
				file_jet_params.close();
			}
			else std::cout << "Unable to open file.";

			// Make shell distribution with input parameters
			bool fluc;
			istringstream(inputs[2]) >> fluc; 
			(*p_jet_shells).linear((*p_model_params).dte, stof(inputs[0]), stof(inputs[1]), fluc);
		}
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Set the source spectrum 
// This class can make the source spectrum using the given emission data, but it can be manually set here.
void SynthGRB::set_source_spectrum(Spectrum *in_source_spectrum)
{
	p_source_spectrum = in_source_spectrum;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Set the source light curve 
// This class can make the source light curve using the given emission data, but it can be manually set here.
void SynthGRB::set_source_light_curve(LightCurve *in_source_light_curve)
{
	p_source_light_curve = in_source_light_curve;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Reset jet simulation variable vectors
void SynthGRB::reset_simulation()
{

	set_jet_shells();

	// Reset thermal variable arrays
	te_therm.resize(0);
	ta_therm.resize(0);
	delt_therm.resize(0);
	T_phot.resize(0);
	L_phot.resize(0);
	r_phot.resize(0);
	shell_ind_th.resize(0);

	// Reset IS variable arrays 
	te_is.resize(0);
	ta_is.resize(0);
	delt_is.resize(0);
	beq_is.resize(0);
	gamma_e_is.resize(0);
	esyn_is.resize(0);
	gamma_r_is.resize(0);
	e_diss_is.resize(0);
	nu_c_is.resize(0);
	nu_m_is.resize(0);
	shell_ind_is.resize(0);
	asyn_is.resize(0);
	tau.resize(0);
	relvel.resize(0);
	// eps_star_is.resize(0);

	// Reset FS variable arrays
	te_fs.resize(0);
	ta_fs.resize(0);
	delt_fs.resize(0);
	beq_fs.resize(0);
	gamma_e_fs.resize(0);
	esyn_fs.resize(0);
	gamma_r_fs.resize(0);
	e_diss_fs.resize(0);
	nu_c_fs.resize(0);
	nu_m_fs.resize(0);

	rad_coll_fs.resize(0);
	rho_fs.resize(0);
	eps_star_fs.resize(0);
	num_swept_e_fs.resize(0);
	theta_fs.resize(0);
	
	// Reset RS variable arrays
	te_rs.resize(0);
	ta_rs.resize(0);
	delt_rs.resize(0);
	beq_rs.resize(0);
	gamma_e_rs.resize(0);
	esyn_rs.resize(0);
	gamma_r_rs.resize(0);
	e_diss_rs.resize(0);
	nu_c_rs.resize(0);
	nu_m_rs.resize(0);
	shell_ind_rs.resize(0);
	// eps_star_rs.resize(0);
	// rho_rs.resize(0);

	mass_te.resize(0);
	mass_ta.resize(0);
	mass_rs.resize(0);
	mass_fs.resize(0);
	mass_ej.resize(0);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Perform simulation of jet dynamics using the loaded jet parameters
void SynthGRB::SimulateJetDynamics()
{
	// Reset physical variable vectors
	reset_simulation();

	/* Calculate Thermal considerations */
	// Thermal considerations do not depend on shell collisions, but only on the energy in a particular shell as it passes the photosphere
	float Phi=0;
	float T0=0;
	te_therm.resize((*p_model_params).numshells);
	ta_therm.resize((*p_model_params).numshells);
	delt_therm.resize((*p_model_params).numshells);
	T_phot.resize((*p_model_params).numshells);
	L_phot.resize((*p_model_params).numshells);
	r_phot.resize((*p_model_params).numshells);
	shell_ind_th.resize((*p_model_params).numshells);

	for(int i=0; i<(*p_model_params).numshells; ++i)
	{
		// If the shell is active
		if((*p_jet_shells).shell_status.at(i) == 1)
		{
			// Calculate the photospheric radius for each jet shell, Equation 9 of Hascoet 2013
			r_phot.at(i) = (kappa_T*(*p_model_params).E_dot_iso / ( (1.+(*p_model_params).sigma)*8.*M_PI*pow(c_cm,4.) *pow( (*p_jet_shells).shell_gamma.at(i) ,3.) ) ); // units of light seconds, cm/speed of light
			// r_phot = 2.9*pow(10,13) *((*p_model_params).E_dot_iso/1e53) / ( cc.c*(1+(*p_model_params).sigma)* pow((*p_jet_shells).shell_gamma/100,3) ) # units of light seconds, cm/(3 * 10^8 cm/s)

			// Times when each shell will cross the photosphere
			te_therm.at(i) = (r_phot.at(i) - (*p_jet_shells).shell_radius.at(i) ) / beta((*p_jet_shells).shell_gamma.at(i)); // sec
			ta_therm.at(i) = (te_therm.at(i) - r_phot.at(i)); // sec, Time emission arrives at observer 
			delt_therm.at(i) = r_phot.at(i)/2./pow((*p_jet_shells).shell_gamma.at(i),2.); // sec, Duration width of the emission 

			// Calculate useful constant for next calculations
			Phi = pow((*p_model_params).theta,-2./3.) * pow(r_phot.at(i)*c_cm,-2./3.) * pow((*p_model_params).r_open,2./3.) * pow((*p_jet_shells).shell_gamma.at(i),2./3.); // Eq. 11 in Hascoet 2013

			// Temperature at photosphere
			T0 = pow(E_dot_therm*pow((*p_model_params).theta,2.) / (4.* M_PI * a * c_cm * pow((*p_model_params).r_open,2.) ),1./4.); // K, rearrangement of Eq. 1 in Hascoet 2013 
			// T0 =  (2/3)*np.pow((*p_model_params).eps_th,1/4)*np.pow((*p_model_params).theta/0.1,1/2)*np.pow((*p_model_params).E_dot_iso/1e53,1/4)*np.pow((*p_model_params).r_open/1e7,-1/2) / (cc.kb_kev/1000) # K, alternative expression, Eq. in Hascoet

			T_phot.at(i) = T0*Phi; // K, Equation 7 in Hascoet 2013, observed temperature ( not corrected for 1+z )
			T_phot.at(i) /= 10.; // This makes spectrum look nice, but idk why its here. 

			// Luminosity at photosphere, from equation 8 in Hascoet 2013
			L_phot.at(i) = (pow((*p_model_params).theta,2.) / 4.) * E_dot_therm * Phi; // erg/s, beamed
			// L_phot.at(i) = pow((*p_jet_shells).shell_gamma.at(i),2.) * a * pow( T0*Phi/(*p_jet_shells).shell_gamma.at(i),4.) * c_cm * (M_PI * pow((*p_model_params).theta,2.) * pow(c_cm * r_phot.at(i),2.) );// erg/s, beamed (alternative expression)
			// L_phot.at(i) = E_dot_therm * Phi; // erg/s, isotropic
			// L_phot.at(i) = E_dot_therm * Phi * (*p_model_params).tw / (*p_model_params).numshells; // erg/s, isotropic


			// Record which shell is producing this thermal emission
			shell_ind_th.at(i) = i;
		}
	}

	/* Calculate Internal and External shock considerations */

	// Declare and initialize flags which indicate whether the jet propagation should continue
	bool ord_lorentz = false; // Indicates whether the internal shells are ordered by their Lorentz factor
	bool fs_relativ = true; // Indicates whether the forward shock is still relativistic.
	bool fs_active = false; // Indicates whether the FS is active or not
	bool rs_active = false; // Indicates whether the RS is active or not
	int num_shocks = 1; // Current number of shocks that have occurred 
	float num_fs_shocks = 0.; // Current number of forward shocks that have occurred times q, when q*num_fs_shocks = 50, increment num_shocks
	bool const_gamma_e = false; // Is the fraction of accelerated electrons proportional to the energy dissipated in a shock? (e.g., is gamma_e constant?)
	bool flag_lat_exp = false; // Indicates if lateral expansion should be used to calculate the opening angle of the forward shell

	// Declare and initialize variables for 
	// Internal Shock processes
	
	int tmp_ind_ds_shell = 0; // Temporary index of the down stream shell
	int tmp_ind_us_shell = 0; // Temporary index of the upstream shell
	double tmp_ds_r = 0.; // Temporary Radius of the down stream shell
	double tmp_ds_b = 0.; // Temporary Beta (i.e., v/c) of the down stream shell
	double tmp_us_r = 0.; // Temporary Radius of the up stream shell
	double tmp_us_b = 0.; // Temporary Beta of the upstream shell

	int ind_s_ds = 0; // Index of the down stream shell which will collide next
	int ind_s_us = 0; // Index of the up stream shell which will collide next
	double shell_ds_m = 0.; // g/M_tot, Downstream shell mass 
	double shell_ds_g = 1.; // Downstream shell Lorentz factor 
	double shell_us_m = 0.; // g/M_tot, Upstream shell mass 
	double shell_us_g = 1.; // Upstream shell Lorentz factor 

	double tmp_t_IS = 0.; // sec, Time until collision of the upstream and downstream shells
	double t_IS_lowest = 0.; // sec, Minimum time until next shell collision	

	double tmp_tau = 0.; // Temporary value of the optical depth at the site of the collision

	double min_mass = 0.; // g/M_tot, mass of the less massive shell involved in the internal shock
	// double gamma_comb = 1.; // Final Lorentz factor of the combined shell (after complete redistribution of momenta)
	double rad_coll = 0.; // light seconds, Radius of the collision

	// Declare and initialize variables for 
	//External Shock processes

	double tmp_ej_m = 0.; // g/M, Temporary mass of the ejecta shell 
	float tmp_ej_g = 1.; // Temporary Lorentz factor of the ejecta shell 
	double tmp_fs_m = 0.; // g/M, Temporary mass of the forward shock shell
	float tmp_fs_g = 1.; // Temporary Lorentz factor of the forward shock shell
	double tmp_fs_r = 0.; // Temporary radius of the forward shock shell
	double tmp_rs_m = 0.; // g/M, Temporary mass of the reverse shock shell
	float tmp_rs_g = 1.; // Temporary Lorentz factor of the reverse shock shell
	double tmp_rs_r = 0.; // Temporary radius of the reverse shock shell

	double t_FS_sweep = 0.; // sec, Time until a mass m_ex is swept up by the FS   
	double t_RS_coll = 0.; // sec, Time until a RS collision

	double tmp_r_before = 0.; // light seconds, Temporary radius before the time step

	int fs_shell_index = 0; // To keep track of which shell is representing the FS
	int rs_shell_index = 0;// To keep track of which shell is representing the RS

	// double is_gamma_int = 1.; // Lorentz factor for internal motion in shocked material
	double fs_gamma_int = 1.; // Internal Lorentz factor for the forward shock
	// double rs_gamma_int = 1.; // Internal Lorentz factor for the forward shock

	float q = 0.01; // Ratio which dictates how much external mass needs to be swept up in order for an internal show to occur
	double m_ex = q * 2 * (*p_jet_shells).shell_mass.at(0) / gamma_bar; // g/M, A value to parameterize when a FS occurs
	double m_swept = 0.; // g/M, Mass currently swept up in the FS
	double m_swept_tot = 0.; // g/M, Mass remaining until the next FS 
	double m_remaining = m_ex; // g/M, Mass remaining until the next FS 

	// Shared variables between the processes
	double tmp_te = 0.; // sec, Temporary value of the emission time
	double tmp_asyn = 0.; // Temporary value of the fraction of energy going to synchrotron electrons
	double tmp_beq = 0.; // Temporary value of the magnetic field
	double tmp_gamma_e = 1.; // Temporary value of the typical Lorentz factor of accelerated electrons 
	double tmp_esyn_kev = 0.; // keV, Temporary value of the typical energy of accelerated electrons
	double tmp_esyn_erg = 0.; // erg, 
	double tmp_gamma_r = 1.; // Temporary value of the combined bulk Lorentz factor two colliding shells  
	double tmp_e_diss = 0.; // erg Temporary value of the energy dissipated during the collision
	double tmp_delt = 0.; // sec, Temporary value of the duration of the current emission event
	double tmp_lum_diss = 0.; // erg/s Temporary value of the energy dissipated during the collision normalized by the duration of the emission
	double tmp_nu_c = 0.; // Hz, Frequency associated with the critical synchrotron Lorentz factor
	double tmp_nu_m = 0.; // Hz, Frequency associated with the minimum Lorentz factor of the accelerated electron population
	double tmp_eps_star = 0.; // erg / g, Average proton factor from the collision of two shells
	double tmp_rho = 0.; // g/cm^3, Comoving proton density
	float tmp_theta = (*p_model_params).theta; // rad, opening angle of the jet, it is initially set to the initial opening angle of the jet
	bool tmp_eff = false; // Is the emission event efficient?
	
	double t_syn = 0.; // sec, Synchrotron time-scale
	double quadr_const = 0.; // Constant in the quadratic equation
	double w = 0.; // Indicates whether we are in the Klein-Nishina regime or Thompson (i.e., w>>1 or w<<1, respectively)
	double Q_IC = 0.; // Inverse Compton parameter
	double alpha_ic = 0.; // Fraction of energy that goes into Inverse Compton electrons

	while( (ord_lorentz == false) | (fs_relativ == true) | (rs_active == true) )
	{
		// Set the time until the next internal, forward, or reverse shock as a high number
		t_IS_lowest = std::numeric_limits<double>::max();
		t_FS_sweep = std::numeric_limits<double>::max();
		t_RS_coll = std::numeric_limits<double>::max();

		// Make list of shells which are currently active:
		std::vector<int> active_inds;
		for(int i=0; i<(*p_model_params).numshells; ++i)
		{
			if((*p_jet_shells).shell_status.at(i) == 1)
			{
				active_inds.push_back(i);
			}
		}

		/* Calculate the time of collision between all adjacent (and active) shells */
		// For each active internal shell, calculate the time it will collide with the shell in front of it
		// Then select the shortest time until the next collision
		if(ord_lorentz == false)
		{
			for(size_t i=0; i<active_inds.size()-1; ++i)
			{
				// The current active shell is the down stream shell
				tmp_ind_ds_shell = active_inds.at(i);
				// The next active shell is the up stream shell
				tmp_ind_us_shell = active_inds.at(i+1);

				// The down stream shell must have a Lorentz factor smaller than the upstream shell or they will never collide
				if( (*p_jet_shells).shell_gamma.at(tmp_ind_ds_shell) < (*p_jet_shells).shell_gamma.at(tmp_ind_us_shell) )
				{
					// Both shells must be launched
					// The shells are launched at a radius of zero, before this they are at "negative" radii
					if( ((*p_jet_shells).shell_radius.at(tmp_ind_ds_shell) > 0) & ((*p_jet_shells).shell_radius.at(tmp_ind_us_shell) > 0) )
					{
						// The down stream shell is farther out, but slower
						// The up stream shell is closer to the central engine, but faster
						tmp_ds_r = (*p_jet_shells).shell_radius.at(tmp_ind_ds_shell);
						tmp_ds_b = beta( (*p_jet_shells).shell_gamma.at(tmp_ind_ds_shell) );
						tmp_us_r = (*p_jet_shells).shell_radius.at(tmp_ind_us_shell);
						tmp_us_b = beta( (*p_jet_shells).shell_gamma.at(tmp_ind_us_shell) );

						// Calculate the time until these shells collide
						tmp_t_IS = (tmp_ds_r - tmp_us_r) / (tmp_us_b - tmp_ds_b);

						// Check if this time is shorter than the previously known shortest time until collision
						if(tmp_t_IS < t_IS_lowest)
						{
							t_IS_lowest = tmp_t_IS;

							// Record which shells these were
							ind_s_ds = tmp_ind_ds_shell;
							ind_s_us = tmp_ind_us_shell;
						}
					}
				}
			} 
		}
	
		/* Calculate time until next FS */
		// If the forward shock has begun, calculate the time until the next FS event		
		if(fs_active == true)
		{
			m_remaining = m_ex - m_swept_tot;
			// Assume a wind density profile for the medium:
			if ((*p_model_params).k_med == 2)
			{
				// Spherical propagation
				t_FS_sweep = m_remaining*m_bar / (4. * M_PI * (*p_model_params).rho_not) / (beta((*p_jet_shells).shell_gamma.at(fs_shell_index))*c_cm); // sec
				// Conical propagation (using jet opening angle)
				// t_FS_sweep = m_remaining*m_bar *2. / (pow(tmp_theta,2.) * (*p_model_params).rho_not) / (beta((*p_jet_shells).shell_gamma.at(fs_shell_index))*c_cm); // sec
			}
			// Assume a constant density profile for the medium:
			else if ((*p_model_params).k_med == 0)
			{
				// Spherical propagation
				t_FS_sweep = ( pow( (3.*m_remaining*m_bar / (4. * M_PI * (*p_model_params).rho_not)) + pow((*p_jet_shells).shell_radius.at(fs_shell_index)*c_cm,3.),1./3.) - ((*p_jet_shells).shell_radius.at(fs_shell_index)*c_cm) )/ (beta((*p_jet_shells).shell_gamma.at(fs_shell_index))*c_cm); // sec
				// Conical propagation (using jet opening angle)
				// t_FS_sweep = ( pow( (3.*m_remaining*m_bar * 2. / ( pow(tmp_theta,2.) * (*p_model_params).rho_not)) + pow((*p_jet_shells).shell_radius.at(fs_shell_index)*c_cm,3.),1./3.) - ((*p_jet_shells).shell_radius.at(fs_shell_index)*c_cm) )/ (beta((*p_jet_shells).shell_gamma.at(fs_shell_index))*c_cm); // sec
				
			}	
		}


		/* Calculate time until next RS */
		// If the reverse shock has begun, calculate the time until the next RS event
		if( (rs_active == true) and ( (*p_jet_shells).shell_gamma.at(active_inds.at(0)) > (*p_jet_shells).shell_gamma.at(rs_shell_index)) )
		{
			
			// Calculate when the outer most shell will collide with the external shock discontinuity, e.g., when a reverse shock will happen
			t_RS_coll = ((*p_jet_shells).shell_radius.at(rs_shell_index) - (*p_jet_shells).shell_radius.at(active_inds.at(0))) / (beta((*p_jet_shells).shell_gamma.at(active_inds.at(0))) - beta((*p_jet_shells).shell_gamma.at(rs_shell_index)));
		}

		/* ## Internal Shock Dynamics and Emission ## */
		if((t_IS_lowest < t_FS_sweep) & (t_IS_lowest < t_RS_coll))
		{
		
			// Keep track of the number of internal shocks that occur, this is used for taking snapshots of the Lorentz distribution
			num_shocks +=1;

			// Move to the emission time of the shell collision (in the rest frame of the jet)
			tmp_te += t_IS_lowest;
			// Move shells forward
			for(int i=0; i<(*p_model_params).numshells; ++i)
			{
				if((*p_jet_shells).shell_status.at(i) == 1)
				{
					(*p_jet_shells).shell_radius.at(i) += beta((*p_jet_shells).shell_gamma.at(i)) * t_IS_lowest;
				}
			}

			// De-active the downstream shell now that it has been merged
			(*p_jet_shells).shell_status.at(ind_s_ds) = 0;

			// Sweep up some mass into the FS
			//  If the FS hasn't begun yet, use the outer most shell to sweep up mass
			if(fs_active == false)
			{
				// Assume a wind density profile for the medium:
				if ((*p_model_params).k_med == 2)
				{
					m_swept = 4. * M_PI * (*p_model_params).rho_not * beta( (*p_jet_shells).shell_gamma.at(active_inds.at(0)) ) * c_cm * t_IS_lowest / m_bar; // g/M
				}
				// Assume a constant density profile for the medium:
				else if ((*p_model_params).k_med == 0)
				{
					tmp_r_before = (*p_jet_shells).shell_radius.at(active_inds.at(0)) - beta((*p_jet_shells).shell_gamma.at(active_inds.at(0))) * t_IS_lowest;
					m_swept = (4./3.) * M_PI * (*p_model_params).rho_not * ( pow((*p_jet_shells).shell_radius.at(active_inds.at(0))*c_cm,3.) - pow(tmp_r_before *c_cm,3.) ) / m_bar; // g/M
				} 
				
				m_swept_tot += m_swept; 

				// If the mass swept up is greater than m_ex, a forward shock begins (and therefor a reverse shock)
				if(m_swept_tot >= m_ex)
				{				

					fs_active = true;
					rs_active = true;

					// Turn the outermost shell into the FS
					fs_shell_index = active_inds.at(0);
					(*p_jet_shells).shell_status.at(fs_shell_index) = 2; // Deactivate the shell so it does not participate in any further IS calculations
					(*p_jet_shells).shell_mass.at(fs_shell_index) += m_swept_tot;

					// Turn the second outermost shell into the RS
					rs_shell_index = active_inds.at(1);
					(*p_jet_shells).shell_status.at(rs_shell_index) = 0; // Deactivate the shell so it does not participate in any further IS calculations

					// Reset swept up mass
					m_swept = 0.;
					m_swept_tot = 0.;
					m_ex = q * ((*p_jet_shells).shell_mass.at(fs_shell_index)+(*p_jet_shells).shell_mass.at(rs_shell_index))/(*p_jet_shells).shell_gamma.at(fs_shell_index);
				}
			}
			// Else, the FS is already active
			else if(fs_active == true)
			{
				// Move the FS forward
				tmp_r_before = (*p_jet_shells).shell_radius.at(fs_shell_index);
				(*p_jet_shells).shell_radius.at(fs_shell_index) += beta((*p_jet_shells).shell_gamma.at(fs_shell_index))*t_IS_lowest;
				(*p_jet_shells).shell_radius.at(rs_shell_index) += beta((*p_jet_shells).shell_gamma.at(rs_shell_index))*t_IS_lowest;

				// Assume a wind density profile for the medium:
				if ((*p_model_params).k_med == 2)
				{
					m_swept = 4. * M_PI * (*p_model_params).rho_not * beta( (*p_jet_shells).shell_gamma.at(fs_shell_index) ) * c_cm * t_IS_lowest/ m_bar; // g/M
				}
				// Assume a constant density profile for the medium:
				else if ((*p_model_params).k_med == 0)
				{
					m_swept = 4. * M_PI * (*p_model_params).rho_not * ( pow((*p_jet_shells).shell_radius.at(fs_shell_index)*c_cm,3.) - pow(tmp_r_before*c_cm ,3.) )/ m_bar; // g/M
				} 
				// Calculate the mass still needed to be swept up until a FS occurs
				m_swept_tot += m_swept;

			}

			// For ease of reading, make variables for the masses and Lorentz factors of the two colliding shells
			// Down stream shell
			shell_ds_m = (*p_jet_shells).shell_mass.at(ind_s_ds);
			shell_ds_g = (*p_jet_shells).shell_gamma.at(ind_s_ds);
			// Upstream shell
			shell_us_m = (*p_jet_shells).shell_mass.at(ind_s_us);
			shell_us_g = (*p_jet_shells).shell_gamma.at(ind_s_us);

			// Record the radius of the collision
			rad_coll = (*p_jet_shells).shell_radius.at(ind_s_us);


			// Calculate optical depth 
			tmp_tau = 0.; // Reset to zero
			if(ind_s_ds > 0) // If this is the the first shell, then there is no material to pass through
			{
				for(int i=0; i<ind_s_ds; ++i)
				{
					if( (*p_jet_shells).shell_status.at(i) == 1)
					{
						tmp_tau += 0.2 * (*p_jet_shells).shell_mass.at(i) * m_bar / (4. * M_PI * pow(c_cm *(*p_jet_shells).shell_radius.at(i) ,2.) ); 
					}
				}
			}

			// Approximate the resulting Lorentz factor from the collision
			// tmp_gamma_r = sqrt(shell_ds_g*shell_us_g);
			// Resulting Lorentz factor of the combined shell (after complete redistribution of momenta)
			tmp_gamma_r = _calc_shell_coll_gamma(shell_ds_g, shell_us_g, shell_ds_m, shell_us_m);

			// Calculate the energy dissipated in this collision, using the smaller mass of the two shells:
			if(shell_us_m > shell_ds_m)
			{
				min_mass = shell_ds_m;
			}
			else
			{
				min_mass = shell_us_m;
			}

			// Energy dissipated in shock
			tmp_e_diss = m_bar * pow(c_cm,2.) * (shell_ds_m*shell_ds_g + shell_us_m*shell_us_g - (shell_us_m + shell_ds_m)*tmp_gamma_r); // erg
			
			// gamma_int = 0.5 *( sqrt(shell_ds_g/shell_us_g) + sqrt(shell_us_g/shell_ds_g) ); // Lorentz factor for internal motion in shocked material
			// tmp_eps_star = (gamma_int - 1) * pow(c_cm,2.); // erg / g, Average proton factor from the collision of two shells
			tmp_eps_star = tmp_e_diss/(2.*m_bar*min_mass*tmp_gamma_r);

			// The typical Lorentz factor of the electron distribution (Gamma_min)
			// Or the minimum Lorentz factor for a population of electrons distributed as a power law with index p
			if(const_gamma_e == false)
			{
				tmp_gamma_e = ((*p_model_params).p_int-2.) * (*p_model_params).eps_e_int * mp * tmp_eps_star / ((*p_model_params).p_int-1.) / (*p_model_params).zeta_int / me / pow(c_cm,2.);
			}
			else
			{
				tmp_gamma_e = 15000; // Constant
			}

			tmp_rho = E_dot_kin / (4. * M_PI * pow(rad_coll*gamma_bar,2.) * pow(c_cm,5.)); // g/cm^3, Comoving proton mass density
			tmp_beq = sqrt( 8.* M_PI * (*p_model_params).eps_b_int * tmp_rho * tmp_eps_star); // (erg/cm^3)^1/2, Magnetic field density, assuming B = Beq  
			tmp_esyn_kev = 50.*(tmp_gamma_r/300.)*(tmp_beq/1e3)*pow((tmp_gamma_e/100.),2.) / 1000.; // keV, Synchrotron energy in the rest frame
			tmp_esyn_erg = kev_to_erg * tmp_esyn_kev; // erg, Synchrotron energy in the rest frame
			
			w = tmp_gamma_e*(tmp_esyn_erg/tmp_gamma_r)/(me * pow(c_cm,2.) ); // Indicates whether we are in the Klein-Nishina regime or Thompson (i.e., w>>1 or w<<1, respectively)
			// w = 33.*(tmp_beq/1e3)*pow((tmp_gamma_e/1e4),3.); // Alternative relation to find w 

			quadr_const = -(3./2.)*(0.2/M_PI/pow(c_cm*rad_coll,2.)) * tmp_gamma_e * 100. * (E_dot_kin/pow(gamma_bar*c_cm,2.)) * pow(tmp_beq/1000.,-2.); // Constant in the quadratic equation
			// Calculate Compton parameter
			Q_IC = (-1. + sqrt(1.-4.*quadr_const))/2.;
			if(w>=1.) // Klein-Nishina regime 
			{
				// Calculate fraction of energy that goes into Inverse Compton electron
				alpha_ic = (Q_IC/w)/(1.+(Q_IC/w));
			}
			else // w < 1, Thompson regime 
			{
				alpha_ic = Q_IC/(1.+Q_IC);
			}

			tmp_asyn = 1. - alpha_ic; // Fraction of energy remaining for synchrotron electrons

			// Duration of the emission event 
			tmp_delt = rad_coll/2./pow(tmp_gamma_r,2.);

			// tmp_lum_diss = tmp_e_diss* (*p_model_params).eps_e_int * tmp_asyn / tmp_delt; // erg / s
			tmp_lum_diss = tmp_e_diss* (*p_model_params).eps_e_int * tmp_asyn; // erg 		

			// Is the emission efficient? 
			t_syn = 6.*pow(tmp_gamma_e/100.,-1.)*pow(tmp_beq/1000.,-2.); // Synchrotron time-scale
			tmp_eff = (t_syn < ((1.+Q_IC)*rad_coll/gamma_bar) );
			// tmp_eff = true;

			// Calculate nu_c
			tmp_nu_c = (18.*M_PI*me*qe*c_cm/pow(sigma_T,2.))/pow(tmp_beq,3.)/pow(tmp_delt,2.);
			// Calculate nu_m
			tmp_nu_m = (1./2./M_PI) * pow( ( (*p_model_params).p_int-2.)*(*p_model_params).eps_e_int/((*p_model_params).p_int-1.)/(*p_model_params).zeta_int,2.) * (qe*pow(mp,2.)/pow(me,3.)/pow(c_cm,5.)) * pow(tmp_eps_star,2.) * tmp_beq;
			

			// Set the Lorentz factor of the merged shells
			(*p_jet_shells).shell_gamma.at(ind_s_us) = tmp_gamma_r;
			// Set the mass of the merged shell (i.e., add the mass of the two shells)
			(*p_jet_shells).shell_mass.at(ind_s_us) += (*p_jet_shells).shell_mass.at(ind_s_ds);


			// If the emission is efficient, add the contribution
			// E.g., if the emission time is less than the shell expansion (the dynamical time scale of the shell)		
			if ( tmp_eff )
			{
				// Calculate the arrival time at the observer and the width of the emission 
				te_is.push_back(tmp_te); // sec, Time of emission (in the rest frame of the jet), is also equal to the time a shell crosses the photosphere
				ta_is.push_back(tmp_te - rad_coll); // sec, Time when the emission arrives at the observer (in the observer frame)
				delt_is.push_back(tmp_delt); // sec, Width of the emission (in observer frame)
				beq_is.push_back(tmp_beq); // (erg/cm^3)^1/2, Magnetic field density, assuming B = Beq 
				gamma_e_is.push_back(tmp_gamma_e); // Lorentz factor of the electron population
				esyn_is.push_back(tmp_esyn_kev); // keV, Synchrotron energy emitted by accelerated electron
				gamma_r_is.push_back(tmp_gamma_r); // Approximation of the combined Lorentz factor of the colliding shells
				e_diss_is.push_back(tmp_lum_diss); // erg, Dissipated energy during the collision 
				nu_c_is.push_back(tmp_nu_c); // Hz, Critical synchrotron frequency
				nu_m_is.push_back(tmp_nu_m); // Hz, Minimum electron frequency

				shell_ind_is.push_back(ind_s_ds); // Downstream shell swept up in the internal shock				
				asyn_is.push_back(tmp_asyn); // Fraction of the energy in electrons which goes into synchrotron (as opposed to Inverse Compton)
				tau.push_back(tmp_tau); // Optical depth at the location of the collision
				relvel.push_back( _rel_vel(shell_us_g,shell_ds_g) ); // Relative velocity between the two shells

				// eps_star_is.push_back(tmp_eps_star); // erg / g,s Internal energy dissipated in a collision 
				// rho_is.push_back(rho); // g cm^-3, Density of the collision region
			}
		}		
		
		/* ## Forward Shock Dynamics and Emission ## */
		else if((fs_active == true ) & (t_FS_sweep < t_RS_coll))
		{				
			// Keep track of the number of forward shocks that occur, this is used for taking snapshots of the Lorentz distribution
			num_fs_shocks += q;
			if(num_fs_shocks > 1)
			{
				num_shocks += 1;
				num_fs_shocks = 0.;
			}

			// Move to the emission time of the shell collision (in the rest frame of the jet)
			tmp_te += t_FS_sweep;

			// Move shells forward
			for(int i=0; i<(*p_model_params).numshells; ++i)
			{
				if((*p_jet_shells).shell_status.at(i) == 1)
				{
					(*p_jet_shells).shell_radius.at(i) += beta((*p_jet_shells).shell_gamma.at(i)) * t_FS_sweep;
				}
			}
			// Move FS and RS forward
			tmp_r_before = (*p_jet_shells).shell_radius.at(fs_shell_index);
			(*p_jet_shells).shell_radius.at(fs_shell_index) += beta((*p_jet_shells).shell_gamma.at(fs_shell_index))*t_FS_sweep;
			(*p_jet_shells).shell_radius.at(rs_shell_index) += beta((*p_jet_shells).shell_gamma.at(rs_shell_index))*t_FS_sweep;

			// Sweep up some mass into the FS
			// Assume a wind density profile for the medium:
			if ((*p_model_params).k_med == 2)
			{
				m_swept = 4. * M_PI * (*p_model_params).rho_not * beta( (*p_jet_shells).shell_gamma.at(fs_shell_index) ) * c_cm * t_FS_sweep/ m_bar; // g/M
			}
			// Assume a constant density profile for the medium:
			else if ((*p_model_params).k_med == 0)
			{
				m_swept = (4./3.) * M_PI * (*p_model_params).rho_not * ( pow((*p_jet_shells).shell_radius.at(fs_shell_index)*c_cm,3.) - pow(tmp_r_before*c_cm ,3.) )/ m_bar; // g/M
			}
			m_swept_tot += m_swept;

			/*
			Forward shock dynamics
			*/ 

			// For ease of reading, make variables for the masses and Lorentz factors of the two colliding shells
			tmp_fs_m = (*p_jet_shells).shell_mass.at(fs_shell_index); // current mass of the forward shock shell
			tmp_fs_g = (*p_jet_shells).shell_gamma.at(fs_shell_index); // current Lorentz factor of the forward shock shell
			tmp_fs_r = (*p_jet_shells).shell_radius.at(fs_shell_index); // current radius of the forward shock shell

			tmp_rs_m = (*p_jet_shells).shell_mass.at(rs_shell_index); // current mass of the reverse shock shell
			tmp_rs_g = (*p_jet_shells).shell_gamma.at(rs_shell_index); // current Lorentz factor of the reverse shock shell
			tmp_rs_r = (*p_jet_shells).shell_radius.at(rs_shell_index); // current radius of the reverse shock shell
			
			// Calculate necessary parameters:
			// The new Lorentz factor at the discontinuity due to forward shock
			tmp_gamma_r = pow( ((tmp_rs_m + tmp_fs_m*fs_gamma_int)*pow(tmp_fs_g,2.) + m_swept_tot*tmp_fs_g ) / ( (tmp_rs_m + tmp_fs_m*fs_gamma_int) + 2.*m_swept_tot*tmp_fs_g ) ,0.5);
			if(tmp_gamma_r < 1)
			{
				tmp_gamma_r = 1.000001;
			}

			// Opening angle of the jet increases due to lateral expansion
			if(flag_lat_exp == true)
			{
				tmp_theta += t_FS_sweep / tmp_gamma_r / tmp_fs_r / sqrt(3.);
			}

			// The energy dissipated by forward shock
			tmp_e_diss = pow(c_cm,2.) * m_bar * ( (tmp_rs_m + tmp_fs_m*fs_gamma_int)*tmp_fs_g + m_swept_tot - (m_swept_tot + tmp_rs_m + tmp_fs_m*fs_gamma_int )*tmp_gamma_r); // erg 
			// tmp_e_diss = pow(c_cm,2.) * m_bar * m_swept_tot * ( tmp_fs_g + 1 - 2*tmp_gamma_r); // erg 

			// tmp_eps_star =  (fs_gamma_int - 1) * pow(c_cm,2.); // erg/g, Average energy dissipated for each proton in the collision of two shells
			tmp_eps_star = tmp_e_diss /(2.*m_bar*m_swept_tot*tmp_gamma_r); // erg/g, Average energy dissipated for each proton in the collision of two shells

			// The typical Lorentz factor of the electron distribution (Gamma_min)
			// Or the minimum Lorentz factor for a population of electrons distributed as a power law with index p
			if(const_gamma_e == false)
			{
				tmp_gamma_e = tmp_gamma_r * ((*p_model_params).p_ext-2.) * (*p_model_params).eps_e_ext * mp / ((*p_model_params).p_ext-1.) /(*p_model_params).zeta_ext / me; 

				// Is this an equivalent expression? Why or why not?
				// tmp_gamma_e = ((*p_model_params).p_ext-2.) * (*p_model_params).eps_e_ext * mp * tmp_eps_star / ((*p_model_params).p_ext-1.) /(*p_model_params).zeta_ext / me / pow(c_cm,2.);
			}
			else
			{
				tmp_gamma_e = 15000; // Constant
			}
			
			// The magnetic field in the forward shock
			tmp_rho = tmp_gamma_r*(*p_model_params).rho_not/pow(tmp_fs_r*c_cm,(*p_model_params).k_med); // g/cm^3, comoving proton density 
			tmp_beq = c_cm*tmp_gamma_r*sqrt(32. * M_PI * (*p_model_params).eps_b_ext * tmp_rho); // (erg/cm^3)^1/2, Magnetic field density, assuming B = Beq  

			// The synchrotron energy of the accelerated electrons 
			tmp_esyn_kev = 50.*(tmp_gamma_r/300.)*(tmp_beq/1e3)*pow((tmp_gamma_e/100.),2.) / 1000.; // keV, Synchrotron energy in the rest frame
			tmp_esyn_erg = kev_to_erg * tmp_esyn_kev; // erg, Synchrotron energy in the rest frame

			w = tmp_gamma_e*(tmp_esyn_erg/tmp_gamma_r)/(me * pow(c_cm,2.) ); // Indicates whether we are in the Klein-Nishina regime or Thompson (i.e., w>>1 or w<<1, respectively)
			quadr_const = -(3./2.)*(0.2/M_PI/pow(c_cm*tmp_fs_r,2.)) * tmp_gamma_e * 100. * (E_dot_kin/pow(gamma_bar*c_cm,2.)) * pow(tmp_beq/1000.,-2.); // Constant in the quadratic equation
			// Calculate Compton parameter
			Q_IC = (-1. + sqrt(1.-4.*quadr_const))/2.;
			if(w>=1.) // Klein-Nishina regime 
			{
				// Calculate fraction of energy that goes into Inverse Compton electron
				alpha_ic = (Q_IC/w)/(1.+(Q_IC/w));
			}
			else // w < 1, Thompson regime 
			{
				alpha_ic = Q_IC/(1.+Q_IC);
			}

			tmp_asyn = 1 - alpha_ic; // Fraction of energy remaining for synchrotron electrons


			// If the jet is still before the the jet break
			if( ( tmp_theta * tmp_gamma_r) >= 1. )
			{
				// Duration of the emission event 
				tmp_delt = tmp_fs_r/2./pow(tmp_gamma_r,2.);
			}
			// If the jet is after the jet break
			else
			{
				// Duration of the emission event 
				tmp_delt = tmp_fs_r * pow(tmp_theta,2.) / 2.;
			}

			// Luminosity of the emission
			// tmp_lum_diss = tmp_e_diss* (*p_model_params).eps_e_ext * tmp_asyn / tmp_delt; // erg / s;  
			tmp_lum_diss = tmp_e_diss* (*p_model_params).eps_e_ext * tmp_asyn; // erg 

			// Calculate nu_c
			tmp_nu_c = (18.*M_PI*me*qe*c_cm/pow(sigma_T,2.))/pow(tmp_beq,3.)/pow(tmp_delt,2.);
			// Calculate nu_m
			tmp_nu_m = (1./2./M_PI) * pow(((*p_model_params).p_ext-2.)*(*p_model_params).eps_e_ext/((*p_model_params).p_ext-1.)/(*p_model_params).zeta_ext,2.) * (qe*pow(mp,2.)/pow(me,3.)/pow(c_cm,5.)) * pow(tmp_eps_star,2.) * tmp_beq;

			// Update the mass of the forward shock shell
			(*p_jet_shells).shell_mass.at(fs_shell_index) += m_swept_tot;

			// Update the Lorentz factor of the forward and reverse shock shells
			(*p_jet_shells).shell_gamma.at(fs_shell_index) = tmp_gamma_r;
			(*p_jet_shells).shell_gamma.at(rs_shell_index) = tmp_gamma_r;

			// The new internal Lorentz factor in the shocked forward medium due to forward shock, if:
			// adiabatic (slow cooling) so FS_gamma_int > 1 and approx = gamma_r  :
			fs_gamma_int = ((tmp_rs_m + tmp_fs_m*fs_gamma_int)*tmp_fs_g + m_swept_tot - tmp_rs_m*tmp_gamma_r)/(tmp_fs_m+m_swept_tot)/tmp_gamma_r;

			// Is the emission efficient? 
			// t_syn = 6.*pow(tmp_gamma_e/100.,-1.)*pow(tmp_beq/1000.,-2.); // Synchrotron time-scale
			// tmp_eff = (t_syn < ((1.+Q_IC)*tmp_fs_r/gamma_bar) );
			tmp_eff = true;
			
			// If the emission is efficient, add the contribution
			// E.g., if the emission time is less than the shell expansion (the dynamical time scale of the shell)		
			if ( tmp_eff )
			{
				// Continuous emission:
				te_fs.push_back(tmp_te); // sec, Time of emission 
				ta_fs.push_back(tmp_te - tmp_fs_r); // sec, Arrival time of emission at observer
				delt_fs.push_back(tmp_delt); // sec, Duration of the emission
				beq_fs.push_back(tmp_beq); // (erg/cm^3)^1/2, Magnetic field density, assuming B = Beq 
				gamma_e_fs.push_back(tmp_gamma_e); // Lorentz factor of the electron population
				esyn_fs.push_back(tmp_esyn_kev); // erg, The typical synchrotron energy of the accelerated electrons
				gamma_r_fs.push_back(tmp_gamma_r); // New Lorentz factor of the FS
				e_diss_fs.push_back(tmp_lum_diss); // erg/s, Dissipated energy of the FS	
				nu_c_fs.push_back(tmp_nu_c); // Hz, Critical synchrotron frequency
				nu_m_fs.push_back(tmp_nu_m); // Hz, Minimum electron frequency

				rad_coll_fs.push_back(tmp_fs_r); // light seconds, Radius of forward shell
				rho_fs.push_back(tmp_rho); // g cm^-3, Density of the collision region
				eps_star_fs.push_back(tmp_eps_star); // Internal energy dissipated in a collision 
				num_swept_e_fs.push_back((*p_model_params).zeta_ext * m_swept_tot / mp); // Number of swept up electrons in the shock
				theta_fs.push_back(tmp_theta); // rad, Jet opening angle, assumes sound speed = c / sqrt(3)

				mass_te.push_back(tmp_te); // s; Source frame time when mass is added to external shock material
				mass_ta.push_back(tmp_te - tmp_fs_r); // s; Observer frame time when mass is added to external shock material
				mass_rs.push_back(tmp_rs_m*m_bar); // g; Mass of material previously passed by reverse shock
				mass_fs.push_back(tmp_fs_m*m_bar); // g; Mass of material previously passed by forward shock
				mass_ej.push_back(0.); // g; Mass of material colliding with reverse shock material
			}


			// Reset the swept up mass
			m_swept = 0.;
			m_swept_tot = 0.;
			m_ex = q * ((*p_jet_shells).shell_mass.at(fs_shell_index)+(*p_jet_shells).shell_mass.at(rs_shell_index))/gamma_bar;

		}

		/* ## Reverse Shock Dynamics and Emission ## */
		else if((rs_active == true ))
		{			
			// Keep track of the number of reverse shocks that occur, this is used for taking snapshots of the Lorentz distribution
			num_shocks +=1;

			// Move to the emission time of the shell collision (in the rest frame of the jet)
			tmp_te += t_RS_coll;

			// Move shells forward
			for(int i=0; i<(*p_model_params).numshells; ++i)
			{
				if((*p_jet_shells).shell_status.at(i) == 1)
				{
					(*p_jet_shells).shell_radius.at(i) += beta((*p_jet_shells).shell_gamma.at(i)) * t_RS_coll;
				}
			}
			// Move FS and RS forward
			tmp_r_before = (*p_jet_shells).shell_radius.at(fs_shell_index);
			(*p_jet_shells).shell_radius.at(fs_shell_index) += beta((*p_jet_shells).shell_gamma.at(fs_shell_index))*t_RS_coll;
			(*p_jet_shells).shell_radius.at(rs_shell_index) += beta((*p_jet_shells).shell_gamma.at(rs_shell_index))*t_RS_coll;

			// Sweep up some mass into the FS
			// Assume a wind density profile for the medium:
			if ((*p_model_params).k_med == 2)
			{
				m_swept = 4. * M_PI * (*p_model_params).rho_not * beta( (*p_jet_shells).shell_gamma.at(fs_shell_index) ) * c_cm * t_RS_coll/m_bar; // g/M
			}
			// Assume a constant density profile for the medium:
			else if ((*p_model_params).k_med == 0)
			{
				m_swept = (4./3.)* M_PI * (*p_model_params).rho_not * ( pow((*p_jet_shells).shell_radius.at(fs_shell_index)*c_cm,3.) - pow(tmp_r_before*c_cm ,3.) )/m_bar; // g/M
			} 
			// Calculate the mass still needed to be swept up until a FS occurs
			m_swept_tot += m_swept; 


			/*
			Reverse shock dynamics 
			*/ 
			
			// For ease of reading, make variables for the masses and Lorentz factors of the two colliding shells
			// Select the ejecta shell that is colliding with the RS 
			tmp_ej_m = (*p_jet_shells).shell_mass.at(active_inds.at(0)); // The mass of the shell that is colliding with the reverse shock 
			tmp_ej_g = (*p_jet_shells).shell_gamma.at(active_inds.at(0)); // The lorentz factor of the shell that is colliding with the reverse shock 
			
			tmp_fs_m = (*p_jet_shells).shell_mass.at(fs_shell_index); // current mass of the forward shock shell
			tmp_fs_g = (*p_jet_shells).shell_gamma.at(fs_shell_index); // current Lorentz factor of the forward shock shell
			tmp_fs_r = (*p_jet_shells).shell_radius.at(fs_shell_index); // current radius of the forward shock shell

			tmp_rs_m = (*p_jet_shells).shell_mass.at(rs_shell_index); // current mass of the reverse shock shell
			tmp_rs_g = (*p_jet_shells).shell_gamma.at(rs_shell_index); // current Lorentz factor of the reverse shock shell
			tmp_rs_r = (*p_jet_shells).shell_radius.at(rs_shell_index); // current radius of the reverse shock shell

			
			// Calculate necessary parameters: 
			// The new Lorentz factor at the discontinuity due to reverse shock 			
			tmp_gamma_r = pow(tmp_rs_g*tmp_ej_g,0.5)*pow( ((tmp_rs_m + tmp_fs_m*fs_gamma_int)*tmp_rs_g + tmp_ej_m*tmp_ej_g )/((tmp_rs_m + tmp_fs_m*fs_gamma_int )*tmp_ej_g + tmp_ej_m*tmp_rs_g ),0.5);
			
			// Calculate the energy dissipated in this collision, using the smaller mass of the two shells:
			if(tmp_rs_m > tmp_ej_m)
			{
				min_mass = tmp_ej_m;
			}
			else
			{
				min_mass = tmp_rs_m;
			}

			// The energy dissipated by reverse shock
			tmp_e_diss = pow(c_cm,2.) * m_bar * ( (tmp_rs_m + tmp_fs_m*fs_gamma_int)*tmp_rs_g + tmp_ej_m*tmp_ej_g - (tmp_ej_m + tmp_rs_m + tmp_fs_m*fs_gamma_int)*tmp_gamma_r); // erg

			// gamma_int = 0.5 *( sqrt(tmp_gamma_r/tmp_ej_g) + sqrt(tmp_ej_g/tmp_gamma_r) ); // Lorentz factor for internal motion in shocked material
			// tmp_eps_star =  (gamma_int - 1.) * pow(c_cm,2.); // erg / g, Average energy per proton within the collided shells

			tmp_eps_star = tmp_e_diss /(2.*m_bar*min_mass*tmp_gamma_r);

			// The typical Lorentz factor of the electron distribution (Gamma_min)
			// Or the minimum Lorentz factor for a population of electrons distributed as a power law with index p
			if(const_gamma_e == false)
			{
				tmp_gamma_e = ((*p_model_params).p_int-2.) * (*p_model_params).eps_e_int * mp * tmp_eps_star / ((*p_model_params).p_int-1.) / (*p_model_params).zeta_int / me / pow(c_cm,2.);
			}
			else
			{
				tmp_gamma_e = 15000; // Constant
			}


			// Density in the reverse shock medium			
			tmp_rho = E_dot_kin / (4. * M_PI * pow(tmp_rs_r*gamma_bar,2.) * pow(c_cm,5.)); // g/cm^3, Comoving proton mass density
			
			// The magnetic field in the reverse shock
			tmp_beq = sqrt(8. * M_PI * (*p_model_params).eps_b_int * tmp_rho * tmp_eps_star); // (erg/cm^3)^1/2, Magnetic field density, assuming B = Beq  

			// The typical energy of an accelerated electron
			tmp_esyn_kev = 50.*(tmp_gamma_r/300.)*(tmp_beq/1e3)*pow((tmp_gamma_e/100.),2.) / 1000.; // keV, Synchrotron energy in the source frame
			tmp_esyn_erg = kev_to_erg * tmp_esyn_kev; // erg, Synchrotron energy in the rest frame 
			
			w = tmp_gamma_e*(tmp_esyn_erg/tmp_gamma_r)/(me * pow(c_cm,2.) ); // Indicates whether we are in the Klein-Nishina regime or Thompson (i.e., w>>1 or w<<1, respectively)
			quadr_const = -(3./2.)*(0.2/M_PI/pow(c_cm*tmp_rs_r,2.)) * tmp_gamma_e * 100. * (E_dot_kin/pow(gamma_bar*c_cm,2.)) * pow(tmp_beq/1000.,-2.); // Constant in the quadratic equation
			// Calculate Compton parameter
			Q_IC = (-1. + sqrt(1.-4.*quadr_const))/2.;
			if(w>=1.) // Klein-Nishina regime 
			{
				// Calculate fraction of energy that goes into Inverse Compton electron
				alpha_ic = (Q_IC/w)/(1.+(Q_IC/w));
			}
			else // w < 1, Thompson regime 
			{
				alpha_ic = Q_IC/(1.+Q_IC);
			}

			tmp_asyn = 1. - alpha_ic; // Fraction of energy remaining for synchrotron electrons

			// Duration of the emission event 
			// If the jet is still before the the jet break
			if( ( (*p_model_params).theta * tmp_gamma_r) >= 1. )
			{
				// Duration of the emission event 
				tmp_delt = tmp_rs_r/2./pow(tmp_gamma_r,2.);
			}
			// If the jet is after the jet break
			else
			{				
				// Duration of the emission event 
				tmp_delt = tmp_rs_r * pow((*p_model_params).theta,2.) / 2.;
			}

			// Luminosity of the emission
			// tmp_lum_diss = tmp_e_diss * (*p_model_params).eps_e_int * tmp_asyn / tmp_delt; // erg / s;
			tmp_lum_diss = tmp_e_diss * (*p_model_params).eps_e_int * tmp_asyn; // erg ;
			
			// Calculate nu_c
			tmp_nu_c = (18.*M_PI*me*qe*c_cm/pow(sigma_T,2.))/pow(tmp_beq,3.)/pow(tmp_delt,2.);
			// Calculate nu_m
			tmp_nu_m = (1./2./M_PI) * pow(((*p_model_params).p_int-2.)*(*p_model_params).eps_e_int/((*p_model_params).p_int-1.)/(*p_model_params).zeta_int,2.) * (qe*pow(mp,2.)/pow(me,3.)/pow(c_cm,5.)) * pow(tmp_eps_star,2.) * tmp_beq;

			t_syn = 6.*pow(tmp_gamma_e/100.,-1.)*pow(tmp_beq/1000.,-2.); // Synchrotron time-scale
			// tmp_eff 	= (t_syn < ((1.+Q_IC)*tmp_rs_r/gamma_bar) );
			tmp_eff = true;


			// Update the mass of the reverse shock shell
			(*p_jet_shells).shell_mass.at(rs_shell_index) += tmp_ej_m;
			// Update the Lorentz factor of the forward and reverse shock shells
			(*p_jet_shells).shell_gamma.at(fs_shell_index) = tmp_gamma_r;
			(*p_jet_shells).shell_gamma.at(rs_shell_index) = tmp_gamma_r;

			// Deactivate internal shell that collided with RS
			(*p_jet_shells).shell_status.at(active_inds.at(0)) = 0;

			// If the emission is efficient, add the contribution
			// E.g., if the emission time is less than the shell expansion (the dynamical time scale of the shell)		
			if ( tmp_eff )
			{
				// Record necessary parameters to calculate the emission
				te_rs.push_back(tmp_te); // sec, Time of emission 
				ta_rs.push_back(tmp_te - tmp_rs_r); // sec, Arrival time of emission at observer
				delt_rs.push_back(tmp_delt); // sec, Duration of the emission
				beq_rs.push_back(tmp_beq); // (erg/cm^3)^1/2, Magnetic field density, assuming B = Beq 
				gamma_e_rs.push_back(tmp_gamma_e); // Lorentz factor of the electron population	
				esyn_rs.push_back(tmp_esyn_kev); // erg, The typical synchrotron energy of the accelerated electrons
				gamma_r_rs.push_back(tmp_gamma_r); // New Lorentz factor of the RS
				e_diss_rs.push_back(tmp_lum_diss); // erg, Dissipated energy of the RS
				nu_c_rs.push_back(tmp_nu_c); // Hz, Critical synchrotron frequency
				nu_m_rs.push_back(tmp_nu_m); // Hz, Minimum electron frequency
				shell_ind_rs.push_back(active_inds.at(0)); // The index of the shell which was crossed by the reverse shock
				// eps_star_rs.push_back(tmp_eps_star); // erg / g, Internal energy dissipated in a collision 
				// rho_rs.push_back(rho); // g cm^-3, Density of the collision region
				
				mass_te.push_back(tmp_te); // s; Source frame time when mass is added to external shock material
				mass_ta.push_back(tmp_te - tmp_fs_r); // s; Observer frame time when mass is added to external shock material
				mass_rs.push_back(tmp_rs_m*m_bar); // g; Mass of material previously passed by reverse shock
				mass_fs.push_back(tmp_fs_m*m_bar); // g; Mass of material previously passed by forward shock
				mass_ej.push_back(tmp_ej_m*m_bar); // g; Mass of material colliding with reverse shock material
			}

		}
		else
		{
			fs_active = true;
			rs_active = true;

			// Turn the outermost shell into the FS
			fs_shell_index = active_inds.at(0);
			(*p_jet_shells).shell_status.at(fs_shell_index) = 2; // Deactivate the shell so it does not participate in any further IS calculations
			(*p_jet_shells).shell_mass.at(fs_shell_index) += m_swept_tot;

			// Turn the second outermost shell into the RS
			rs_shell_index = active_inds.at(1);
			(*p_jet_shells).shell_status.at(rs_shell_index) = 0; // Deactivate the shell so it does not participate in any further IS calculations

			// Reset swept up mass
			m_swept = 0.;
			m_swept_tot = 0.;
			m_ex = q * ((*p_jet_shells).shell_mass.at(fs_shell_index)+(*p_jet_shells).shell_mass.at(rs_shell_index))/(*p_jet_shells).shell_gamma.at(fs_shell_index);
		}


		// Check if any more collisions will occur
		if(_check_if_sorted())
		{
			ord_lorentz = true;
			// std::cout << "At time t=" << tmp_te << " s, all shells have been launched and Lorentz factors are ordered.\n";
		}		
		// Check if the FS is no longer relativistic 
		if((*p_jet_shells).shell_gamma.at(fs_shell_index) < 2)
		{
			fs_relativ = false;
			// std::cout << "At time t=" << tmp_te << " s, the forward shock is no longer relativistic.\n";
		}
		// Check if all internal shells have been passed through by the reverse shock
		if(active_inds.size()<=1)
		{
			rs_active = false;
			// std::cout << "At time t=" << tmp_te << " s, all shells have passed through the reverse shock.\n";
		}
		if( (anim_lor_dist == true ) & (num_shocks % 100 == 0))
		{
			// Append the current Lorentz distribution to the data file
			(*p_jet_shells).WriteToTXT("data-file-dir/synthGRB_shell_dist.txt", tmp_te, true);
			num_shocks +=1;
		}
	}	
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make the source spectrum using the emission data
void SynthGRB::make_source_spectrum(float energ_min, float energ_max, int num_energ_bins, float tmin, float tmax, std::string comp)
{
	// Release previous spectrum if it exists
	if(p_source_spectrum != NULL)
	{
		delete p_source_spectrum;
	}
	// Make a Spectrum object
	p_source_spectrum = new Spectrum(energ_min, energ_max, num_energ_bins);
	(*p_source_spectrum).ZeroSpectrum(); // Set spectrum to zero.


	if(comp.compare("all") == 0)
	{
		MakeThermalSpec(p_source_spectrum, tmin, tmax);
		MakeISSpec(p_source_spectrum, tmin, tmax);		
		MakeExtShockSpec(p_source_spectrum, tmin, tmax);
	}
	else if(comp.compare("TH") == 0)
	{
		MakeThermalSpec(p_source_spectrum, tmin, tmax);
	}
	else if(comp.compare("IS") == 0)
	{
		MakeISSpec(p_source_spectrum, tmin, tmax);
	}
	else if(comp.compare("FS") == 0)
	{
		MakeFSSpec(p_source_spectrum, tmin, tmax);
	}
	else if(comp.compare("RS") == 0)
	{
		MakeRSSpec(p_source_spectrum, tmin, tmax);
	}
	else
	{
		cout << "Please enter a valid component to plot a light curve for." << endl;
		return;
	}

	// Normalize by time bin
	for(size_t i=0; i< (*p_source_spectrum).spectrum_rate.size();++i)
	{
		(*p_source_spectrum).spectrum_rate.at(i) /= (tmax - tmin);
	}
	(*p_source_spectrum).spectrum_sum /= (tmax - tmin);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Thermal spectrum member functions*/

// Calls function to calculate the thermal spectrum rate for each energy bin
void SynthGRB::MakeThermalSpec(Spectrum * therm_spectrum, float tmin, float tmax)
{
	float num = 0;
	// For each thermal emission event, calculate the emitted spectrum
	for(size_t i=0; i < te_therm.size(); ++i)
	{
		// We only want to take the emission that occurs between the specified tmin and tmax. 
		// The emission occurs between (ta+delt), if any of it overlaps with tmin and tmax, calculate its contribution.
		if ( ta_therm.at(i) <= tmax and (ta_therm.at(i)+delt_therm.at(i)) >= tmin )
		{
			// Integrated emission profile factor
			// double profile_factor = _calc_pulse_profile_factor(ta_therm.at(i), delt_therm.at(i), tmin, tmax, (*p_model_params).theta, (te_therm.at(i) - ta_therm.at(i)));


			// Call function to calculate thermal spectrum rate
			CalcThermalContribution(therm_spectrum, T_phot.at(i),L_phot.at(i));
			num+=1;
		}
		// Else, if this emission occurred outside of the time interval, don't add it to the spectrum.
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Calculate the thermal spectrum from the given temperature and flux of the emission
void SynthGRB::CalcThermalContribution(Spectrum * therm_spectrum, float temp, double flux)
{
	// Calculate the normalization
	double norm = 0.; // Set normalization to zero
	double e_hi = 5e3;
	double e_lo = 1e-6;
	int norm_num_bin= 20.*log10(e_hi/e_lo); // Number of energy bins to use to calculate normalization
	std::vector<float> norm_energy_axis(norm_num_bin+1); // Make vector to store energy axis (initialized to zero's) 
	_make_en_axis(norm_energy_axis,e_lo,e_hi,norm_num_bin+1); // Make energy axis

	float en_curr=0.; // Current energy to evaluate the addition to the normalization
	// For each energy bin along the normalization energy axis, calculate the addition to the normalization and add it. 
	for(int i=0;i < norm_num_bin; ++i)
	{
		en_curr = (norm_energy_axis.at(i+1) + norm_energy_axis.at(i) ) / 2.; // Set the current energy
		if(en_curr < 5.*pow(10.,3.))
		{
			// Calculate the contribution according to a Left-Riemann-Sum
			norm += (norm_energy_axis.at(i+1)-norm_energy_axis.at(i)) * ThermalSpec(en_curr,temp);
		}
		// else: don't contribute to the normalization

	}

	// Calculate the rate of the thermal spectrum
	double tmp_val=0.; // Temporary value to store spectrum contributions
	// For each energy bin along the energy axis
	for(int i=0; i < (*therm_spectrum).num_energ_bins; ++i)
	{
		en_curr = (*therm_spectrum).energ_mid.at(i); // Set the current energy
		
		// The current temp and energy bin define the count rate, the normalization found above is applied.
		// This must still be multiplied by the flux of the source.
		// If the current energy bin is > 2 MeV, then thermal radiation can be ignored.
		if(en_curr < 5.*pow(10.,3.))
		{
			tmp_val = flux * ThermalSpec(en_curr,temp) / norm;
		}
		else{tmp_val=0;}

		// Add the contribution to the total spectrum according to a Left-Riemann-Sum
		(*therm_spectrum).spectrum_sum += ((*therm_spectrum).energ_hi.at(i) - (*therm_spectrum).energ_lo.at(i)) * tmp_val;
		// Check if the spectrum rate per energy bin is requested and store it if so. 
		(*therm_spectrum).spectrum_rate.at(i) += tmp_val;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Thermal spectrum function form based on a modified Planck function. 
double SynthGRB::ThermalSpec(float energy, float temp, float alpha)
{
	double param_list[3] = {temp, alpha, 1.};
	return BB(energy,param_list);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Internal Shock spectrum member functions*/

// Calls function to calculate the spectrum rate of the internal shocks for each energy bin
void SynthGRB::MakeISSpec(Spectrum * intsh_spectrum, float tmin, float tmax)
{
	// For each emission event, calculate the emitted spectrum
	for(size_t i=0; i < te_is.size(); ++i)
	{
		// We only want to take the emission that occurs between the specified Tmin and Tmax. 
		// The emission occurs between (ta+delt), if any of it overlaps with Tmin and Tmax, calculate its contribution.
		if ( (ta_is.at(i) <= tmax) and ( (ta_is.at(i)+ 7.*delt_is.at(i)) >= tmin) )
		{
			// The emission will only be observable if the relativistic velocity is great than the local sound speed v_s/c = 0.1
			// And if the wind is transparent to the radiation
			if ( relvel.at(i) > 0.1 and tau.at(i) < 1)
			{
				// Integrated emission profile factor
				double profile_factor = _calc_pulse_profile_factor(ta_is.at(i), delt_is.at(i), tmin, tmax, (*p_model_params).theta, (te_is.at(i) - ta_is.at(i)));

				// Call function to calculate spectrum count rate, assuming synchrotron emission
				CalcSynchContribution(intsh_spectrum, esyn_is.at(i),e_diss_is.at(i),delt_is.at(i), nu_c_is.at(i), nu_m_is.at(i), (*p_model_params).p_int, beq_is.at(i), gamma_r_is.at(i), profile_factor);
			}	
		}
		// Else, if this emission occurred outside of the time interval, don't add it to the spectrum. 
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/* External shock spectrum member functions*/

// Calls function to calculate the external shock spectrum rate for each energy bin
void SynthGRB::MakeExtShockSpec(Spectrum * extsh_spectrum, float tmin, float tmax)
{
	// Calls function to calculate forward shock spectrum
	MakeFSSpec(extsh_spectrum, tmin, tmax);
	// Calls function to calculate reverse shock spectrum
	MakeRSSpec(extsh_spectrum, tmin, tmax);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Calls function to calculate the external shock spectrum rate for each energy bin
void SynthGRB::MakeFSSpec(Spectrum * extsh_spectrum, float tmin, float tmax)
{
	// For each emission event, calculate the emitted spectrum
	for(size_t i=0; i < te_fs.size(); ++i)
	{
		// We only want to take the emission that occurs between the specified Tmin and Tmax. 
		// The emission occurs between (ta+delt), if any of it overlaps with Tmin and Tmax, calculate its contribution.
		if ( (ta_fs.at(i) <= tmax ) and ( (ta_fs.at(i)+ 7.*delt_fs.at(i)) >= tmin) )
		{
			// Integrated emission profile factor
			double profile_factor = _calc_pulse_profile_factor(ta_fs.at(i), delt_fs.at(i), tmin, tmax, theta_fs.at(i), rad_coll_fs.at(i) );

			// Call function to calculate spectrum count rate, assuming synchrotron emission
			CalcSynchContribution(extsh_spectrum, esyn_fs.at(i),e_diss_fs.at(i),delt_fs.at(i), nu_c_fs.at(i), nu_m_fs.at(i),(*p_model_params).p_ext, beq_fs.at(i), gamma_r_fs.at(i), profile_factor);
		}
		// Else, if this emission occurred outside of the time interval, don't add it to the spectrum. 
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Calls function to calculate the external shock spectrum rate for each energy bin
void SynthGRB::MakeRSSpec(Spectrum * extsh_spectrum, float tmin, float tmax)
{
	// For each emission event, calculate the emitted spectrum
	for(size_t i=0; i < te_rs.size(); ++i)
	{
		// We only want to take the emission that occurs between the specified Tmin and Tmax. 
		// The emission occurs between (ta+delt), if any of it overlaps with Tmin and Tmax, calculate its contribution.
		if ( (ta_rs.at(i) <= tmax) and ( (ta_rs.at(i)+ 7.*delt_rs.at(i)) >= tmin) )
		{
			// Integrated emission profile factor
			double profile_factor = _calc_pulse_profile_factor(ta_rs.at(i), delt_rs.at(i), tmin, tmax, (*p_model_params).theta, (te_rs.at(i) - ta_rs.at(i)));

			// Call function to calculate spectrum count rate, assuming synchrotron emission
			CalcSynchContribution(extsh_spectrum, esyn_rs.at(i),e_diss_rs.at(i),delt_rs.at(i), nu_c_rs.at(i), nu_m_rs.at(i),(*p_model_params).p_int, beq_rs.at(i), gamma_r_rs.at(i), profile_factor);
		}
		// Else, if this emission occurred outside of the time interval, don't add it to the spectrum. 
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Synchrotron spectrum member functions*/

// Calculate the synchrotron spectrum from the synchrotron energy and flux of the emission
void SynthGRB::CalcSynchContribution(Spectrum * synch_spectrum, double esyn, double e_diss, double delt, float nu_c, float nu_m, float p, double B, float Gamma, double profile_factor)
{
	// Calculate the normalization
	double norm = 0.; // Set normalization to zero
	double e_hi = 1e8;
	double e_lo = 1e-6;
	int norm_num_bin= 20.*log10(e_hi/e_lo); // Number of energy bins to use to calculate normalization
	std::vector<float> norm_energy_axis(norm_num_bin+1); // Make vector to store energy axis (initialized to zero's) 
	_make_en_axis(norm_energy_axis,e_lo,e_hi,norm_num_bin+1); // Make energy axis

	float en_curr=0.; // Current energy to evaluate the addition to the normalization
	// For each energy bin along the normalization energy axis, calculate the addition to the normalization and add it. 
	for(int i=0;i < norm_num_bin; ++i)
	{
		en_curr = (norm_energy_axis.at(i+1) + norm_energy_axis.at(i) ) / 2.; // Set the current energy
		// Calculate the contribution according to a Left-Riemann-Sum
		norm += (norm_energy_axis.at(i+1)-norm_energy_axis.at(i)) * SynchSpec(en_curr, esyn, nu_c, nu_m, p, B, Gamma);
	}

	// Calculate the rate of the thermal spectrum
	double tmp_val=0; // Temporary value to store spectrum contributions
	// For each energy bin along the energy axis
	for(int i=0; i < (*synch_spectrum).num_energ_bins; ++i)
	{
		en_curr = (*synch_spectrum).energ_mid.at(i); // Set the current energy

		// The current temp and energy bin define the count rate, the normalization found above is applied.
		// This must still be multiplied by the energy dissipated during the emission event.
		// The energy dissipated can be turned into Flux by dividing the energy dissipated by the observed emission duration (delt).
		tmp_val = e_diss * profile_factor * SynchSpec(en_curr, esyn, nu_c, nu_m, p, B, Gamma) / norm / delt;
		
		// Add the contribution to the total spectrum according to a Center-Riemann-Sum
		(*synch_spectrum).spectrum_sum += ((*synch_spectrum).energ_hi.at(i) - (*synch_spectrum).energ_lo.at(i)) * tmp_val;

		// Check if the spectrum rate per energy bin is requested and store it if so. 
		(*synch_spectrum).spectrum_rate.at(i) += tmp_val;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Synchrotron spectrum function form,
double SynthGRB::SynchSpec(float energy, double esyn, float nu_c, float nu_m, float p, double B, float Gamma)
{
	// Synchrotron 
	double param_list[6] = {nu_c*Gamma, nu_m*Gamma, p, B*Gamma, Gamma, 1.};
	return Synchrotron(energy,  param_list);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Broken Power Law spectrum function form,
// Defaults alpha = -1.5, beta = -2.5
double SynthGRB::BPLSpec(float energy, double esyn, float alpha, float beta)
{
	// BPL 
	double param_list[4] = {esyn, alpha, beta, 1.};
	return (1./esyn) * BPL(energy,  param_list);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make the source light curve using the emission data and specify the time step of the light curve
void SynthGRB::make_source_light_curve(float energ_min, float energ_max, float Tstart, float Tend, float dt, std::string comp, bool logscale)
{
	// Release previous spectrum if it exists
	if(p_source_light_curve != NULL)
	{
		delete p_source_light_curve;
	} 
	// Make a LightCurve object
	p_source_light_curve = new LightCurve(Tstart , Tend, dt, logscale);
	(*p_source_light_curve).ZeroLightCurve(); // Reset spectrum

	// Define the number of energy bins
	float tmp_num_energ_bins = log10(energ_max/energ_min)*20.;

	// For each time bin, calculate the photon rate.
	for(int i=0; i < (*p_source_light_curve).num_time_bins-1; ++i)
	{
		// Find the spectrum sum for each emission event which occurs between (light_curve_time[i], light_curve_time[i+1]) 
		make_source_spectrum(energ_min, energ_max, tmp_num_energ_bins, (*p_source_light_curve).lc_time.at(i), (*p_source_light_curve).lc_time.at(i+1), comp);		
		(*p_source_light_curve).lc_rate.at(i) = (*p_source_spectrum).spectrum_sum; // Because we use ediss/delt to calculate the spectra, we don't need to normalize by the time bin of the light curve

		// Apply distance corrections
		// (*p_source_light_curve).lc_rate.at(i) /= 4 * M_PI * pow(lum_dist(z),2);
	}	
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load SynthGRB from file
void SynthGRB::LoadFromFile(std::string file_name)
{
	// Do some GRB loading
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to write out jet parameters found in the jet simulation
void SynthGRB::write_out_jet_params(std::string dir_path_name)
{
	// Write out thermal params
	ofstream thermal_params_file; // Construct file 
	thermal_params_file.open(dir_path_name+"synthGRB_jet_params_TH.txt"); // Open text file with this name
	size_t i=0;
	thermal_params_file << "# Thermal Dynamics" << endl;
	thermal_params_file << "# t_e (s) \t t_a (s) \t delt_a (s) \t T_phot (K) \t L_phot (erg/s) \t R_phot (cm/c) \t shell_index" << endl << endl;
	// For each time bin, write the time and count rate to the file.
	while ( i < te_therm.size())
	{
		thermal_params_file << te_therm.at(i);
		thermal_params_file << " \t ";
		thermal_params_file << ta_therm.at(i);
		thermal_params_file << " \t ";
		thermal_params_file << delt_therm.at(i);
		thermal_params_file << " \t ";
		thermal_params_file << T_phot.at(i);
		thermal_params_file << " \t ";
		thermal_params_file << L_phot.at(i);
		thermal_params_file << " \t ";
		thermal_params_file << r_phot.at(i);
		thermal_params_file << " \t ";	
		thermal_params_file << shell_ind_th.at(i);
		thermal_params_file << endl;	
		++i;
	}
	thermal_params_file.close(); // Close file

	// Write out internal shock params
	ofstream is_params_file; // Construct file 
	is_params_file.open(dir_path_name+"synthGRB_jet_params_IS.txt"); // Open text file with this name
	i=0;
	is_params_file << "# Internal Shock Dynamics" << endl;
	is_params_file << "# t_e (s) \t t_a (s) \t delt_a (s) \t B_eq (G) \t Gamma_e \t E_syn (keV) \t Gamma_r \t L_diss (erg/s) \t nu_c (Hz) \t nu_m (Hz) \t shell_index \t a_syn \t tau \t relvel" << endl << endl;
	// For each time bin, write the time and count rate to the file.
	while ( i < te_is.size())
	{
		is_params_file << te_is.at(i);
		is_params_file << " \t ";
		is_params_file << ta_is.at(i);
		is_params_file << " \t ";
		is_params_file << delt_is.at(i);
		is_params_file << " \t ";
		is_params_file << beq_is.at(i);
		is_params_file << " \t ";
		is_params_file << gamma_e_is.at(i);
		is_params_file << " \t ";
		is_params_file << esyn_is.at(i);
		is_params_file << " \t ";
		is_params_file << gamma_r_is.at(i);
		is_params_file << " \t ";
		is_params_file << e_diss_is.at(i);
		is_params_file << " \t ";
		is_params_file << nu_c_is.at(i);
		is_params_file << " \t ";
		is_params_file << nu_m_is.at(i);
		is_params_file << " \t ";
		is_params_file << shell_ind_is.at(i);
		is_params_file << " \t ";
		is_params_file << asyn_is.at(i);
		is_params_file << " \t ";
		is_params_file << tau.at(i);
		is_params_file << " \t ";
		is_params_file << relvel.at(i);
		is_params_file << endl;	
		++i;
	}
	is_params_file.close(); // Close file

	// Write out external shock params
	ofstream fs_params_file; // Construct file 
	fs_params_file.open(dir_path_name+"synthGRB_jet_params_FS.txt"); // Open text file with this name
	i=0;
	fs_params_file << "# Forward Shock Dynamics" << endl;
	fs_params_file << "# t_e (s) \t t_a (s) \t delt_a (s) \t B_eq (G) \t Gamma_e \t E_syn (keV) \t Gamma_r \t L_diss (erg/s) \t nu_c (Hz) \t nu_m (Hz) \t theta_open (rad) \t shell_ind" << endl << endl;
	// For each time bin, write the time and count rate to the file.
	while ( i < te_fs.size())
	{
		fs_params_file << te_fs.at(i);
		fs_params_file << " \t ";
		fs_params_file << ta_fs.at(i);
		fs_params_file << " \t ";
		fs_params_file << delt_fs.at(i);
		fs_params_file << " \t ";
		fs_params_file << beq_fs.at(i);
		fs_params_file << " \t ";
		fs_params_file << gamma_e_fs.at(i);
		fs_params_file << " \t ";
		fs_params_file << esyn_fs.at(i);
		fs_params_file << " \t ";
		fs_params_file << gamma_r_fs.at(i);
		fs_params_file << " \t ";
		fs_params_file << e_diss_fs.at(i);
		fs_params_file << " \t ";
		fs_params_file << nu_c_fs.at(i);
		fs_params_file << " \t ";
		fs_params_file << nu_m_fs.at(i);
		fs_params_file << " \t ";
		fs_params_file << theta_fs.at(i);
		fs_params_file << " \t ";
		fs_params_file << 0;
		// fs_params_file << " \t ";
		// fs_params_file << eps_star_fs.at(i);
		// fs_params_file << " \t ";
		// fs_params_file << rho_fs.at(i);
		fs_params_file << endl;	
		++i;
	}
	fs_params_file.close(); // Close file

	ofstream rs_params_file; // Construct file 
	rs_params_file.open(dir_path_name+"synthGRB_jet_params_RS.txt"); // Open text file with this name
	i=0;
	rs_params_file << "# Reverse Shock Dynamics" << endl;
	rs_params_file << "# t_e (s) \t t_a (s) \t delt_a (s) \t B_eq (G) \t Gamma_e \t E_syn (keV) \t Gamma_r \t L_diss (erg/s) \t nu_c (Hz) \t nu_m (Hz) \t shell_ind" << endl << endl;
	// For each time bin, write the time and count rate to the file.
	while ( i < te_rs.size())
	{
		rs_params_file << te_rs.at(i);
		rs_params_file << " \t ";
		rs_params_file << ta_rs.at(i);
		rs_params_file << " \t ";
		rs_params_file << delt_rs.at(i);
		rs_params_file << " \t ";
		rs_params_file << beq_rs.at(i);
		rs_params_file << " \t ";
		rs_params_file << gamma_e_rs.at(i);
		rs_params_file << " \t ";
		rs_params_file << esyn_rs.at(i);
		rs_params_file << " \t ";
		rs_params_file << gamma_r_rs.at(i);
		rs_params_file << " \t ";
		rs_params_file << e_diss_rs.at(i);
		rs_params_file << " \t ";
		rs_params_file << nu_c_rs.at(i);
		rs_params_file << " \t ";
		rs_params_file << nu_m_rs.at(i);
		rs_params_file << " \t ";
		rs_params_file << shell_ind_rs.at(i);
		// rs_params_file << " \t ";
		// rs_params_file << eps_star_rs.at(i);
		// rs_params_file << " \t ";
		// rs_params_file << rho_rs.at(i);
		rs_params_file << endl;	
		++i;
	}
	rs_params_file.close(); // Close file


	ofstream mass_params_file; // Construct file 
	mass_params_file.open(dir_path_name+"synthGRB_jet_params_MA.txt"); // Open text file with this name
	i=0;
	mass_params_file << "# Shell Mass Evolution" << endl;
	mass_params_file << "# t_e (s) \t t_a (s) \t Mass FS (g) \t Mass RS (g) \t Mass EJ (g)" << endl << endl;
	// For each time bin, write the time and count rate to the file.
	while ( i < mass_te.size())
	{
		mass_params_file << mass_te.at(i);
		mass_params_file << " \t ";
		mass_params_file << mass_ta.at(i);
		mass_params_file << " \t ";
		mass_params_file << mass_fs.at(i);
		mass_params_file << " \t ";
		mass_params_file << mass_rs.at(i);
		mass_params_file << " \t ";
		mass_params_file << mass_ej.at(i);
		mass_params_file << endl;	
		++i;
	}
	mass_params_file.close(); // Close file
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Save spectrum to text file
void SynthGRB::WriteSpectrumToTXT(std::string out_file_name)
{
	// Write the light curve to a text file
	ofstream spec_file; // Construct file 
	spec_file.open(out_file_name); // Open text file with this name

	spec_file << "# Spectrum Count Rate" << endl;
	spec_file << "# Energy bin (keV) \t Rate (c/s) \t Uncertainty (c/s)" << endl;

	size_t i=0;
	// For each time bin, write the time and count rate to the file.
	while ( i < (*p_source_spectrum).spectrum_rate.size())
	{
		spec_file << (*p_source_spectrum).energ_mid.at(i);
		spec_file << " ";
		spec_file << (*p_source_spectrum).spectrum_rate.at(i);
		spec_file << " ";
		spec_file << (*p_source_spectrum).spectrum_unc.at(i);
		spec_file << "\n";
		++i;
	}
	spec_file.close(); // Close file
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Save spectrum to FITS file
void SynthGRB::WriteSpectrumToFITS(std::string out_file_name)
{

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Save light curve to text file
void SynthGRB::WriteLightCurveToTXT(std::string out_file_name)
{
	// Write the light curve to a text file
	ofstream light_curve_file; // Construct file 
	light_curve_file.open(out_file_name); // Open text file with this name

	light_curve_file << "# Light Curve Count Rate" << endl;
	light_curve_file << "# Time (sec) \t Rate (c/s)" << endl;

	size_t i=0;
	// For each time bin, write the time and count rate to the file.
	while ( i < (*p_source_light_curve).lc_rate.size())
	{
		light_curve_file << (*p_source_light_curve).lc_time.at(i);
		light_curve_file << " ";
		light_curve_file << (*p_source_light_curve).lc_rate.at(i);
		light_curve_file << "\n";		
		++i;
	}
	light_curve_file.close(); // Close file
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Save light curve to FITS file
void SynthGRB::WriteLightCurveToFITS(std::string out_file_name)
{

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Write out file to be read by the amazing afterglow calculation code, Boulodrome. 
void SynthGRB::WriteBoulodromeTXT(std::string out_file_name, float te_start)
{
	ofstream boul_file; // Construct file 
	
	boul_file.open(out_file_name); // Open text file with this name

	boul_file << "# Boulodrome Data File" << endl;
	boul_file << "# Radius (light second) \t Time (sec) \t Gamma \t Beta \t Num Density (cm^-3) \t Specific Energy (erg cm^-3) \t Num. Acc. e- \t Jet Opening Angle (rad)" << endl;

	size_t j=0;
	size_t i=0;

	// Find the index where the FS emits at a time = tstart
	while ( j < te_fs.size())
	{
		if(te_fs.at(j) >= te_start)
		{
			i = j;
			j = te_fs.size()+1;
		}
		else{ ++j; }
	}

	// For each time bin, write the time and count rate to the file.
	while ( i < te_fs.size())
	{
		boul_file << rad_coll_fs.at(i); // light seconds, Radius of forward shell
		boul_file << " \t ";
		boul_file << te_fs.at(i); // sec, Emission time for this information
		boul_file << " \t ";
		boul_file << gamma_r_fs.at(i); // Lorentz factor
		boul_file << " \t ";
		boul_file << beta(gamma_r_fs.at(i)); // v/c
		boul_file << " \t ";
		boul_file << rho_fs.at(i)*pow(rad_coll_fs.at(i),(*p_model_params).k_med)/gamma_r_fs.at(i)/mp; // g cm^-3, Mass density
		boul_file << " \t ";
		boul_file << eps_star_fs.at(i); // erg cm^-3, Specific energy density
		boul_file << " \t ";
		boul_file << num_swept_e_fs.at(i); // Number of swept up electrons in the forward shell
		boul_file << " \t ";
		boul_file << theta_fs.at(i); // rad, Jet opening angle
		boul_file << "\n";
		++i;
	}
	boul_file.close(); // Close file
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Private Member Functions */

// Find the resulting Lorentz factor of two collided shells from the Lorentz factor and Mass of shell 1 and 2
double SynthGRB::_calc_shell_coll_gamma(float s1g, float s2g, float s1m, float s2m)
{
	// The approximate Lorentz factor of the combined layer after complete redistribution of momentum and energy
	return sqrt( s1g*s2g * (s1m*s1g + s2m*s2g) / (s1m*s2g + s2m*s1g)  );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Calculate relative velocity between two shells with Lorentz factors gamma_1 and gamma_2, where gamma_1 > gamma_2
double SynthGRB::_rel_vel(float gamma_1, float gamma_2)
{
	return (pow(gamma_1,2.) - pow(gamma_2,2.)) / (pow(gamma_1,2.) + pow(gamma_2,2.));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Check if shells are ordered
int SynthGRB::_check_if_sorted()
{
	// Make list of shells which are currently active:
	std::vector<int> active_inds;
	for(int i=0; i<(*p_model_params).numshells; ++i)
	{
		if((*p_jet_shells).shell_status.at(i) == 1)
		{
			active_inds.push_back(i);
		}
	}

	// Check if any shells are active
	if(active_inds.size() == 0)
	{
		return 1; // True, if there are no active shells, then its sorted, right?
	}

	float tmp_diff = 0;
	int shell_ds = 0.;
	int shell_us = 0.;
	for(size_t i=0; i<active_inds.size()-1; ++i)
	{
		shell_ds = active_inds.at(i);
		shell_us = active_inds.at(i+1);

		tmp_diff = (*p_jet_shells).shell_gamma.at(shell_us) - (*p_jet_shells).shell_gamma.at(shell_ds);
		if(tmp_diff > 0) // Then the upstream shell has a higher lorentz factor than the down stream shell (so a collision will occur)
		{
			return 0; // False, the Lorentz distribution is not sorted and there will still be a collision that will occur
		}
	
	}
	return 1; // True, the Lorentz distribution is sorted and there will no longer be any collisions
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make the energy axis of a spectrum, bounded by Emin and Emax with a number of bins = num_en_bins
int SynthGRB::_make_en_axis(std::vector<float> & energy_axis, float emin, float emax, int num_en_bins)
{
	// Move to log space to define the energy interval with equally spaced points (in log space)
	float log_emin = log10(emin);
	float log_emax = log10(emax);
	float log_de = (log_emax - log_emin) / num_en_bins;

	// For each bin, calculate the energy axis
	for(int i=0;i<num_en_bins;++i)
	{
		energy_axis.at(i) = pow( 10., log_emin + (i*log_de));
	}

	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Used to calculate the fraction of an emission time interval overlaps with a time bin
float SynthGRB::_fraction_of_interval_in_time_bin(float t_a, float del_t_a, float tbin_min, float tbin_max )
{
	float profile_factor = 1.;

	// Case 1: If the emission starts before the time bin and ends after this time bin
	if( (t_a < tbin_min ) and ( (t_a + del_t_a) > tbin_max) )
	{
		profile_factor = (tbin_max - tbin_min) / del_t_a;
	}
	// Case 2: If the emission starts before the time bin and ends within the time bin
	else if(t_a < tbin_min)
	{
		profile_factor = ( (t_a + del_t_a) - tbin_min) / del_t_a;
	}
	// Case 3: If the emission starts within the time bin and ends after the time bin
	else if( (t_a+del_t_a) > tbin_max)
	{
		profile_factor = (tbin_max - t_a) / del_t_a;
	}
	// Case 4: If the emission starts and ends within the time bin, then the fraction is = 1. 

	return profile_factor;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Calculates the integral under the pulse profile curve integrated between t_1 and t_2 (see Equation 5 in Genet et al 2007)

double SynthGRB::_calc_pulse_profile_factor(float t_a, float del_t_a, float t_low, float t_hi, float theta, float Rsh)
{
	if(t_low < t_a)
	{
		t_low = t_a;
	}

	if(t_hi > (t_a +  ( (1. - cos(theta))*Rsh)) )
	{
		t_hi = (t_a +  ( (1. - cos(theta))*Rsh));
	}

	if(t_low > (t_a +  ( (1. - cos(theta))*Rsh)) )
	{
		return 0; 
	}

	// Evaluating lower bound 
	double bound_low = 1./pow(t_low - t_a + del_t_a,2.);

	// Evaluating upper bound
	double bound_hi = 1./pow(t_hi - t_a + del_t_a,2.);

	// Analytical integral of Equation 5 in Genet et al 2007
	double pulse_profile_factor = - pow(del_t_a,3.) * (bound_hi - bound_low);

	return pulse_profile_factor;
}