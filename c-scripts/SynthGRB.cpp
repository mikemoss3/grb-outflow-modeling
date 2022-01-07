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

	// Assign default values for the jet parameters
	/*
	(*p_model_params).tw = 10; // sec, Duration of the wind
	(*p_model_params).dte = 0.002; // sec, Time between successive shell launches
	(*p_model_params).eps_e = 1./3.; // Fraction of dissipated energy that goes into the electrons 
	(*p_model_params).eps_b = 1./3.; // Fraction of dissipated energy that goes into the magnetic field 
	(*p_model_params).zeta = 1e-3; // Fraction of electrons which are accelerated 
	(*p_model_params).E_dot_iso = 1e53; // erg/s, Injected isotropic energy rate
	(*p_model_params).theta = 0.1; // radians, Half-opening angle of the jet
	(*p_model_params).r_open = 1e6; // cm, Opening radius of the jet
	(*p_model_params).eps_th = 0.03; // Fraction of energy in the outflow in the form of thermal energy 
	(*p_model_params).sigma = 0.1; // Magnetization of the outflow 
	// (*p_model_params).LorentzDist = "step"; // Distribution of the jet shells
	(*p_model_params).LorentzDist = "oscillatory"; // Distribution of the jet shells
	(*p_model_params).ShellDistParamsFile = "Default"; // File that contains the parameters to create the distribution of jet shells. If Default is used, default params are passed
	*/ 

	// Initialize jet based on current jet parameters and shell distribution
	InitializeJet();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Overloaded SynthGRB constructor, uses given arguments 
SynthGRB::SynthGRB(float tw, float dte, float eps_e, float eps_b, float zeta, double E_dot_iso, float theta, float r_open, float eps_th, float sigma, std::string LorentzDist, std::string ShellDistParamsFile)
{
	// Assign pointers to NULL
	p_source_spectrum = NULL;
	p_source_light_curve = NULL;
	p_jet_shells = NULL;
	p_model_params = new ModelParams(tw, dte, eps_e, eps_b, zeta, E_dot_iso, theta, r_open, eps_th, sigma, LorentzDist, ShellDistParamsFile);

	// Assign default values for the jet parameters
	/*
	(*p_model_params).tw = 10; // sec, Duration of the wind
	(*p_model_params).dte = 0.002; // sec, Time between successive shell launches
	(*p_model_params).eps_e = eps_e; // Fraction of dissipated energy that goes into the electrons 
	(*p_model_params).eps_b = eps_b; // Fraction of dissipated energy that goes into the magnetic field 
	(*p_model_params).zeta = zeta; // Fraction of electrons which are accelerated 
	(*p_model_params).E_dot_iso = E_dot_iso; // erg/s, Injected isotropic energy rate
	(*p_model_params).theta = theta; // radians, Half-opening angle of the jet
	(*p_model_params).r_open = r_open; // cm, Opening radius of the jet
	(*p_model_params).eps_th = eps_th; // Fraction of energy in the outflow in the form of thermal energy 
	(*p_model_params).sigma = sigma; // Magnetization of the outflow 
	// (*p_model_params).LorentzDist = LorentzDist; // Distribution of the jet shells
	(*p_model_params).LorentzDist = LorentzDist; // Distribution of the jet shells
	(*p_model_params).ShellDistParamsFile = ShellDistParamsFile; // File that contains the parameters to create the distribution of jet shells. If Default is used, default params are passed
	*/ 

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
	m_bar = (*p_model_params).E_dot_iso*(*p_model_params).dte/gamma_bar/pow(c_cm,2.);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load jet parameters from a file
void SynthGRB::LoadJetParamsFromTXT(std::string file_name)
{
	// Load in the file
	ifstream file_jet_params(file_name);
	string line_jet_params;
	if ( file_jet_params.is_open() ) 
	{
		// tw
		getline( file_jet_params, line_jet_params);
		float col1_val = stof(line_jet_params);
		// dte
		getline( file_jet_params, line_jet_params);
		float col2_val = stof(line_jet_params);	
		// eps_e
		getline( file_jet_params, line_jet_params);
		float col3_val = stof(line_jet_params);
		// eps_b
		getline( file_jet_params, line_jet_params);
		float col4_val = stof(line_jet_params);
		// zeta
		getline( file_jet_params, line_jet_params);
		float col5_val = stof(line_jet_params);
		// E_dot_iso
		getline( file_jet_params, line_jet_params);
		double col6_val = stod(line_jet_params);
		// theta
		getline( file_jet_params, line_jet_params);
		float col7_val = stof(line_jet_params);
		// r_open
		getline( file_jet_params, line_jet_params);
		float col8_val = stof(line_jet_params);
		// eps_th
		getline( file_jet_params, line_jet_params);
		float col9_val = stof(line_jet_params);
		// sigma
		getline( file_jet_params, line_jet_params);
		float col10_val = stof(line_jet_params);
		// LorentzDist
		getline( file_jet_params, line_jet_params);
		string col11_val = line_jet_params;
		// ShellDistParamsFile
		getline( file_jet_params, line_jet_params);
		string col12_val = line_jet_params;
		
		p_model_params = new ModelParams(
			col1_val, 
			col2_val, 
			col3_val, 
			col4_val, 
			col5_val, 
			col6_val, 
			col7_val, 
			col8_val, 
			col9_val, 
			col10_val,
			col11_val, 
			col12_val);
		
		// Initialize jet based on current jet parameters and shell distribution
		InitializeJet();

		// Close files and free memory 
		file_jet_params.close();
	}
	else std::cout << "Unable to open file.";


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
			float g1, g2, mfrac;
			bool fluctuations;
			// Load from designated file
			ifstream file_jet_params((*p_model_params).ShellDistParamsFile);
			string line_jet_params;
			if ( file_jet_params.is_open() ) 
			{

				getline( file_jet_params, line_jet_params);
				g1 = stof(line_jet_params);

				getline( file_jet_params, line_jet_params);
				g2 = stof(line_jet_params);
				
				getline( file_jet_params, line_jet_params);
				mfrac = stof(line_jet_params);
				
				getline( file_jet_params, line_jet_params);
				istringstream(line_jet_params) >> fluctuations;


				// Make shell distribution with input parameters
				(*p_jet_shells).step((*p_model_params).dte, g1, g2, mfrac, fluctuations);	

				// Close files and free memory 
				file_jet_params.close();
			}
			else std::cout << "Unable to open file.";
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
			float median, amp, freq, decay;
			bool fluctuations;
			// Load from designated file
			ifstream file_jet_params((*p_model_params).ShellDistParamsFile);
			string line_jet_params;
			if ( file_jet_params.is_open() ) 
			{
				getline( file_jet_params, line_jet_params);
				median = stof(line_jet_params);

				getline( file_jet_params, line_jet_params);
				amp = stof(line_jet_params);

				getline( file_jet_params, line_jet_params);
				freq = stof(line_jet_params);

				getline( file_jet_params, line_jet_params);
				decay = stof(line_jet_params);

				getline( file_jet_params, line_jet_params);
				istringstream(line_jet_params) >> fluctuations;

				// Make shell distribution with input parameters
				(*p_jet_shells).oscillatory((*p_model_params).dte, median, amp, freq, decay, fluctuations);	

				// Close files and free memory 
				file_jet_params.close();
			}
			else std::cout << "Unable to open file.";
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

	// Reset Synchrotron variable arrays 
	te_synch.resize(0);
	ta_synch.resize(0);
	delt_synch.resize(0);
	asyn.resize(0);
	beq.resize(0);
	gamma_e.resize(0);
	esyn.resize(0);
	gamma_r.resize(0);
	e_diss.resize(0);
	tau.resize(0);
	relvel.resize(0);

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

	for(int i=0; i<(*p_model_params).numshells; ++i)
	{
		// Calculate the photospheric radius for each jet shell, Equation 9 of Hascoet 2013
		r_phot.at(i) = (0.2*(*p_model_params).E_dot_iso / ( (1+(*p_model_params).sigma)*8.*M_PI*pow(c_cm,4.) *pow( (*p_jet_shells).shell_gamma.at(i) ,3.) ) ); // units of lightseconds, cm/speed of light
		// r_phot = 2.9*pow(10,13) *((*p_model_params).E_dot_iso/1e53) / ( cc.c*(1+(*p_model_params).sigma)* pow((*p_jet_shells).shell_gamma/100,3) ) # units of lightseconds, cm/(cm/s)

		// Times when each shell will cross the photosphere
		te_therm.at(i) = (r_phot.at(i) - (*p_jet_shells).shell_radius.at(i) ) / beta((*p_jet_shells).shell_gamma.at(i));
		ta_therm.at(i) = (te_therm.at(i) - r_phot.at(i)); // Time emission arrives at observer 

		delt_therm.at(i) = r_phot.at(i)/2./pow((*p_jet_shells).shell_gamma.at(i),2.); // Duration width of the emission 
		// Calculate useful constant for next calculations
		Phi = pow((*p_model_params).theta,-2./3.) * pow(r_phot.at(i)*c_cm,-2./3.) * pow((*p_model_params).r_open,2./3.) * pow((*p_jet_shells).shell_gamma.at(i),2./3.); // Eq. 11 in Hascoet 2013

		// Temperature at photosphere
		T0 = pow((*p_model_params).E_dot_iso*pow((*p_model_params).theta,2.) / (4*M_PI * a * c_cm * pow((*p_model_params).r_open,2.) ),1./4.); // K, rearrangement Eq. 1 of Hascoet 2013 
		// T0 =  (2/3)*np.pow((*p_model_params).eps_th,1/4)*np.pow((*p_model_params).theta/0.1,1/2)*np.pow((*p_model_params).E_dot_iso/1e53,1/4)*np.pow((*p_model_params).r_open/1e7,-1/2) / (cc.kb_kev/1000) # K, alternative expression, Eq. in Hascoet

		T_phot.at(i) = T0*Phi; // K, Equation 7 in Hascoet 2013

		// Luminosity at photosphere, from equation 8 in Hascoet 2013
		L_phot.at(i) = (pow((*p_model_params).theta,2.) / 4.) * E_dot_therm * Phi; // erg/s, beamed
		// L_phot.at(i) = E_dot_therm * Phi; // erg/s, isotropic
		// L_phot.at(i) = pow((*p_jet_shells).shell_gamma.at(i),2.) * a * pow( T0*Phi/(*p_jet_shells).shell_gamma.at(i),4.) * c_cm * (M_PI * pow((*p_model_params).theta,2.) * pow(c_cm * r_phot.at(i),2) ); // erg/s, alternative expression
	}

	/* Calculate Internal and External shock considerations */
	// Set flag to mark if all the Lorentz factors of the jet shells are in order
	bool ord_lorentz = false;
	
	// Initialize variables for internal shocks
	int ind_ds_shell_tmp = 0; // Temporary index of the down stream shell
	int ind_us_shell_tmp = 0; // Temporary index of the upstream shell
	double tmp_ds_r = 0.; // Temporary Radius of the down stream shell
	double tmp_ds_b = 0.; // Temporary Beta (i.e., v/c) of the down stream shell
	double tmp_us_r = 0.; // Temporary Radius of the up stream shell
	double tmp_us_b = 0.; // Temporary Beta of the upstream shell
	double t_coll_tmp = 0.; // Time until collision of the upstream and downstream shells
	
	int ind_s_ds = 0; // Index of the down stream shell which will collide next
	int ind_s_us = 0; // Index of the up stream shell which will collide next
	double t_coll_lowest = 0.; // Mini(*p_model_params).mum time until next shell collision

	double shell_ds_m = 0.; // Downstream shell mass 
	double shell_ds_g = 0.; // Downstream shell Lorentz factor 
	double shell_us_m = 0.; // Upstream shell mass 
	double shell_us_g = 0.; // Upstream shell Lorentz factor 

	double gamma_comb = 0.; // Final Lorentz factor of the combined shell (after complete redistribution of momenta)
	double rad_coll = 0.; // Radius of the collision
	double gamma_int = 0.; // Lorentz factor for internal motion in shocked material
	double eps = 0.; // erg, Average proton factor from the collision of two shells
	double n = 0.; // 1/cm^3, Comoving proton density
	double t_syn = 0.; // Synchrotron time-scale
	double quadr_const = 0.; // Constant in the quadratic equation
	double w = 0.; // Indicates whether we are in the Klein-Nishina regime or Thompson (i.e., w>>1 or w<<1, respectively)
	double Q_IC = 0.; // Inverse Compton parameter
	double alpha_ic = 0.; // Fraction of energy that goes into Inverse Compton electrons
	
	double tmp_te = 0.; // Temp values
	double tmp_asyn = 0.; // Temp values
	double tmp_beq = 0.; // Temp values
	double tmp_gamma_e = 0.; // Temp values
	double tmp_esyn_kev = 0.; // Temp values
	double tmp_esyn_erg = 0.; // Temp values
	double tmp_gamma_r = 0.; // Temp values
	double tmp_e_diss = 0.; // Temp values
	double tmp_tau = 0.; // Temp values

	// Initialize variable relevant for external shocks
	bool external_shock = false; // Flag to indicate whether external shock has started or not
	bool all_inactve = false; // Flag to indicate whether all internal shells have been passed by the reverse shock.
	float R_ext0 = 5*pow(10,15); // cm, radius at which the external medium becomes relevant  
	float k = 2.;// Wind density profile for external medium
	// float k = 0.;// Constant density profile for external medium
	float n0 = 0.1; // 1/cm^3, particle number density at R_ext0
	float A = mp*n0*pow(R_ext0,k); // Normalization of the density
	float q = 0.01; // Ratio which dictates how much external mass needs to be swept up in order for an internal show to occur
	float m_ex = 0; // g, Mass limit for a forward shock to occur 
	float FS_t_sweep = 0.; // sec, Time until a mass m_ex is swept up by the FS   
	float RS_t_coll = 0.; // sec, Time until a RS collision
	float gamma_r_ext = 0; // Resulting Lorentz factor from a FS or RS

	float FS_shell_mass = 0.; // g, FS shell mass
	float FS_shell_gamma = 0.; // FS shell Lorentz factor
	float FS_shell_gamma_int = 0.; // FS shell Lorentz factor
	float FS_shell_radius = 0.; // cm, FS shell position
	float RS_shell_mass = 0.; // g, RS shell mass
	float RS_shell_gamma = 0.; // RS shell Lorentz factor
	float RS_shell_radius = 0.; // cm, RS shell position

	while( (ord_lorentz==false) & (all_inactve==false) )
	{
		/* Calculate the time of collision between all adjacent (and active) shells */
		// Initialize a high number as the mini(*p_model_params).mum collision time
		t_coll_lowest = std::numeric_limits<double>::max();

		// Make list of shells which are currently active:
		std::vector<int> active_inds;
		for(int i=0; i<(*p_model_params).numshells; ++i)
		{
			if((*p_jet_shells).shell_status.at(i) == 1)
			{
				active_inds.push_back(i);
			}
		}

		// For each active shell
		for(size_t i=0; i<active_inds.size()-1; ++i)
		{
			// The current active shell is the down stream shell
			ind_ds_shell_tmp = active_inds.at(i);
			// The next active shell is the up stream shell
			ind_us_shell_tmp = active_inds.at(i+1);

			// The down stream shell (*p_model_params).must have a Lorentz factor smaller than the upstream shell or they will never collide
			if( (*p_jet_shells).shell_gamma.at(ind_ds_shell_tmp) < (*p_jet_shells).shell_gamma.at(ind_us_shell_tmp) )
			{
				// The down stream shell is farther out, but slower
				// The up stream shell is closer to the central engine, but faster
				tmp_ds_r = (*p_jet_shells).shell_radius.at(ind_ds_shell_tmp);
				tmp_ds_b = beta( (*p_jet_shells).shell_gamma.at(ind_ds_shell_tmp) );
				tmp_us_r = (*p_jet_shells).shell_radius.at(ind_us_shell_tmp);
				tmp_us_b = beta( (*p_jet_shells).shell_gamma.at(ind_us_shell_tmp) );

				// Calculate the time until these shells collide
				t_coll_tmp = (tmp_ds_r - tmp_us_r) / (tmp_us_b - tmp_ds_b);

				// Check if this time is shorter than the previously known shortest time until collision
				if(t_coll_tmp < t_coll_lowest)
				{
					t_coll_lowest = t_coll_tmp;

					// Record which shells these were
					ind_s_ds = ind_ds_shell_tmp;
					ind_s_us = ind_us_shell_tmp;
				}
			}
		} 

		// If the external shock has not started yet, check if any shell has reached the external medium, begin the external shock
		if( (external_shock == false) & ((*p_jet_shells).shell_radius.at(active_inds.at(0)) >= R_ext0))
		{
			external_shock = true;

			FS_shell_mass = (*p_jet_shells).shell_mass.at(active_inds.at(0));
			FS_shell_gamma = (*p_jet_shells).shell_gamma.at(active_inds.at(0));
			FS_shell_radius = (*p_jet_shells).shell_radius.at(active_inds.at(0));
			RS_shell_mass = (*p_jet_shells).shell_mass.at(active_inds.at(0));
			// RS_shell_gamma = (*p_jet_shells).shell_gamma.at(active_inds.at(0)); // Same as FS
			// RS_shell_radius = (*p_jet_shells).shell_radius.at(active_inds.at(0)); // Same as FS

		}
		else
		{
			/* Calculate time until next FS */
			// Calculate when the wind has swept up external medium with mass m_ex = q*M/Gamma, e.g., when a forward shock will occur
			m_ex = q*(FS_shell_mass + RS_shell_mass)/FS_shell_gamma;
			if (k == 2)
			{
				// Wind density profile for the medium:
				FS_t_sweep = beta(FS_shell_gamma)*c_cm * q*(FS_shell_mass + RS_shell_mass) / (4 * M_PI * A * FS_shell_gamma);	
			}
			else if (k==0)
			{
				// Constant density profile for the medium:
				FS_t_sweep = beta(FS_shell_gamma)*c_cm *  ( pow(3*q*(FS_shell_mass + RS_shell_mass) / (4 * M_PI * A * FS_shell_gamma) + pow(FS_shell_radius,3),1/3) - FS_shell_radius);
			}

			/* Calculate time until next RS */
			// Calculate when a wind ejecta shell will collide with the discontinuity, e.g., a reverse shock
			// Select farthest active shell
			float outer_shell_radius = (*p_jet_shells).shell_radius.at(active_inds.at(0));
			float outer_shell_gamma = (*p_jet_shells).shell_gamma.at(active_inds.at(0));
			RS_t_coll = (FS_shell_radius - outer_shell_radius) / (beta(outer_shell_gamma) - beta(FS_shell_gamma));
		}
		
		// If t_coll_lowest < FS_t_sweep & < RS_t_coll: calculate internal shock contribution
		if((t_coll_lowest < FS_t_sweep) & (t_coll_lowest < RS_t_coll))
		{
			// Move to the emission time of the shell collision (in the rest frame of the jet)
			tmp_te += t_coll_lowest;
			// Move shells forward
			for(int i=0; i<(*p_model_params).numshells; ++i)
			{
				if((*p_jet_shells).shell_status.at(i) == 1)
				{
					(*p_jet_shells).shell_radius.at(i) += beta((*p_jet_shells).shell_gamma.at(i)) * t_coll_lowest;
				}
			}

			// For ease of reading, make variables for the masses and Lorentz factors of the two colliding shells
			// Down stream shell
			shell_ds_m = (*p_jet_shells).shell_mass.at(ind_s_ds);
			shell_ds_g = (*p_jet_shells).shell_gamma.at(ind_s_ds);
			// Upstream shell
			shell_us_m = (*p_jet_shells).shell_mass.at(ind_s_us);
			shell_us_g = (*p_jet_shells).shell_gamma.at(ind_s_us);

			// Approximate the resulting Lorentz factor from the collision
			tmp_gamma_r = sqrt(shell_ds_g*shell_us_g);
			// Final Lorentz factor of the combined shell (after complete redistribution of momenta)
			gamma_comb = _calc_shell_coll_gamma(shell_ds_g, shell_us_g, shell_ds_m, shell_us_m);

			// Set the Lorentz factor of the merged shells
			(*p_jet_shells).shell_gamma.at(ind_s_us) = gamma_comb;
			// Set the mass of the merged shell (i.e., add the mass of the two shells)
			(*p_jet_shells).shell_mass.at(ind_s_us) += (*p_jet_shells).shell_mass.at(ind_s_ds);
			
			// Record the radius of the collision
			rad_coll = (*p_jet_shells).shell_radius.at(ind_s_us);

			// De-active the downstream shell now that it has been merged
			(*p_jet_shells).shell_status.at(ind_s_ds) = 0;

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

			gamma_int = 0.5 *( sqrt(shell_ds_g/shell_us_g) + sqrt(shell_us_g/shell_ds_g) ); // Lorentz factor for internal motion in shocked material
			eps = (gamma_int - 1) *mp * pow(c_cm,2.); // erg, Average proton factor from the collision of two shells

			// Characteristic Lorentz factor of accelerated electrons
			// tmp_gamma_e = (*p_model_params).eps_e * eps / (me * pow(c_cm,2) ); // Assuming electron are in equipartition with protons
			// tmp_gamma_e =  pow( ((*p_model_params).alpha_m/(*p_model_params).zeta) * (eps/me/pow(c_cm,2)) ,1/(3-(*p_model_params).mu)); // For electrons directly producing synchrotron radiation (via turbulent magnetic fields)
			tmp_gamma_e = 1e4; // Constant

			n = E_dot_kin / (4. * M_PI * mp * pow(rad_coll*gamma_bar,2.) * pow(c_cm,5.)); // 1/cm^3, Comoving proton density
			tmp_beq = sqrt( 8.* M_PI * (*p_model_params).eps_b * n * eps); // (erg/cm^3)^1/2, Magnetic field density, assuming B = Beq  
			tmp_esyn_kev = 50.*(tmp_gamma_r/300.)*(tmp_beq/1e3)*pow((tmp_gamma_e/100.),2.) / 1000.; // keV, Synchrotron energy in the rest frame
			tmp_esyn_erg = kev_to_erg * tmp_esyn_kev; // erg, Synchrotron energy in the rest frame
			
			w = tmp_gamma_e*(tmp_esyn_erg/tmp_gamma_r)/(me * pow(c_cm,2.) ); // Indicates whether we are in the Klein-Nishina regime or Thompson (i.e., w>>1 or w<<1, respectively)
			// w = 33.*(tmp_beq/1e3)*pow((tmp_gamma_e/1e4),3.); // Alternative relation to find w 

			t_syn = 6.*pow(tmp_gamma_e/100.,-1.)*pow(tmp_beq/1000.,-2.); // Synchrotron time-scale
			quadr_const = -(3./2.)*(0.2/M_PI/pow(c_cm*rad_coll,2.)) * tmp_gamma_e * 100. * (E_dot_kin/pow(gamma_bar*c_cm,2.)) * pow(tmp_beq/1000.,-2.); // Constant in the quadratic equation
			// Calculate Compton parameter
			Q_IC = (-1 + sqrt(1-4.*quadr_const))/2.;
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

			// Calculate the energy dissipated in this collision, using the smaller mass of the two shells:
			if(shell_us_m > shell_ds_m)
			{
				tmp_e_diss = shell_ds_m * m_bar * pow(c_cm,2.) * (shell_ds_g + shell_us_g - 2.*tmp_gamma_r) * (*p_model_params).eps_e * tmp_asyn / (rad_coll/2./pow(tmp_gamma_r,2.)); // erg / s
				// tmp_e_diss = shell_ds_m * m_bar * pow(c_cm,2.) * (shell_ds_g + shell_us_g - 2.*tmp_gamma_r) * (*p_model_params).eps_e * tmp_asyn; // erg
			}
			else
			{
				tmp_e_diss = shell_us_m * m_bar * pow(c_cm,2.) * (shell_ds_g + shell_us_g - 2.*tmp_gamma_r) * (*p_model_params).eps_e * tmp_asyn / (rad_coll/2./pow(tmp_gamma_r,2.)); // erg / s
				// tmp_e_diss = shell_us_m * m_bar * pow(c_cm,2.) * (shell_ds_g + shell_us_g - 2.*tmp_gamma_r) * (*p_model_params).eps_e * tmp_asyn; // erg
			}

			// If the emission is efficient, add the contribution
			// E.g., if the emission time is less than the shell expansion (the dynamical time scale of the shell)		
			if ( t_syn < ((1+Q_IC)*rad_coll/gamma_bar) )
			{
				// Calculate the arrival time at the observer and the width of the emission 
				te_synch.push_back(tmp_te); // sec, Time of emission (in the rest frame of the jet), is also equal to the time a shell crosses the photosphere
				ta_synch.push_back(tmp_te - rad_coll); // sec, Time when the emission arrives at the observer (in the observer frame)
				delt_synch.push_back(rad_coll/2./pow(tmp_gamma_r,2.)); // sec, Width of the emission (in observer frame)
				asyn.push_back(tmp_asyn); // Fraction of the energy in electrons which goes into synchrotron (as opposed to Inverse Compton)
				beq.push_back(tmp_beq); // (erg/cm^3)^1/2, Magnetic field density, assuming B = Beq 
				gamma_e.push_back(tmp_gamma_e); // Lorentz factor of the electron population
				esyn.push_back(tmp_esyn_kev); // keV, Synchrotron energy emitted by accelerated electron
				gamma_r.push_back(tmp_gamma_r); // Approximation of the combined Lorentz factor of the colliding shells
				e_diss.push_back(tmp_e_diss); // erg, Dissipated energy during the collision 
				tau.push_back(tmp_tau); // Optical depth at the location of the collision
				relvel.push_back( _rel_vel(shell_us_g,shell_ds_g) ); // Relative velocity between the two shells
			}
		}
		// If RS_t_coll < t_coll_lowest & < FS_t_sweep: calculate RS contribution
		else if( (RS_t_coll < t_coll_lowest) & (RS_t_coll < FS_t_sweep))
		{
			// Move to the emission time of the shell collision (in the rest frame of the jet)
			tmp_te += RS_t_coll;
			// Move shells forward
			for(int i=0; i<(*p_model_params).numshells; ++i)
			{
				if((*p_jet_shells).shell_status.at(i) == 1)
				{
					(*p_jet_shells).shell_radius.at(i) += beta((*p_jet_shells).shell_gamma.at(i)) * RS_t_coll;
				}
			}
			// Move FS and RS forward
			FS_shell_radius += beta(FS_shell_gamma)*RS_t_coll;
			RS_shell_radius += beta(RS_shell_gamma)*RS_t_coll;

			/*
			Reverse Shock
			*/ 
			// Calculate new Lorentz factor at the discontinuity due to reverse shock 
			float m_ej = (*p_jet_shells).shell_gamma.at(active_inds.at(0)); // The mass of the shell that is colliding with the reverse shock 
			float gamma_ej = (*p_jet_shells).shell_gamma.at(active_inds.at(0)); // The lorentz factor of the shell that is colliding with the reverse shock 
			gamma_r_ext = pow(FS_shell_gamma*gamma_ej,0.5)*pow( ((RS_shell_mass + FS_shell_mass*FS_shell_gamma_int )*FS_shell_gamma + m_ej*gamma_ej )/((RS_shell_mass + FS_shell_mass*FS_shell_gamma_int )*gamma_ej + m_ej*FS_shell_gamma ),0.5);

			// Record the emission contribution and update the RS shell
			// Calculate energy dissipated by reverse shock
			RS_e_diss.push_back( pow(c_cm,2) * ( (RS_shell_mass + FS_shell_mass*FS_shell_gamma_int)*FS_shell_gamma + m_ej*gamma_ej - (RS_shell_mass + m_ej+ FS_shell_mass*FS_shell_gamma_int )*FS_shell_gamma) );

			// Update the mass of the revers shock shell
			RS_shell_mass += m_ej ;
			// Update the Lorentz factor of the forward and reverse shock shells
			FS_shell_gamma = gamma_r_ext;

			// Record the time of the emission
			RS_te.push_back(tmp_te);
			RS_ta.push_back(tmp_te - RS_shell_radius);
			RS_delt.push_back(RS_shell_radius/2./pow(gamma_r_ext,2.));
			// Record the new Lorentz factor
			RS_gamma_r.push_back(gamma_r_ext);

			// Deactivate internal shell that collided with RS
			(*p_jet_shells).shell_status.at(active_inds.at(0)) = 0;
		}
		// else (FS is shortest time) calculate FS contribution
		else
		{
			// Move to the emission time of the shell collision (in the rest frame of the jet)
			tmp_te += FS_t_sweep;
			// Move shells forward
			for(int i=0; i<(*p_model_params).numshells; ++i)
			{
				if((*p_jet_shells).shell_status.at(i) == 1)
				{
					(*p_jet_shells).shell_radius.at(i) += beta((*p_jet_shells).shell_gamma.at(i)) * FS_t_sweep;
				}
			}
			// Move FS and RS forward
			FS_shell_radius += beta(FS_shell_gamma)*FS_t_sweep;
			RS_shell_radius += beta(RS_shell_gamma)*FS_t_sweep;

			/*
			Forward  Shock
			*/ 
			// Calculate new Lorentz factor at the discontinuity due to forward shock
			gamma_r_ext = pow( ((FS_shell_mass + RS_shell_mass)*pow(FS_shell_gamma,2) + m_ex*FS_shell_gamma ) / ( (RS_shell_mass + FS_shell_mass*FS_shell_gamma_int) + 2*m_ex*FS_shell_gamma ) ,0.5);

			// Calculate new internal Lorentz factor in the shocked forward medium due to forward shock 
			FS_shell_gamma_int = ((RS_shell_mass + FS_shell_mass*FS_shell_gamma_int)*FS_shell_gamma + m_ex - RS_shell_mass*gamma_r_ext)/(FS_shell_mass+m_ex)/gamma_r_ext ;

			// Calculate energy dissipated by forward shock
			FS_e_diss.push_back( gamma_r_ext*FS_shell_gamma_int*m_ex*pow(c_cm,2) );
			// FS_e_diss.push_back( gamma_r_ext*(FS_shell_gamma_int-1)*m_ex*pow(c_cm,2) );

			// Update the mass of the forward shock shell
			FS_shell_mass = m_ex + FS_shell_mass;
			// Update the Lorentz factor of the forward and reverse shock shells
			FS_shell_gamma = gamma_r_ext;
	
			// Record the time of the collision
			FS_te.push_back(tmp_te);
			FS_ta.push_back(tmp_te - FS_shell_radius);
			FS_delt.push_back(RS_shell_radius/2./pow(gamma_r_ext,2.));

			// Record the new Lorentz factor
			FS_gamma_r.push_back(gamma_r_ext);
		}	

		// Check if any more collisions will occur
		if(_check_if_sorted())
		{
			ord_lorentz = true;
			// std::cout << "At time t=" << tmp_te << " s, all shells have been launched and Lorentz factors are ordered.\n";
		}
		// Check if all internal shells have been passed through by the reverse shock
		if(active_inds.size()==0)
		{
			all_inactve = true;
		}
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make the source spectrum using the emission data
void SynthGRB::make_source_spectrum(float energ_min, float energ_max, int num_energ_bins, float tmin, float tmax)
{
	// Release previous spectrum if it exists
	if(p_source_spectrum != NULL)
	{
		delete p_source_spectrum;
	}
	// Make a Spectrum object
	p_source_spectrum = new Spectrum(energ_min, energ_max, num_energ_bins);
	(*p_source_spectrum).ZeroSpectrum(); // Reset spectrum
	

	MakeThermalSpec(p_source_spectrum, tmin, tmax);

	MakeSynchSpec(p_source_spectrum, tmin, tmax);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Thermal spectrum member functions*/

// Calls function to calculate the thermal spectrum rate for each energy bin
void SynthGRB::MakeThermalSpec(Spectrum * therm_spectrum, float tmin, float tmax)
{

	// For each thermal emission event, calculate the emitted spectrum
	for(size_t i=0; i < te_therm.size(); ++i)
	{
		// We only want to take the emission that occurs between the specified tmin and tmax. 
		// The emission occurs between (ta+delt), if any of it overlaps with tmin and tmax, calculate its contribution.
		if ( ta_therm.at(i) <= tmax and (ta_therm.at(i)+delt_therm.at(i)) >= tmin )
		{
			// Call function to calculate thermal spectrum rate
			CalcThermalContribution(therm_spectrum, T_phot.at(i),L_phot.at(i));
		}
		// Else, if this emission occurred outside of the time interval, don't add it to the spectrum.
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Calculate the thermal spectrum from the given temperature and flux of the emission
void SynthGRB::CalcThermalContribution(Spectrum * therm_spectrum, float temp, double flux)
{
	// Calculate the normalization
	double norm=0.; // Set normalization to zero
	int norm_num_bin=30.*log10(E_bol_max/E_bol_min); // Number of energy bins to use to calculate normalization
	std::vector<float> norm_energy_axis(norm_num_bin+1); // Make vector to store energy axis (initialized to zero's) 
	_make_en_axis(norm_energy_axis,E_bol_min,E_bol_max,norm_num_bin+1); // Fill in energy values

	float en_curr=0.; // Current energy to evaluate the addition to the normalization
	// For each energy bin along the normalization energy axis, calculate the addition to the normalization and add it. 
	for(int i=0;i < norm_num_bin; ++i)
	{
		en_curr = (norm_energy_axis.at(i+1) + norm_energy_axis.at(i) ) / 2.; // Set the current energy
		// if(en_curr < 2.*pow(10.,3.))
		// {
			// Calculate the contribution according to a Left-Riemann-Sum
			// norm += (norm_energy_axis.at(i+1)-norm_energy_axis.at(i)) * en_curr * ThermalSpec(en_curr,temp);
		// }
		// else: don't contribute to the normalization
		norm += (norm_energy_axis.at(i+1)-norm_energy_axis.at(i)) * en_curr * ThermalSpec(en_curr,temp);

	}

	// Calculate the rate of the thermal spectrum
	double tmp_val=0.; // Temporary value to store spectrum contributions
	// For each energy bin along the energy axis
	for(int i=0; i < (*therm_spectrum).num_energ_bins; ++i)
	{
		en_curr = (*therm_spectrum).energ_mid.at(i); // Set the current energy
		
		// The current temp and energy bin define the count rate, the normalization found above is applied.
		// This (*p_model_params).must still be (*p_model_params).multiplied by the flux of the source.
		// Also, if the current energy bin is > 2 MeV, then thermal radiation can be ignored.
		// if(en_curr < 2.*pow(10.,3.))
		// {
		// 	tmp_val = flux * ThermalSpec(en_curr,temp) / norm;
		// }
		// else{tmp_val=0;}
		tmp_val = flux * ThermalSpec(en_curr,temp) / norm;

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
	return pow(energy/(kb_kev*temp),1.+alpha)/(exp(energy/(kb_kev*temp)) - 1.);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Synchrotron spectrum member functions*/

// Calls function to calculate the synchrotron spectrum rate for each energy bin
void SynthGRB::MakeSynchSpec(Spectrum * synch_spectrum, float tmin, float tmax)
{
	// For each emission event, calculate the emitted spectrum
	for(size_t i=0; i < te_synch.size(); ++i)
	{
		// We only want to take the emission that occurs between the specified Tmin and Tmax. 
		// The emission occurs between (ta+delt), if any of it overlaps with Tmin and Tmax, calculate its contribution.
		if ( ta_synch.at(i) <= tmax and (ta_synch.at(i)+delt_synch.at(i)) >= tmin )
		{
			// The emission will only be observable if the relativistic velocity is great than the local sound speed v_s/c = 0.1
			// And if the wind is transparent to the radiation
			if ( relvel.at(i) > 0.1 and tau.at(i) < 1)
			{
				// Call function to calculate thermal spectrum rate
				CalcSynchContribution(synch_spectrum, esyn.at(i),e_diss.at(i),delt_synch.at(i));
			}	
		}
		// Else, if this emission occurred outside of the time interval, don't add it to the spectrum. 
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Calculate the synchrotron spectrum from the synchrotron energy and flux of the emission
void SynthGRB::CalcSynchContribution(Spectrum * synch_spectrum, double esyn, double e_diss, double delt)
{
	// Calculate the normalization
	double norm = 0.; // Set normalization to zero
	int norm_num_bin= 30.*log10(E_bol_max/E_bol_min); // Number of energy bins to use to calculate normalization
	std::vector<float> norm_energy_axis(norm_num_bin+1); // Make vector to store energy axis (initialized to zero's) 
	_make_en_axis(norm_energy_axis,E_bol_min,E_bol_max,norm_num_bin+1); // Make energy axis

	float en_curr=0.; // Current energy to evaluate the addition to the normalization
	// For each energy bin along the normalization energy axis, calculate the addition to the normalization and add it. 
	for(int i=0;i < norm_num_bin; ++i)
	{
		en_curr = (norm_energy_axis.at(i+1) + norm_energy_axis.at(i) ) / 2.; // Set the current energy
		// Calculate the contribution according to a Left-Riemann-Sum
		norm += (norm_energy_axis.at(i+1)-en_curr) * en_curr * SynchSpec(en_curr, esyn);
	}

	// Calculate the rate of the thermal spectrum
	double tmp_val=0; // Temporary value to store spectrum contributions
	// For each energy bin along the energy axis
	for(int i=0; i < (*synch_spectrum).num_energ_bins; ++i)
	{
		en_curr = (*synch_spectrum).energ_mid.at(i); // Set the current energy

		// The current temp and energy bin define the count rate, the normalization found above is applied.
		// This (*p_model_params).must still be (*p_model_params).multiplied by the energy dissipated during the emission event.
		// The energy dissipated can be turned into Flux by dividing the energy dissipated by the observed emission duration (delt).
		tmp_val = e_diss * SynchSpec(en_curr, esyn) / norm;
		
		// Add the contribution to the total spectrum according to a Center-Riemann-Sum
		(*synch_spectrum).spectrum_sum += ((*synch_spectrum).energ_hi.at(i) - (*synch_spectrum).energ_lo.at(i)) * tmp_val;

		// Check if the spectrum rate per energy bin is requested and store it if so. 
		(*synch_spectrum).spectrum_rate.at(i) += tmp_val;
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Synchrotron spectrum function form,
// Defaults alpha = -1, beta = -2.5
double SynthGRB::SynchSpec(float energy, double esyn, float alpha, float beta)
{
	double param_list[3] = {esyn, alpha, beta};
	// return (1./esyn) * Band(energy, param_list );
	return (1./esyn) * BPL(energy,  param_list);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make the source light curve using the emission data
void SynthGRB::make_source_light_curve(float energ_min, float energ_max, float Tstart, float Tend, float dt)
{
	// Release previous spectrum if it exists
	if(p_source_light_curve != NULL)
	{
		delete p_source_light_curve;
	} 
	// Make a LightCurve object
	p_source_light_curve = new LightCurve(Tstart , Tend, dt);
	(*p_source_light_curve).ZeroLightCurve(); // Reset spectrum

	// Define the number of energy bins
	float tmp_num_energ_bins = log10(energ_max/energ_min)*20;

	// For each time bin, calculate the photon rate.
	for(int i=0; i < (*p_source_light_curve).num_time_bins-1; ++i)
	{
		// Find the spectrum sum for each emission event which occurs between (light_curve_time[i], light_curve_time[i+1]) 
		make_source_spectrum(energ_min, energ_max, tmp_num_energ_bins, (*p_source_light_curve).lc_time.at(i), (*p_source_light_curve).lc_time.at(i+1) );
		
		// Ensure that the light curve rate is set to zero before adding values to it.
		(*p_source_light_curve).lc_rate.at(i) = (*p_source_spectrum).spectrum_sum;

		// Convert units from erg / s to keV / s
		// (*p_source_light_curve).lc_rate.at(i) *= erg_to_kev;
		// Normalize by the time bin size 
		// (*p_source_light_curve).lc_rate.at(i) /= ((*p_source_light_curve).lc_time.at(1)-(*p_source_light_curve).lc_time.at(0));
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
	thermal_params_file.open(dir_path_name+"synthGRB_jet_params_therm.txt"); // Open text file with this name
	size_t i=0;
	// For each time bin, write the time and count rate to the file.
	while ( i < te_therm.size())
	{
		thermal_params_file << te_therm.at(i);
		thermal_params_file << " ";
		thermal_params_file << ta_therm.at(i);
		thermal_params_file << " ";
		thermal_params_file << delt_therm.at(i);
		thermal_params_file << " ";
		thermal_params_file << T_phot.at(i);
		thermal_params_file << " ";
		thermal_params_file << L_phot.at(i);
		thermal_params_file << " ";
		thermal_params_file << r_phot.at(i);
		thermal_params_file << "\n";		
		++i;
	}
	thermal_params_file.close(); // Close file

	// Write out synchrotron params
	ofstream synch_params_file; // Construct file 
	synch_params_file.open(dir_path_name+"synthGRB_jet_params_synch.txt"); // Open text file with this name
	i=0;
	// For each time bin, write the time and count rate to the file.
	while ( i < te_synch.size())
	{
		synch_params_file << te_synch.at(i);
		synch_params_file << " ";
		synch_params_file << ta_synch.at(i);
		synch_params_file << " ";
		synch_params_file << delt_synch.at(i);
		synch_params_file << " ";
		synch_params_file << asyn.at(i);
		synch_params_file << " ";
		synch_params_file << beq.at(i);
		synch_params_file << " ";
		synch_params_file << gamma_e.at(i);
		synch_params_file << " ";
		synch_params_file << esyn.at(i);
		synch_params_file << " ";
		synch_params_file << gamma_r.at(i);
		synch_params_file << " ";
		synch_params_file << e_diss.at(i);
		synch_params_file << " ";
		synch_params_file << tau.at(i);
		synch_params_file << " ";
		synch_params_file << relvel.at(i);
		synch_params_file << "\n";	
		++i;
	}
	synch_params_file.close(); // Close file

	// Write out external shock params
	ofstream fs_params_file; // Construct file 
	fs_params_file.open(dir_path_name+"synthGRB_jet_params_fs.txt"); // Open text file with this name
	i=0;
	// For each time bin, write the time and count rate to the file.
	while ( i < FS_te.size())
	{
		fs_params_file << FS_te.at(i);
		fs_params_file << " ";
		fs_params_file << FS_ta.at(i);
		fs_params_file << " ";
		fs_params_file << FS_delt.at(i);
		fs_params_file << " ";
		fs_params_file << FS_e_diss.at(i);
		fs_params_file << "\n";	
		++i;
	}
	fs_params_file.close(); // Close file

	ofstream rs_params_file; // Construct file 
	rs_params_file.open(dir_path_name+"synthGRB_jet_params_rs.txt"); // Open text file with this name
	i=0;
	// For each time bin, write the time and count rate to the file.
	while ( i < RS_te.size())
	{
		rs_params_file << RS_te.at(i);
		rs_params_file << " ";
		rs_params_file << RS_ta.at(i);
		rs_params_file << " ";
		rs_params_file << RS_delt.at(i);
		rs_params_file << " ";
		rs_params_file << RS_e_diss.at(i);
		rs_params_file << "\n";	
		++i;
	}
	rs_params_file.close(); // Close file
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Save spectrum to text file
void SynthGRB::WriteSpectrumToTXT(std::string out_file_name)
{
	// Write the light curve to a text file
	ofstream spec_file; // Construct file 
	spec_file.open(out_file_name); // Open text file with this name
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

