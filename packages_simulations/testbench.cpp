#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <tuple>
#include <cfitsio/fitsio.h>
#include <list>

using namespace std;

// Import Custom Libraries
// #include "SynthGRB.hpp"
// #include "SynthGRBLibrary.hpp"
// #include "ObsGRB.hpp"
// #include "Response.hpp"
// #include "LightCurve.hpp"
// #include "Spectrum.hpp"
#include "cosmology.hpp"
#include "utilfuncs.hpp"
// #include "DataAnalysis.hpp"
// #include "TTEs.hpp"
#include "ModelParams.hpp"
#include "ShellDist.hpp"
// #include "GammaFunction.hpp"

using namespace std;


// Check if shells are ordered
int _check_if_sorted(ShellDist *p_jet_shells)
{
	// Make list of shells which are currently active:
	std::vector<int> active_inds;
	for(int i=0; i<5000; ++i)
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

// Find the resulting Lorentz factor of two collided shells from the Lorentz factor and Mass of shell 1 and 2
double _calc_shell_coll_gamma(float s1g, float s2g, float s1m, float s2m)
{
	// The approximate Lorentz factor of the combined layer after complete redistribution of momentum and energy
	return sqrt( s1g*s2g * (s1m*s1g + s2m*s2g) / (s1m*s2g + s2m*s1g)  );
}

int main(int argc, char const *argv[])
{

	float tw = 10; 
	float dte = 0.002;
	float E_dot_iso = 1e+52;
	float theta = 0.07;
	float r_open = 3000000.0;
	float eps_th = 0.02;
	float sigma = 0.1;
	float eps_e_int = 0.3;
	float eps_b_int = 0.07;
	float zeta_int = 0.001;
	float p_int = 2.2;
	float eps_e_ext = 0.1;
	float eps_b_ext = 0.0001;
	float zeta_ext = 1;
	float p_ext = 2.2;
	float k_med = 0;
	float rho_not = 1.672e-24;
	std::string LorentzDist = "step";
	std::string ShellDistParamsFile = "Default";

	ModelParams * p_model_params; // Stores all model parameters	
	p_model_params = new ModelParams(tw, dte, 
	E_dot_iso, theta, r_open, eps_th, sigma,
	eps_e_int, eps_b_int, zeta_int, p_int,
	eps_e_ext, eps_b_ext, zeta_ext, p_ext,
	k_med, rho_not,
	LorentzDist, ShellDistParamsFile);

	ShellDist * p_jet_shells; // Pointer to the jet shells 
	p_jet_shells = NULL;
	p_jet_shells = new ShellDist(5000);
	(*p_jet_shells).step(0.002, 100, 400, 0.5, false);


	// Average Lorentz factor of all jet shells  
	float gamma_bar; // Average Lorentz factor of all jet shells 
	float m_bar; // Average mass of each shell
	// double m_tot; // Average mass of each shell
	gamma_bar = 0.0;
	for(int i=0; i<(*p_model_params).numshells; ++i) 
	{
	     gamma_bar += (*p_jet_shells).shell_gamma.at(i);
	}
	gamma_bar /= (*p_model_params).numshells;
	// Average mass of each shell
	m_bar = (*p_model_params).E_dot_iso*(*p_model_params).tw/gamma_bar/pow(c_cm,2.)/(*p_model_params).numshells;
	// m_tot = (*p_model_params).E_dot_iso*(*p_model_params).tw/gamma_bar/pow(c_cm,2.);


	for(int i=0; i<(*p_model_params).numshells; ++i) 
	{
	    (*p_jet_shells).shell_radius.at(i) += 0.3e2;
	}

	bool ord_lorentz = false; // Indicates whether the internal shells are ordered by their Lorentz factor
	double min_mass = 0.; // g/M_tot, mass of the less massive shell involved in the internal shock
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

	double tmp_te = 0.; // sec, Temporary value of the emission time
	double tmp_t_IS = 0.; // sec, Time until collision of the upstream and downstream shells
	// double t_IS_lowest = 0.; // sec, Minimum time until next shell collision	
	double t_IS_lowest = std::numeric_limits<double>::max();

	double rad_coll = 0.; // light seconds, Radius of the collision
	// double tmp_r_before = 0.; // light seconds, Temporary radius before the time step

	double tmp_gamma_r = 1.; // Temporary value of the combined bulk Lorentz factor two colliding shells  
	double tmp_e_diss = 0.; // erg Temporary value of the energy dissipated during the collision
	double tmp_eps_star = 0.; // erg / g, Average proton factor from the collision of two shells

	while(ord_lorentz == false)
	{
	t_IS_lowest = std::numeric_limits<double>::max();

		// Make list of shells which are currently active:
		std::vector<int> active_inds;
		for(int i=0; i<(*p_model_params).numshells; ++i)
		{
			if((*p_jet_shells).shell_status.at(i) == 1)
			{
				active_inds.push_back(i);
			}
		}

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

		std::cout << "tmp_ds_r = " << (*p_jet_shells).shell_radius.at(ind_s_ds) << std::endl;
		std::cout << "tmp_ds_b = " << beta( (*p_jet_shells).shell_gamma.at(ind_s_ds) ) << std::endl;
		std::cout << "tmp_us_r = " << (*p_jet_shells).shell_radius.at(ind_s_us) << std::endl;
		std::cout << "tmp_us_b = " << beta( (*p_jet_shells).shell_gamma.at(ind_s_us) ) << std::endl;

		// Move to the emission time of the shell collision (in the rest frame of the central engine)
		tmp_te += t_IS_lowest;
		// Move shells forward
		for(int i=0; i<5000; ++i)
		{
			if((*p_jet_shells).shell_status.at(i) == 1)
			{
				(*p_jet_shells).shell_radius.at(i) += beta((*p_jet_shells).shell_gamma.at(i)) * t_IS_lowest;
			}
		}

		// De-active the downstream shell now that it has been merged
		(*p_jet_shells).shell_status.at(ind_s_ds) = 0;

		// For ease of reading, make variables for the masses and Lorentz factors of the two colliding shells
		// Down stream shell
		shell_ds_m = (*p_jet_shells).shell_mass.at(ind_s_ds);
		shell_ds_g = (*p_jet_shells).shell_gamma.at(ind_s_ds);
		// Upstream shell
		shell_us_m = (*p_jet_shells).shell_mass.at(ind_s_us);
		shell_us_g = (*p_jet_shells).shell_gamma.at(ind_s_us);

		// Record the radius of the collision
		rad_coll = (*p_jet_shells).shell_radius.at(ind_s_us);


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


		std::cout << "t_IS_lowest = " << t_IS_lowest << std::endl;
		std::cout << "tmp_te = " << tmp_te << std::endl;
		std::cout << "rad_coll = " << rad_coll << std::endl;
		std::cout << "min m = " << min_mass << std::endl;
		std::cout << "US Gam = " << shell_us_g << std::endl;
		std::cout << "US Mass = " << shell_us_m << std::endl;
		std::cout << "DS Gam = " << shell_ds_g << std::endl;
		std::cout << "DS Mass = " << shell_ds_m << std::endl;
		std::cout << "Gam_r = " << tmp_gamma_r << std::endl;
		std::cout << " " << std::endl;


		// Energy dissipated in shock
		tmp_e_diss = m_bar * pow(c_cm,2.) * (shell_ds_m*shell_ds_g + shell_us_m*shell_us_g - (shell_us_m + shell_ds_m)*tmp_gamma_r); // erg
		
		// gamma_int = 0.5 *( sqrt(shell_ds_g/shell_us_g) + sqrt(shell_us_g/shell_ds_g) ); // Lorentz factor for internal motion in shocked material
		// tmp_eps_star = (gamma_int - 1) * pow(c_cm,2.); // erg / g, Average proton factor from the collision of two shells
		tmp_eps_star = tmp_e_diss/(2.*m_bar*min_mass*tmp_gamma_r);

		// Set the Lorentz factor of the merged shells
		(*p_jet_shells).shell_gamma.at(ind_s_us) = tmp_gamma_r;
		// Set the mass of the merged shell (i.e., add the mass of the two shells)
		(*p_jet_shells).shell_mass.at(ind_s_us) += (*p_jet_shells).shell_mass.at(ind_s_ds);


		// Check if any more collisions will occur
		if(_check_if_sorted(p_jet_shells))
		{
			ord_lorentz = true;
		}	


	}

	return 0;

	// for(int i=0; i<5000; ++i)
	// {
	// 	std::cout << (*p_jet_shells).shell_gamma.at(i) <<std::endl;

	// }


}