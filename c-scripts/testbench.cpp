#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <tuple>
#include "fitsio.h"

using namespace std;

// Import Custom Libraries
#include "SynthGRB.hpp"
#include "SynthGRBLibrary.hpp"
#include "ObsGRB.hpp"
#include "Response.hpp"
#include "LightCurve.hpp"
#include "Spectrum.hpp"
#include "cosmology.hpp"
#include "utilfuncs.hpp"
#include "DataAnalysis.hpp"
#include "TTEs.hpp"
#include "ModelParams.hpp"
#include "ShellDist.hpp"

int main(int argc, char const *argv[])
{

	return 0;
}

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

	// Declare and initialize flags which indicate whether the jet propagation should continue
	bool ord_lorentz = false; // Indicates whether the internal shells are ordered by their Lorentz factor
	bool fs_relativ = true; // Indicates whether the forward shock is still relativistic.
	bool fs_active = false; // Indicates whether the FS is active or not
	bool rs_active = false; // Indicates whether the RS is active or not

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
	double shell_ds_g = 0.; // Downstream shell Lorentz factor 
	double shell_us_m = 0.; // g/M_tot, Upstream shell mass 
	double shell_us_g = 0.; // Upstream shell Lorentz factor 

	double tmp_t_IS = 0.; // sec, Time until collision of the upstream and downstream shells
	double t_IS_lowest = 0.; // sec, Minimum time until next shell collision	

	double tmp_tau = 0.; // Temporary value of the optical depth at the site of the collision

	double t_syn = 0.; // sec, Synchrotron time-scale
	double gamma_comb = 0.; // Final Lorentz factor of the combined shell (after complete redistribution of momenta)
	double rad_coll = 0.; // light seconds, Radius of the collision


	// Declare and initialize variables for 
	//External Shock processes

	double t_FS_sweep = 0.; // sec, Time until a mass m_ex is swept up by the FS   
	double t_RS_coll = 0.; // sec, Time until a RS collision

	double tmp_ej_m = 0.; // g/M, Temporary mass of the ejecta shell 
	float tmp_ej_g = 0.; // Temporary Lorentz factor of the ejecta shell 
	double tmp_fs_m = 0.; // g/M, Temporary mass of the forward shock shell
	float tmp_fs_g = 0.; // Temporary Lorentz factor of the forward shock shell
	double tmp_rs_m = 0.; // g/M, Temporary mass of the reverse shock shell
	float tmp_rs_g = 0.; // Temporary Lorentz factor of the reverse shock shell

	double tmp_r_before = 0.; // light seconds, Temporary radius before the time step
	double tmp_m_ex = 0; // g/M_tot, Mass limit for a forward shock to occur 

	int fs_shell_index = NULL; // To keep track of which shell is representing the FS
	int rs_shell_index = NULL;// To keep track of which shell is representing the RS

	double fs_gamma_int = 0.; // Internal Lorentz factor for the forward shock
	// double rs_gamma_int = 0.; // Internal Lorentz factor for the forward shock

	float R_ext0 = 6.9*pow(10.,10.)/c_cm; // cm/speed of light = light seconds, radius of the sun. (Approximate radius of a wolf-rayet star?)  

	// float k = 2.;// Wind density profile for external medium
	// float Astar = 1; // Normalization of the stellar wind parameter 
	// float rho_not = Astar*6.6*pow(10.,12.); // g cm^-3, Normalization of external the density

	float k = 0.;// Constant density profile for external medium
	float n0 = 1; // cm^-3, particle number density at R_ext0
	double rho_not = n0 * mp; // g cm^-3, Normalization of the external density

	float q = 0.01; // Ratio which dictates how much external mass needs to be swept up in order for an internal show to occur
	double m_ex = q / gamma_bar; // A value to parameterize when a FS occurs


	// Shared variables between the processes
	double tmp_te = 0.; // sec, Temporary value of the emission time
	double tmp_asyn = 0.; // Temporary value of the fraction of energy going to synchrotron electrons
	double tmp_beq = 0.; // Temporary value of the magnetic field
	double tmp_gamma_e = 0.; // Temporary value of the typical Lorentz factor of accelerated electrons 
	double tmp_esyn_kev = 0.; // keV, Temporary value of the typical energy of accelerated electrons
	double tmp_esyn_erg = 0.; // erg, 
	double tmp_gamma_r = 0.; // Temporary value of the combined bulk Lorentz factor two colliding shells  
	double tmp_e_diss = 0.; // erg/s Temporary value of the energy dissipated during the collision

	double gamma_int = 0.; // Lorentz factor for internal motion in shocked material
	double eps = 0.; // erg, Average proton factor from the collision of two shells
	double rho = 0.; // g/cm^3, Comoving proton density
	double quadr_const = 0.; // Constant in the quadratic equation
	double w = 0.; // Indicates whether we are in the Klein-Nishina regime or Thompson (i.e., w>>1 or w<<1, respectively)
	double Q_IC = 0.; // Inverse Compton parameter
	double alpha_ic = 0.; // Fraction of energy that goes into Inverse Compton electrons

	while( (ord_lorentz==false) | (all_inactve==false) | (fs_relativ == true) )
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

		/* Calculate time until next FS */
		// A forward shock will occur every time a mass of m_ex is swept up, which occurs in a time t_FS_sweep
		// m_swept = q*M_tot/Gamma_FS; // g/M, m_ex is a fraction q of the total external shock mass

		// Mass can only be swept up when the outermost shell has passed the start of the external medium
		if( (*p_jet_shells).shell_radius.at(active_inds.at(0)) > R_ext0 )
		{
			// Below we calculate the mass swept up between the shell position before and after the time step (i.e., R_2 - R_1 = beta * t)
			// Assume a wind density profile for the medium:
			if (k == 2)
			{
				m_swept += 4 * M_PI * rho_not * beta( (*p_jet_shells).shell_gamma.at(active_inds.at(0)) ) * t_coll_lowest; // g/M
			}
			// Assume a constant density profile for the medium:
			else if (k == 0)
			{
				r_before = (*p_jet_shells).shell_radius.at(active_inds.at(0)) - beta((*p_jet_shells).shell_gamma.at(active_inds.at(0))) * t_coll_lowest;
				m_swept += 4 * M_PI * rho_not * ( pow(r_curr ,3.) - pow(r_before ,3.) ); // g/M
			} 

			// Calculate the mass still needed to be swept up until a FS occurs
			m_remaining = m_ex - m_swept; 
		}

		// If the forward shock has begun, calculate the time until the next FS event		
		if(fs_active == true)
		{
			// Assume a wind density profile for the medium:
			if (k == 2)
			{
				FS_t_sweep = m_remaining*m_bar / (4. * M_PI * rho_not) / (beta((*p_jet_shells).shell_gamma.at(fs_shell_index))*c_cm); // sec
			}
			// Assume a constant density profile for the medium:
			else if (k == 0)
			{
				FS_t_sweep = ( pow( 3.*m_remaining*m_bar / (4. * M_PI * rho_not) + pow((*p_jet_shells).shell_gamma.at(fs_shell_index)*c_cm,3.),1./3.) - (*p_jet_shells).shell_gamma.at(fs_shell_index)*c_cm)/ (beta((*p_jet_shells).shell_gamma.at(fs_shell_index))*c_cm); // sec
			}			
		}

		/* Calculate time until next RS */
		// If the reverse shock has begun, calculate the time until the next RS event
		if(rs_active == true)
		{
			// Calculate when the outer most shell will collide with the external shock discontinuity, e.g., when a reverse shock will happen
			t_RS_coll = ((*p_jet_shells).shell_radius.at(rs_shell_index) - (*p_jet_shells).shell_radius.at(active_inds.at(0))) / (beta((*p_jet_shells).shell_gamma.at(active_inds.at(0))) - beta((*p_jet_shells).shell_gamma.at(rs_shell_index)));
		}

		/* Internal Shock Dynamics and Emission*/
		if((t_coll_lowest < t_FS_sweep) & (t_coll_lowest < t_RS_coll))
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

			// Sweep up some mass into the FS
			//  If the FS hasn't begun yet, use the forward use the outer most shell to sweep up mass
			if(fs_active == false)
			{
				// Assume a wind density profile for the medium:
				if (k == 2)
				{
					m_swept += 4 * M_PI * rho_not * beta( (*p_jet_shells).shell_gamma.at(active_inds.at(0)) ) * t_coll_lowest; // g/M
				}
				// Assume a constant density profile for the medium:
				else if (k == 0)
				{
					tmp_r_before = (*p_jet_shells).shell_radius.at(active_inds.at(0)) - beta((*p_jet_shells).shell_gamma.at(active_inds.at(0))) * t_coll_lowest;
					m_swept += 4 * M_PI * rho_not * ( pow((*p_jet_shells).shell_radius.at(active_inds.at(0)),3.) - pow(tmp_r_before ,3.) ); // g/M
				} 

				// If the mass swept up is greater than m_ex, a forward shock begins (and therefor a reverse shock)
				if(m_swept >= m_ex)
				{				
					fs_active = true;
					rs_active = true;

					// Turn the outermost shell into the FS
					fs_shell_index = active_inds.at(0);
					(*p_jet_shells).shell_status.at(fs_shell_index) = 0; // Deactivate the shell so it does not participate in any further IS calculations
					// Turn the second outermost shell into the RS
					rs_shell_index = active_inds.at(1);
					(*p_jet_shells).shell_status.at(rs_shell_index) = 0; // Deactivate the shell so it does not participate in any further IS calculations
				}
			}
			// Else, the FS is already active
			else
			{
				// Move the FS forward
				(*p_jet_shells).shell_radius.at(fs_shell_index) += beta((*p_jet_shells).shell_gamma.at(fs_shell_index))*t_IS_lowest;

				// Assume a wind density profile for the medium:
				if (k == 2)
				{
					m_swept += 4 * M_PI * rho_not * beta( (*p_jet_shells).shell_gamma.at(fs_shell_index) ) * t_IS_lowest; // g/M
				}
				// Assume a constant density profile for the medium:
				else if (k == 0)
				{
					tmp_r_before = (*p_jet_shells).shell_radius.at(fs_shell_index) - beta((*p_jet_shells).shell_gamma.at(fs_shell_index)) * t_IS_lowest;
					m_swept += 4 * M_PI * rho_not * ( pow((*p_jet_shells).shell_radius.at(fs_shell_index),3.) - pow(tmp_r_before ,3.) ); // g/M
				} 
				// Calculate the mass still needed to be swept up until a FS occurs
				m_remaining = m_ex - m_swept; 
			}
			if(rs_active == true)
			{
				(*p_jet_shells).shell_radius.at(rs_shell_index) += beta((*p_jet_shells).shell_gamma.at(rs_shell_index))*t_IS_lowest;

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
			eps = (gamma_int - 1) * pow(c_cm,2.); // erg / g, Average proton factor from the collision of two shells

			// Characteristic Lorentz factor of accelerated electrons
			// tmp_gamma_e = ((*p_model_params).p-2) * (*p_model_params).eps_e * mp * eps / ((*p_model_params).p-1) / zeta / me / pow(c_cm,2); // Minimum Lorentz factor for a population of electrons distributed as a power law with index p
			tmp_gamma_e = 1e4; // Constant

			rho = E_dot_kin / (4. * M_PI * pow(rad_coll*gamma_bar,2.) * pow(c_cm,5.)); // g/cm^3, Comoving proton mass density
			tmp_beq = sqrt( 8.* M_PI * (*p_model_params).eps_b * rho * eps); // (erg/cm^3)^1/2, Magnetic field density, assuming B = Beq  
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
		
		/* Forward Shock Dynamics and Emission*/
		else if((fs_active == true ) & (t_FS_sweep < t_RS_coll))
		{

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
			(*p_jet_shells).shell_radius.at(fs_shell_index) += beta((*p_jet_shells).shell_gamma.at(fs_shell_index))*t_FS_sweep;
			(*p_jet_shells).shell_radius.at(rs_shell_index) += beta((*p_jet_shells).shell_gamma.at(rs_shell_index))*t_FS_sweep;

			// Sweep up some mass into the FS
			// Assume a wind density profile for the medium:
			if (k == 2)
			{
				m_swept += 4 * M_PI * rho_not * beta( (*p_jet_shells).shell_gamma.at(fs_shell_index) ) * t_FS_sweep; // g/M
			}
			// Assume a constant density profile for the medium:
			else if (k == 0)
			{
				tmp_r_before = (*p_jet_shells).shell_radius.at(fs_shell_index) - beta((*p_jet_shells).shell_gamma.at(fs_shell_index)) * t_FS_sweep;
				m_swept += 4 * M_PI * rho_not * ( pow((*p_jet_shells).shell_radius.at(fs_shell_index),3.) - pow(tmp_r_before ,3.) ); // g/M
			}

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
			tmp_gamma_r = pow( ((tmp_rs_m + tmp_fs_m*fs_gamma_int)*pow(tmp_fs_g,2.) + m_swept*tmp_fs_g ) / ( (tmp_rs_m + tmp_fs_m*fs_gamma_int) + 2*m_swept*tmp_fs_g ) ,0.5);

			// The new internal Lorentz factor in the shocked forward medium due to forward shock, if:
			// adiabatic (slow cooling) so FS_gamma_int > 1 and approx = gamma_r  :
			fs_gamma_int = ((tmp_rs_m + tmp_fs_m*fs_gamma_int)*tmp_fs_g + m_swept - tmp_rs_m*tmp_gamma_r)/(tmp_fs_m+m_swept)/tmp_gamma_r;

			eps =  (fs_gamma_int - 1) * pow(c_cm,2.); // erg/g, Average energy dissipated for each proton in the collision of two shells

			// The typical Lorentz factor of the electron distribution (Gamma_min)
			// If we assume zeta is approximately = 1 in the FS
			tmp_gamma_e = ((*p_model_params).p-2.) * (*p_model_params).eps_e * mp * eps / ((*p_model_params).p-1.) / me / pow(c_cm,2.); // Minimum Lorentz factor for a population of electrons distributed as a power law with index p
			
			// The magnetic field in the forward shock
			tmp_rho = rho_not/pow(tmp_fs_r*c_cm,k); // g/cm^3, comoving proton density 
			tmp_beq = sqrt(8 * M_PI * (*p_model_params).eps_b * tmp_rho * eps);

			// The synchrotron energy of the accelerated electrons 
			tmp_esyn_kev = 50.*(tmp_gamma_r/300.)*(tmp_beq/1e3)*pow((tmp_gamma_e/100.),2.) / 1000.; // keV, Synchrotron energy in the rest frame
			tmp_esyn_erg = kev_to_erg * tmp_esyn_kev; // erg, Synchrotron energy in the rest frame

			w = tmp_gamma_e*(tmp_esyn_erg/tmp_gamma_r)/(me * pow(c_cm,2.) ); // Indicates whether we are in the Klein-Nishina regime or Thompson (i.e., w>>1 or w<<1, respectively)
			quadr_const = -(3./2.)*(0.2/M_PI/pow(c_cm*tmp_fs_r,2.)) * tmp_gamma_e * 100. * (E_dot_kin/pow(gamma_bar*c_cm,2.)) * pow(tmp_beq/1000.,-2.); // Constant in the quadratic equation
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

			// The energy dissipated by forward shock
			// tmp_e_diss = pow(c_cm,2) * m_bar * ( (tmp_rs_m + tmp_fs_m*fs_gamma_int)*tmp_fs_g + m_swept - (tmp_rs_m + m_swept + tmp_fs_m*fs_gamma_int )*tmp_gamma_r) * (*p_model_params).eps_e * tmp_asyn ; // erg
			tmp_e_diss = pow(c_cm,2.) * m_bar * ( (tmp_rs_m + tmp_fs_m*fs_gamma_int)*tmp_fs_g + m_swept - (tmp_rs_m + m_swept + tmp_fs_m*fs_gamma_int )*tmp_gamma_r) * (*p_model_params).eps_e * tmp_asyn / (tmp_fs_r/2./pow(tmp_gamma_r,2.)); // erg / s;  

			// Update the Lorentz factor of the forward and reverse shock shells
			(*p_jet_shells).shell_gamma.at(fs_shell_index) = tmp_gamma_r;
			(*p_jet_shells).shell_gamma.at(rs_shell_index) = tmp_gamma_r;

			// Update the mass of the forward shock shell
			(*p_jet_shells).shell_mass.at(fs_shell_index) += m_swept;

			// Reset the swept up mass
			m_swept = 0.;

			// Continuous emission:
			FS_te.push_back(tmp_te); // sec, Time of emission 
			FS_ta.push_back(tmp_te - RS_shell_radius); // sec, Arrival time of emission at observer
			FS_delt.push_back(tmp_fs_r/2./pow(tmp_gamma_r,2.)); // sec, Duration of the emission
			FS_gamma_r.push_back(tmp_gamma_r); // New Lorentz factor of the FS
			FS_esyn.push_back(tmp_esyn_kev); // erg, The typical synchrotron energy of the accelerated electrons
			FS_e_diss.push_back(tmp_e_diss); // erg/s, Dissipated energy of the FS	

		}

		/* Reverse Shock Dynamics and Emission*/
		else if((rs_active == true ))
		{
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
			(*p_jet_shells).shell_radius.at(fs_shell_index) += beta((*p_jet_shells).shell_gamma.at(fs_shell_index))*t_RS_coll;
			(*p_jet_shells).shell_radius.at(rs_shell_index) += beta((*p_jet_shells).shell_gamma.at(rs_shell_index))*t_RS_coll;

			// Sweep up some mass into the FS
			// Assume a wind density profile for the medium:
			if (k == 2)
			{
				m_swept += 4 * M_PI * rho_not * beta( (*p_jet_shells).shell_gamma.at(fs_shell_index) ) * t_RS_coll; // g/M
			}
			// Assume a constant density profile for the medium:
			else if (k == 0)
			{
				tmp_r_before = (*p_jet_shells).shell_radius.at(fs_shell_index) - beta((*p_jet_shells).shell_gamma.at(fs_shell_index)) * t_RS_coll;
				m_swept += 4 * M_PI * rho_not * ( pow((*p_jet_shells).shell_radius.at(fs_shell_index),3.) - pow(tmp_r_before ,3.) ); // g/M
			} 
			// Calculate the mass still needed to be swept up until a FS occurs
			m_remaining = m_ex - m_swept; 


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
			
			gamma_int = 0.5 *( sqrt(tmp_rs_g/tmp_ej_g) + sqrt(tmp_ej_g/tmp_rs_g) ); // Lorentz factor for internal motion in shocked material
			eps =  (gamma_int - 1) * pow(c_cm,2.); // erg / g, Average energy per proton within the collided shells

			// The typical Lorentz factor of the electron distribution (Gamma_min)
			tmp_gamma_e = ((*p_model_params).p-2) * (*p_model_params).eps_e * mp * eps / ((*p_model_params).p-1) / (*p_model_params).zeta / me / pow(c_cm,2.); // Minimum Lorentz factor for a population of electrons distributed as a power law with index p
			
			// Density in the reverse shock medium			
			tmp_rho = E_dot_kin / (4. * M_PI * pow(tmp_rs_r*gamma_bar,2.) * pow(c_cm,5.)); // g/cm^3, Comoving proton mass density
			
			// The magnetic field in the reverse shock
			tmp_beq = sqrt(8. * M_PI * (*p_model_params).eps_b * tmp_rho * eps); // (erg/cm^3)^1/2, Magnetic field density, assuming B = Beq  

			// The typical energy of an accelerated electron
			tmp_esyn_kev = 50.*(tmp_gamma_r/300.)*(tmp_beq/1e3)*pow((tmp_gamma_e/100.),2.) / 1000.; // keV, Synchrotron energy in the source frame
			tmp_esyn_erg = kev_to_erg * tmp_esyn_kev; // erg, Synchrotron energy in the rest frame 

			w = tmp_gamma_e*(tmp_esyn_erg/tmp_gamma_r)/(me * pow(c_cm,2.) ); // Indicates whether we are in the Klein-Nishina regime or Thompson (i.e., w>>1 or w<<1, respectively)
			quadr_const = -(3./2.)*(0.2/M_PI/pow(c_cm*tmp_rs_r,2.)) * tmp_gamma_e * 100. * (E_dot_kin/pow(gamma_bar*c_cm,2.)) * pow(tmp_beq/1000.,-2.); // Constant in the quadratic equation
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

			// The energy dissipated by reverse shock
			// tmp_e_diss = pow(c_cm,2) * m_bar * ( (tmp_rs_m + tmp_fs_m*fs_gamma_int)*tmp_fs_m + tmp_ej_m*tmp_ej_g - (tmp_rs_m + tmp_ej_m+ tmp_fs_m*fs_gamma_int )*tmp_gamma_r) * (*p_model_params).eps_e * tmp_asyn; // erg
			tmp_e_diss = pow(c_cm,2.) * m_bar * ( (tmp_rs_m + tmp_fs_m*fs_gamma_int)*tmp_rs_g + tmp_ej_m*tmp_ej_g - (tmp_rs_m + tmp_ej_m+ tmp_fs_m*fs_gamma_int )*tmp_gamma_r) * (*p_model_params).eps_e * tmp_asyn / (tmp_rs_r/2./pow(tmp_gamma_r,2.)); // erg / s;

			// Update the mass of the reverse shock shell
			(*p_jet_shells).shell_mass.at(rs_shell_index) += tmp_ej_m;
			// Update the Lorentz factor of the forward and reverse shock shells
			(*p_jet_shells).shell_gamma.at(fs_shell_index) = tmp_gamma_r;
			(*p_jet_shells).shell_gamma.at(rs_shell_index) = tmp_gamma_r;

			// Record necessary parameters to calculate the emission
			RS_te.push_back(tmp_te); // sec, Time of emission 
			RS_ta.push_back(tmp_te - tmp_rs_r); // sec, Arrival time of emission at observer
			RS_delt.push_back(tmp_rs_r/2./pow(tmp_gamma_r,2.)); // sec, Duration of the emission
			RS_gamma_r.push_back(tmp_gamma_r); // New Lorentz factor of the RS
			RS_esyn.push_back(tmp_esyn_kev); // erg, The typical synchrotron energy of the accelerated electrons
			RS_e_diss.push_back(tmp_e_diss); // erg, Dissipated energy of the RS

			// Deactivate internal shell that collided with RS
			(*p_jet_shells).shell_status.at(active_inds.at(0)) = 0;
		}
		else
		{
			cout << "This shouldn't be possible" << endl;
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
			rs_active = false;
			// std::cout << "At time t=" << tmp_te << " s, all shells have passed through the reverse shock.\n";
		}
		// Check if the FS is still relativistic 
		if((*p_jet_shells).shell_gamma.at(fs_shell_index) < 2)
		{
			fs_relativ = false;
			// std::cout << "At time t=" << tmp_te << " s, the forward shock is no longer relativistic.\n";
		}
		else{ fs_relativ = true;}

	}


//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Perform simulation of jet dynamics using the loaded jet parameters
void SynthGRB::SimulateJetDynamics()
{
		
		// else, RS happens before IS, but we still require some internal shells to be active 
		// else if( (active_inds.size()>0) & (RS_active == 1) )
		else if( (RS_t_coll < FS_t_sweep) & (active_inds.size()>0) )
		{
			
		}
		// No RS or IS remain, but the FS may still propagate into the external medium
		else
		{
			
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
			// std::cout << "At time t=" << tmp_te << " s, all shells have passed through the reverse shock.\n";
		}
		// Check if the FS is still relativistic 
		if(tmp_fs_g < 2)
		{
			fs_relativ = false;
			// std::cout << "At time t=" << tmp_te << " s, the forward shock is no longer relativistic.\n";
		}
		else
		{
			fs_relativ = true;
		}
	}

}