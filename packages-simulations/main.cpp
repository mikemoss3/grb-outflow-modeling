
/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-20

Main function to create and use Spectra, Light Curves, and Response functions

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
#include "SynthGRB.hpp"
#include "ObsGRB.hpp"
#include "Response.hpp"
// #include "LightCurve.hpp"
// #include "Spectrum.hpp"
#include "cosmology.hpp"
#include "utilfuncs.hpp"
#include "DataAnalysis.hpp"
#include "TTEs.hpp"
#include "ModelParams.hpp"

using namespace std;

int main(int argc, char const *argv[])
{

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* 
	Called during interactive plotting manipulation

	argv[1] = indicates if the update changes the plotted time interval or the energy band

	argv[2] = indicates if all the components should be calculated or just the summed total 
	argv[3] = indicates if a log scale should be used to calculate the light curve
	argv[4], argv[5] = input time interval
	argv[6], argv[7] = input energy interval

	*/

	if(argv[1]!=NULL)
	{

		bool logscale = false;
		if( (strcmp(argv[3],"True")==0) | (strcmp(argv[3],"true")==0) | (strcmp(argv[3],"1")==0) )
		{
			logscale = true;
		}
		
		if( strcmp(argv[1], "timechange")==0 )
		{

			float tmin = stof(argv[4]);
			float tmax = stof(argv[5]);

			float energ_min = stof(argv[6]);
			float energ_max = stof(argv[7]);
			float num_energ_bins = 10.*log10(energ_max/energ_min);

			SynthGRB synth_grb = SynthGRB();
			synth_grb.LoadJetParamsFromTXT("input-files/jet-params.txt");
			synth_grb.SimulateJetDynamics();

			// Total spectrum
			synth_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tmin, tmax);
			synth_grb.WriteSpectrumToTXT("data-file-dir/quickplot_spectrum.txt");

			// Make component spectra
			if(strcmp(argv[2], "true")==0)
			{
				// Thermal Component
				synth_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tmin, tmax, "TH");
				synth_grb.WriteSpectrumToTXT("data-file-dir/quickplot_spectrum_TH.txt");

				// Internal Shock Component
				synth_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tmin, tmax, "IS");
				synth_grb.WriteSpectrumToTXT("data-file-dir/quickplot_spectrum_IS.txt");

				// Forward Shock Component
				synth_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tmin, tmax, "FS");
				synth_grb.WriteSpectrumToTXT("data-file-dir/quickplot_spectrum_FS.txt");

				// Reverse Shock Component
				synth_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tmin, tmax, "RS");
				synth_grb.WriteSpectrumToTXT("data-file-dir/quickplot_spectrum_RS.txt");
			}
			
			return 0; 
		}
		if( strcmp(argv[1], "energychange")==0 )
		{

			float tmin = stof(argv[4]);
			float tmax = stof(argv[5]);
			float dt = stof(argv[6]);

			float energ_min = stof(argv[7]);
			float energ_max = stof(argv[8]);

			SynthGRB synth_grb = SynthGRB();
			synth_grb.LoadJetParamsFromTXT("input-files/jet-params.txt");
			synth_grb.SimulateJetDynamics();

			synth_grb.make_source_light_curve(energ_min, energ_max, tmin, tmax, dt, "all", logscale);
			synth_grb.WriteLightCurveToTXT("data-file-dir/quickplot_light_curve.txt");
			

			// Make component spectra and light curves
			if(strcmp(argv[2], "true")==0)
			{
				// Plot the other light curves components 
				synth_grb.make_source_light_curve(energ_min, energ_max, tmin, tmax, dt, "TH", logscale);
				synth_grb.WriteLightCurveToTXT("data-file-dir/quickplot_light_curve_TH.txt");
				synth_grb.make_source_light_curve(energ_min, energ_max, tmin, tmax, dt, "IS", logscale);
				synth_grb.WriteLightCurveToTXT("data-file-dir/quickplot_light_curve_IS.txt");
				synth_grb.make_source_light_curve(energ_min, energ_max, tmin, tmax, dt, "FS", logscale);
				synth_grb.WriteLightCurveToTXT("data-file-dir/quickplot_light_curve_FS.txt");
				synth_grb.make_source_light_curve(energ_min, energ_max, tmin, tmax, dt, "RS", logscale);
				synth_grb.WriteLightCurveToTXT("data-file-dir/quickplot_light_curve_RS.txt");				
			}

			return 0; 
		}
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* Testing SynthGRB default Light Curve and Spectrum making */
	

	float energ_min = 8;
	float energ_max = 5.e4;
	float energ_min_lc = energ_min; 
	float energ_max_lc = energ_max;
	float num_energ_bins = 400;

	// 15 - 5000 keV = Konus-Wind energy band

	// Light curve time interval
	float tmin = 0.;
	float tmax = 10.;
	float dt = 0.05;

	SynthGRB test_grb = SynthGRB();
	test_grb.LoadJetParamsFromTXT("input-files/jet-params.txt");
	
	test_grb.anim_lor_dist = true;
	(*test_grb.p_jet_shells).WriteToTXT("data-file-dir/synthGRB_shell_dist.txt");

	test_grb.SimulateJetDynamics();
	test_grb.write_out_jet_params("./data-file-dir/");

	// // Spectrum time interval
	float tlo = 0.5; // The T90? essentially.
	float thi = 1.;

	// // Total spectrum
	test_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi);
	test_grb.WriteSpectrumToTXT("data-file-dir/synthGRB_spec_TOT.txt");

	// // Component spectrum
	test_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi, "TH");
	test_grb.WriteSpectrumToTXT("data-file-dir/synthGRB_spec_TH.txt");
	test_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi, "IS");
	test_grb.WriteSpectrumToTXT("data-file-dir/synthGRB_spec_IS.txt");
	test_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi, "FS");
	test_grb.WriteSpectrumToTXT("data-file-dir/synthGRB_spec_FS.txt");
	test_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi, "RS");
	test_grb.WriteSpectrumToTXT("data-file-dir/synthGRB_spec_RS.txt");

	// Total light curve
	test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, tmin, tmax, dt);
	test_grb.WriteLightCurveToTXT("data-file-dir/synthGRB_light_curve.txt");

	// Component light curves
	test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, tmin, tmax, dt, "TH", false);
	test_grb.WriteLightCurveToTXT("data-file-dir/synthGRB_light_curve_TH.txt");
	test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, tmin, tmax, dt, "IS", false);
	test_grb.WriteLightCurveToTXT("data-file-dir/synthGRB_light_curve_IS.txt");
	test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, tmin, tmax, dt, "FS", false);
	test_grb.WriteLightCurveToTXT("data-file-dir/synthGRB_light_curve_FS.txt");
	test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, tmin, tmax, dt, "RS", false);
	test_grb.WriteLightCurveToTXT("data-file-dir/synthGRB_light_curve_RS.txt");

	// Looking at after glows 
	// test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, 15, 4e5, dt, "FS", true);
	// test_grb.WriteLightCurveToTXT("data-file-dir/synthGRB_light_curve_afterglow_gbm.txt");
	// test_grb.make_source_light_curve(0.2, 10, 1e2, 1e6, dt, "FS", true);
	// test_grb.WriteLightCurveToTXT("data-file-dir/synthGRB_light_curve_afterglow_xrt.txt");
	// test_grb.make_source_light_curve(1e-3, 5e-3, 1e2, 1e6, dt, "FS", true);
	// test_grb.WriteLightCurveToTXT("data-file-dir/synthGRB_light_curve_afterglow_opt.txt");

	return 0;
}
