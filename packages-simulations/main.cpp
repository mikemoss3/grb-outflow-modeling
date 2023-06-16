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
	// Directory path to data storage directory
	std::string dir="../files-data/synthetic-data/";

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
	test_grb.LoadJetParamsFromTXT("files-input/jet-params.txt");
	
	test_grb.anim_lor_dist = true;
	(*test_grb.p_jet_shells).WriteToTXT(dir+"synthGRB_shell_dist.txt");

	test_grb.SimulateJetDynamics();
	test_grb.write_out_jet_params(dir);

	// // Spectrum time interval
	float tlo = 0.5; // The T90? essentially.
	float thi = 1.;

	// // Total spectrum
	test_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi);
	test_grb.WriteSpectrumToTXT(dir+"synthGRB_spec_TOT.txt");

	// // Component spectrum
	test_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi, "TH");
	test_grb.WriteSpectrumToTXT(dir+"synthGRB_spec_TH.txt");
	test_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi, "IS");
	test_grb.WriteSpectrumToTXT(dir+"synthGRB_spec_IS.txt");
	test_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi, "FS");
	test_grb.WriteSpectrumToTXT(dir+"synthGRB_spec_FS.txt");
	test_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi, "RS");
	test_grb.WriteSpectrumToTXT(dir+"synthGRB_spec_RS.txt");

	// Total light curve
	test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, tmin, tmax, dt);
	test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve.txt");

	// Component light curves
	test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, tmin, tmax, dt, "TH", false);
	test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve_TH.txt");
	test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, tmin, tmax, dt, "IS", false);
	test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve_IS.txt");
	test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, tmin, tmax, dt, "FS", false);
	test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve_FS.txt");
	test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, tmin, tmax, dt, "RS", false);
	test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve_RS.txt");

	// Looking at after glows 
	// test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, 15, 4e5, dt, "FS", true);
	// test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve_afterglow_gbm.txt");
	// test_grb.make_source_light_curve(0.2, 10, 1e2, 1e6, dt, "FS", true);
	// test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve_afterglow_xrt.txt");
	// test_grb.make_source_light_curve(1e-3, 5e-3, 1e2, 1e6, dt, "FS", true);
	// test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve_afterglow_opt.txt");

	return 0;
}
