#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <tuple>
#include "fitsio.h"
#include <list>

using namespace std;

// Import Custom Libraries
#include "SynthGRB.hpp"
// #include "SynthGRBLibrary.hpp"
// #include "ObsGRB.hpp"
// #include "Response.hpp"
#include "LightCurve.hpp"
#include "Spectrum.hpp"
#include "cosmology.hpp"
#include "utilfuncs.hpp"
// #include "DataAnalysis.hpp"
// #include "TTEs.hpp"
#include "ModelParams.hpp"
#include "ShellDist.hpp"
#include "GammaFunction.hpp"

int main(int argc, char const *argv[])
{
	SynthGRB test_grb = SynthGRB();
	test_grb.LoadJetParamsFromTXT("input-files/jet-params.txt");
	
	(*test_grb.p_jet_shells).WriteToTXT("data-file-dir/shell_dist.txt");

	test_grb.SimulateJetDynamics();
	test_grb.write_out_jet_params("./data-file-dir/");

	
	return 0;
}

