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
// #include "utilfuncs.hpp"
// #include "DataAnalysis.hpp"
// #include "TTEs.hpp"
// #include "ModelParams.hpp"
// #include "ShellDist.hpp"
// #include "GammaFunction.hpp"

using namespace std;

int main(int argc, char const *argv[])
{

	float en = 20.;
	float gamma = 100;
	float temp = gamma*20./kb_kev;
	float alpha=0.4;
	float norm = 1.;
	double param_list[3] = {temp, alpha, norm};
	float val = BB(en, param_list);
	
	std::cout << val << endl;

	return 0;
}