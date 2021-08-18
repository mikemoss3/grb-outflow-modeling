/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-11

Script to make response function

*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <random>
#include "fitsio.h"

using namespace std;

// Include custom libraries and headers
#include "Response.hpp"
#include "response_funcs.hpp"
#include "Spectrum.hpp"
#include "spectrum_funcs.hpp"

int main(int argc, char const *argv[])
{

	char file_name[] = "/Users/mjmoss/Research/instr-rsp-investigate/example-instr-files/sw00883832000b_1speak.rsp"
	Response test_rsp = Response();
	test_rsp.load_rsp_from_file(file_name);


	return 0;
}