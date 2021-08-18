/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-11

Header file for the Response.hpp class

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef RSPCLASS_H
#define RSPCLASS_H

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "fitsio.h"

// Import Custom Libraries

using namespace std;

// Declare Response class
class Response
{
public:
	// Response constructor
	Response();

	// Response member variables 
	float E_phot_min;
	float E_phot_max;
	int num_phot_bins;

	float E_chan_min;
	float E_chan_max;
	int num_chans;

	std::vector<float> ENERG_LO;
	std::vector<float> ENERG_MID;
	std::vector<float> ENERG_HI;

	std::vector<float> ECHAN_LO;
	std::vector<float> ECHAN_MID;
	std::vector<float> ECHAN_HI;

	std::vector<float> N_GRP;
	std::vector<float> F_CHAN;
	std::vector<float> N_CHAN;

	std::vector< vector<float> > MATRIX;

	// Response member functions

	void set_ENERG(float E_phot_min, float E_phot_max, int num_phot_bins, bool logscale = true);
	void set_ECHAN(float E_chan_min,float E_chan_max, int num_chans, bool logscale = true);

	void set_N_GRP();
	void set_F_CHAN();
	void set_N_CHAN();

	void make_empty_resp();
	void identity();
	void overDeltaE(float alpha);
	int load_rsp_from_file(char file_name[]);


// private:

};

#endif 
