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
	float phot_energ_min;
	float phot_energ_max;
	int num_phot_bins;

	float chan_energ_min;
	float chan_energ_max;
	int num_chans;

	std::vector<float> phot_energ_lo;
	std::vector<float> phot_energ_mid;
	std::vector<float> phot_energ_hi;

	std::vector<float> chan_energ_lo;
	std::vector<float> chan_energ_mid;
	std::vector<float> chan_energ_hi;

	std::vector<float> n_grp;
	std::vector<float> f_chan;
	std::vector<float> n_chan;

	std::vector< vector<float> > prob_matrix;

	// Response member functions

	void set_phot_energ(float phot_energ_min, float phot_energ_max, int num_phot_bins, bool logscale = true);
	void set_chan_energ(float chan_energ_min,float chan_energ_max, int num_chans, bool logscale = true);

	void set_n_grp();
	void set_f_chan();
	void set_n_chan();

	void make_empty_resp();
	void Identity();
	void OverDeltaE(float alpha);
	int LoadRespFromFile(std::string file_name);

	void WriteToFits(std::string out_file_name[]);


// private:

};

#endif 
