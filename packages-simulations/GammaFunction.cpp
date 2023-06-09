/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2022-02-28

Defines the necessary scripts to calculate the incomplete gamma function. 
This script follows the procedures outlined in Numerical Recipes by Press, Teukolsky, Vetterling, and Flannery.

*/

#include <math.h>
#include <cmath>
#include <iostream>

using namespace std;


#define ITMAX 100 // Maximum allowed number of iterations
#define EPS 3.0e-7 // Relative accuracy
#define FPMIN 1.0e-30 // Number near the smaller representable floating-point number


// Returns the value for the natural log of the gamma function, Gamma(a)
float gamma_func(float xx)
{
	double x,y,tmp,ser;
	static double cof[6] = {76.18009172947146,-86.50532032941677, 
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y = x = xx;
	tmp = x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for(j=0; j<=5; j++){ ser+= cof[j]/++y; }
	return exp(-tmp+log(2.5066282746310005*ser/x) );
}


// Returns the incomplete gamma function P(a,x)
float gamma_inc(float a, float x)
{
	void gser(float *gamser, float a, float x, float *gamf);
	void gcf(float *gammacf, float a, float x, float *gamf);
	// void nrerror(char error_text[]);
	float gamser, gammacf, gamf;

	if (x < 0.0 || a <= 0.0 ){ std::cout << "Invalid arguments in routine gammap." << std::endl; }
	if ( x < ( a + 1. ))
	{
		// Use series representation
		gser(&gamser, a, x, &gamf);
		return gamser;
	}
	else
	{
		// Use the continued fraction representation
		gcf(&gammacf, a, x, &gamf);
		return gammacf;
	}
}




// Returns the incomplete gamma function Gamma(a,x) evaluated by its series representation
// Also return Gamma(a) as gamf.
void gser(float *gamser, float a, float x, float *gamf)
{
	float gamma_func(float xx);
	// void nrerror(char error_text[]);
	int n;
	float sum,del,ap;

	(*gamf) = gamma_func(a);
	if (x <= 0.0 )
	{
		if (x < 0.0){ std::cout << "x less than 0 in routine gser." << std::endl;}
		(*gamser) = 0.0;
		return;
	}
	else
	{
		ap = a;
		del = sum = 1.0 / a;
		for ( n = 1; n<=ITMAX; n++)
		{
			++ap;
			del *= x/ap;
			sum += del;

			if ( fabs(del) < fabs(sum)*EPS)
			{
				// (*gamser) = sum*exp(-x + (a*log(x)) - (*gamf));
				(*gamser) = (*gamf) - sum*exp(-x + (a*log(x)));
				return;
			}
		}
		std::cout << " 'a' too large, ITMAX too small in routine gser." << std::endl;
		return;
	}
}

// Returns the complement of the incomplete gamma function, evaluated by is complement, using the Lentz's method
// Also return Gamma(a) as gamf.
void gcf(float *gammacf, float a, float x, float *gamf)
{
	float gamma_func(float xx);
	// void nrerror(char error_text[]);
	int i;
	float an, b, c, d, del, h;

	(*gamf) = gamma_func(a);
	b = x+1.0-a;
	c = 1.0 / FPMIN;
	d = 1.0 / b;
	h=d;
	for( i = 1; i<=ITMAX; i++)
	{
		an = -i*(i-a);
		b += 2.0;
		d = an*d+b;
		if (fabs(d) < FPMIN){ d = FPMIN; }
		c = b+an/c;
		if (fabs(c) < FPMIN){ c = FPMIN; }
		d = 1.0/d;
		del = d*c;
		h *= del;
		if(fabs(del - 1.0) < EPS ){ break; }
	}
	if (i > ITMAX){ std::cout << " 'a' too large, ITMAX too small in routine gcf." << std::endl;}
	// (*gammacf) = exp(-x+ (a*log(x)) - (*gamf))*h;
	(*gammacf) = exp(-x+ (a*log(x))) * h;
}