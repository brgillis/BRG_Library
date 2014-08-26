/**********************************************************************\
brg_random_functions.h
 -----------

 If this header is used, the source file brg_functions.cpp must be included
 and compiled with the project. This file automatically includes
 brg_functions.hpp, which contains all template and inline functions.
 More complex functions are declared in this file and implemented in
 brg_functions.cpp.

 This file includes various classes and functions for general-purpose
 use. The file is split into two primary sections:

 -Class definitions
 -Global function declarations

 These sections are explained in further detail in their respective
 documentation blocks.

 Everything in this file is declared in the namespace brgastro.

 \**********************************************************************/

#ifndef __BRG_RANDOM_FUNCTIONS_H_INCLUDED__
#define __BRG_RANDOM_FUNCTIONS_H_INCLUDED__

#include <stdlib.h>

#include "brg_global.h"

namespace brgastro
{

/** Global function declarations **/
#if (1)

// Generates a random double between min and max
inline const double drand()
{
	return drand48();
}
inline const double drand( double min, double max )
{
	return min + (max-min)*drand48();
} // double drand(double min, double max)

// Returns a random variable from a Gaussian distribution
const double Gaus_rand();
const double Gaus_rand( const double mean, const double stddev = 0 );

// Returns a random variable from a Gaussian distribution in log space
// Note that "mean" here is the desired mean, NOT the peak of the function (which differ in log space). If you want to use
// the peak, simply use the standard Gaus_rand version instead.
const double log10Gaus_rand(); // Implicit mean=1, stddev=1
const double log10Gaus_rand( const double mean, const double stddev = 1 );

// Returns a random variable from a Rayleigh distribution
const double Rayleigh_rand();
const double Rayleigh_rand( const double sigma );

// Returns a Poisson random variable.
const int Pois_rand( const double lambda=1 );

#endif // End global function declarations

}

#endif
