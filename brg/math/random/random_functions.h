/**********************************************************************\
  @file random_functions.h

 **********************************************************************

 Copyright (C) 2014  Bryan R. Gillis

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

\**********************************************************************/

// body file: random_functions.cpp

#ifndef _BRG_RANDOM_FUNCTIONS_H_INCLUDED_
#define _BRG_RANDOM_FUNCTIONS_H_INCLUDED_

#include <stdlib.h>

#include "brg/global.h"

namespace brgastro
{

/** Global function declarations **/
#if (1)

// Generates a random double between min and max
inline double drand()
{
	return drand48();
}
inline double drand( double min, double max )
{
	return min + (max-min)*drand48();
} // double drand(double min, double max)

// Returns a random variable from a Gaussian distribution
double Gaus_rand();
double Gaus_rand( double mean, double stddev = 0 );

// Returns a random variable from a Gaussian distribution in log space
// Note that "mean" here is the desired mean, NOT the peak of the function (which differ in log space). If you want to use
// the peak, simply use the standard Gaus_rand version instead.
double log10Gaus_rand(); // Implicit mean=1, stddev=1
double log10Gaus_rand( double mean, double stddev = 1 );

// Returns a random variable from a Rayleigh distribution
double Rayleigh_rand();
double Rayleigh_rand( double sigma );

// Returns a Poisson random variable.
int Pois_rand( double lambda=1 );

#endif // End global function declarations

}

#endif
