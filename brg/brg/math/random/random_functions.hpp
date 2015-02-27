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

#include <cmath>
#include <stdlib.h>

#include "brg/global.h"
#include "brg/math/misc_math.hpp"

namespace brgastro
{

/** Global function declarations **/
#if (1)

// Generates a random double between min and max
inline double drand()
{
	return ( (double)rand() ) / RAND_MAX ;
}
inline double drand( const double & min, const double & max )
{
	return min + (max-min)*drand();
} // double drand(double min, double max)

// Returns a random variable from a Gaussian distribution
inline double Gaus_rand()
{
	double x1, x2, w;

	do
	{
		x1 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		x2 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = std::sqrt( ( -2.0 * std::log( w ) ) / w );
	return x1 * w;
} // double Gaus_rand()
inline double Gaus_rand( const double & mean, const double & stddev = 1 )
{
	if ( stddev <= 0 )
		return mean;
	return ( mean + Gaus_rand() * stddev );

} // double Gaus_rand(double mean, double stddev)

// Returns a random variable from a Gaussian distribution in log space
// Note that "mean" here is the desired mean, NOT the peak of the function (which differ in log space). If you want to use
// the peak, simply use the standard Gaus_rand version instead.
inline double log10Gaus_rand()
{
	double x1, x2, w, fact;

	fact = std::exp( -square( std::log( 10 ) ) / 2 );

	do
	{
		x1 = 2.0 * ( ( (double)std::rand() ) / RAND_MAX ) - 1.0;
		x2 = 2.0 * ( ( (double)std::rand() ) / RAND_MAX ) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = std::sqrt( ( -2.0 * std::log( w ) ) / w );
	return ( fact * pow(10., x1 * w ) );
} // double lnGaus_rand(double mean, double stddev)
inline double log10Gaus_rand( const double & mean, double stddev = 1 )
{
	double x1, x2, w, fact;

	if ( stddev <= 0 )
		return mean;

	stddev *= std::log( 10 ); // Converts dex to natural log
	fact = std::exp( -square( stddev ) / 2 );

	do
	{
		x1 = 2.0 * ( ( (double)std::rand() ) / RAND_MAX ) - 1.0;
		x2 = 2.0 * ( ( (double)std::rand() ) / RAND_MAX ) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = std::sqrt( ( -2.0 * std::log( w ) ) / w );
	return ( mean * fact * std::exp( x1 * w * stddev ) );
} // double lnGaus_rand(double mean, double stddev)

// Returns a random variable from a Rayleigh distribution
inline double Rayleigh_rand()
{
	return std::sqrt(-2.*std::log(drand()));
}
inline double Rayleigh_rand( const double & sigma )
{
	return sigma*Rayleigh_rand();
}

// Returns a Poisson random variable.
inline int Pois_rand( const double & lambda=1 )
{
	double L, p;
	int k;
	L = std::exp( -lambda );
	k = 0;
	p = 1;
	do
	{
		k++;
		p *= drand( 0, 1 );
	} while ( p > L );

	return ( k - 1 );
} // int Pois_rand(double lambda)

#endif // End global function declarations

}

#endif
