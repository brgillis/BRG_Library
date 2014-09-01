/**********************************************************************\
  @file random_functions.cpp

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

#include <cstdlib>
#include <cmath>

#include "brg/global.h"

#include "brg/math/misc_math.hpp"
#include "brg/math/random/random_functions.h"
#include "brg/utility.hpp"

using namespace std;

/** Global function implementations **/
#if (1)

const double brgastro::Gaus_rand()
{
	double x1, x2, w;

	do
	{
		x1 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		x2 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( ( -2.0 * log( w ) ) / w );
	return x1 * w;
}

const double brgastro::Gaus_rand( double mean, double stddev )
{
	if ( stddev <= 0 )
		return mean;
	return ( mean + Gaus_rand() * stddev );

} // double Gaus_rand(double mean, double stddev)

const double brgastro::log10Gaus_rand()
{
	double x1, x2, w, fact;

	fact = exp( -square( log( 10 ) ) / 2 );

	do
	{
		x1 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		x2 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( ( -2.0 * log( w ) ) / w );
	return ( fact * pow(10., x1 * w ) );
} // double lnGaus_rand(double mean, double stddev)
const double brgastro::log10Gaus_rand( double mean, double stddev )
{
	double x1, x2, w, fact;

	if ( stddev <= 0 )
		return mean;

	stddev *= log( 10 ); // Converts dex to natural log
	fact = exp( -square( stddev ) / 2 );

	do
	{
		x1 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		x2 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( ( -2.0 * log( w ) ) / w );
	return ( mean * fact * exp( x1 * w * stddev ) );
} // double lnGaus_rand(double mean, double stddev)

// Returns a random variable from a Rayleigh distribution
const double brgastro::Rayleigh_rand()
{
	return sqrt(-2.*log(drand()));
}
const double brgastro::Rayleigh_rand( const double sigma )
{
	return sigma*Rayleigh_rand();
}

const int brgastro::Pois_rand( double lambda )
{
	double L, p;
	int k;
	L = exp( -lambda );
	k = 0;
	p = 1;
	do
	{
		k++;
		p *= drand( 0, 1 );
	} while ( p > L );

	return ( k - 1 );
} // int Pois_rand(double lambda)

#endif // end global function implementations
