/**********************************************************************\
  @file random_functions.hpp

 **********************************************************************

 Copyright (C) 2014, 2015  Bryan R. Gillis

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

#ifndef _BRG_RANDOM_FUNCTIONS_HPP_INCLUDED_
#define _BRG_RANDOM_FUNCTIONS_HPP_INCLUDED_

#include <cmath>
#include <random>

#include "IceBRG_main/common.h"
#include "IceBRG_main/math/misc_math.hpp"

namespace IceBRG
{

extern std::ranlux48 rng; // Initialised in random_functions.cpp

/** Global function declarations **/
#if (1)

// Generates a random int_t between min and max, inclusive of min, exclusive of max
template< typename T=int_t, typename T_in=int_t, typename T_gen=decltype(rng) >
T irand( T_in && min, T_in && max, T_gen & gen=rng )
{
	return std::uniform_int_distribution<T>(std::forward<T_in>(min),std::forward<T_in>(max)-1)(gen);
} // flt_t drand(flt_t min, flt_t max)

// Generates a random flt_t between min and max
template< typename T=flt_t, typename T_gen=decltype(rng) >
T drand( T_gen & gen=rng )
{
	return std::uniform_real_distribution<T>()(gen);
}
template< typename T=flt_t, typename T_in=flt_t, typename T_gen=decltype(rng) >
T drand( T_in && min, T_in && max, T_gen & gen=rng )
{
	return std::uniform_real_distribution<T>(std::forward<T_in>(min),std::forward<T_in>(max))(gen);
} // flt_t drand(flt_t min, flt_t max)

// Returns a random variable from a Gaussian distribution
template< typename T=flt_t, typename T_gen=decltype(rng) >
T Gaus_rand( T_gen & gen=rng )
{
	return std::normal_distribution<T>()(gen);
} // flt_t Gaus_rand()
template< typename T=flt_t, typename T_in=flt_t, typename T_gen=decltype(rng) >
T Gaus_rand( T_in && mean, T_in && stddev = 1.0, T_gen & gen=rng )
{
	return std::normal_distribution<T>(std::forward<T_in>(mean),std::forward<T_in>(stddev))(gen);

} // flt_t Gaus_rand(flt_t mean, flt_t stddev)

// Returns a random variable from a Gaussian distribution, truncated between min and max
template< typename T=flt_t, typename T_in=flt_t, typename T_gen=decltype(rng) >
T trunc_Gaus_rand( T_in && mean, T_in && stddev, T_in && min, T_in && max, T_gen & gen=rng )
{
	assert(max>min);

	// Try values until one works
	bool good_value = false;
	int attempt_counter = 0;

	while( (!good_value) and (attempt_counter < 1000) )
	{
		T test_result = Gaus_rand(mean,stddev,gen);
		if((test_result >= min)and(test_result <= max))
		{
			return test_result;
		}
		else
		{
			++attempt_counter;
		}
	}

	// Failsafe
	return (min+max)/2.;

} // T trunc_Gaus_rand( T_in && mean, T_in && stddev, T_in && min, T_in && max, gen_t & gen=rng )

// Returns a random variable from a Gaussian distribution in log space
// Note that "mean" here is the desired mean, NOT the peak of the function (which differ in log space). If you want to use
// the peak, simply use the standard Gaus_rand version instead.
template< typename T=flt_t, typename T_gen=decltype(rng) >
T log10Gaus_rand( T_gen & gen=rng )
{
	const flt_t & fact = std::exp( -square( std::log( 10. ) ) / 2 );

	return ( fact * std::pow(10., Gaus_rand<T>(gen) ) );
} // flt_t log10Gaus_rand()
template< typename T=flt_t, typename T_in=flt_t, typename T_gen=decltype(rng) >
T log10Gaus_rand( T_in && mean, T_in && stddev = 1., T_gen & gen=rng )
{
	const flt_t & fact = std::exp( -square( stddev*std::log(10.) ) / 2 );

	return ( fact * std::pow(10., Gaus_rand<T,T_in>(std::forward<T_in>(mean),std::forward<T_in>(stddev),gen) ) );
} // flt_t log10Gaus_rand(flt_t mean, flt_t stddev)

// Returns a random variable from a log10_Gaussian distribution, truncated between min and max
template< typename T=flt_t, typename T_in=flt_t, typename T_gen=decltype(rng) >
T trunc_log10Gaus_rand( T_in && mean, T_in && stddev, T_in && min, T_in && max, T_gen & gen=rng )
{
	assert(max>min);

	// Try values until one works
	bool good_value = false;
	int attempt_counter = 0;

	T p10_min = std::pow(10.,min);
	T p10_max = std::pow(10.,max);

	while( (!good_value) and (attempt_counter < 1000) )
	{
		T test_result = log10Gaus_rand(mean,stddev,gen);
		if((test_result >= p10_min)and(test_result <= p10_max))
		{
			return test_result;
		}
		else
		{
			++attempt_counter;
		}
	}

	// Failsafe
	return std::pow(10.,(min+max)/2.);

} // T trunc_Gaus_rand( T_in && mean, T_in && stddev, T_in && min, T_in && max, gen_t & gen=rng )

// Returns a random variable from a Rayleigh distribution
template< typename T=flt_t, typename T_in=flt_t, typename T_gen=decltype(rng) >
T Rayleigh_rand( T_in && sigma=1., T_gen & gen=rng )
{
	return std::forward<T_in>(sigma)*std::sqrt(-2.*std::log(drand<T>(gen)));
}

// Returns a random variable from a Gaussian distribution, truncated between min and max
template< typename T=flt_t, typename T_in=flt_t, typename T_gen=decltype(rng) >
T trunc_Rayleigh_rand( T_in && sigma, T_in && max, T_gen & gen=rng )
{
	assert(max>0.);

	// Try values until one works
	bool good_value = false;
	int attempt_counter = 0;

	while( (!good_value) and (attempt_counter < 1000) )
	{
		T test_result = Rayleigh_rand(sigma,gen);
		if(test_result <= max)
		{
			return test_result;
		}
		else
		{
			++attempt_counter;
		}
	}

	// Failsafe
	return max/2.;

} // T trunc_Gaus_rand( T_in && mean, T_in && stddev, T_in && min, T_in && max, gen_t & gen=rng )

// Get a Rayleigh random variable, smoothly contracted to always be less than max
template< typename T=flt_t, typename T_in=flt_t, typename T_gen=decltype(rng) >
T contracted_Rayleigh_rand( T_in && sigma, T_in && max, T_in && p, T_gen & gen=rng )
{
	// Generate an initial random Rayleigh variable
	T first_result = Rayleigh_rand(sigma);

	// Transform it via Bryan's formula to rein in large values to be less than the max_val
	return (first_result / std::pow(1 + std::pow(first_result / max, p), 1.0 / p));
}

// Returns a Poisson random variable.
template< typename T=int_t, typename T_gen=decltype(rng) >
T Pois_rand( T_gen & gen=rng )
{
	return std::poisson_distribution<T>()(gen);
} // T Pois_rand( T_gen & gen=rng )
template< typename T=int_t, typename T_in=flt_t, typename T_gen=decltype(rng) >
T Pois_rand( T_in && lambda=1., T_gen & gen=rng )
{
	return std::poisson_distribution<T>(std::forward<T_in>(lambda))(gen);
} // T Pois_rand( T_in && lambda=1., T_gen & gen=rng )

template< typename T=flt_t, typename T_gen=decltype(rng) >
T rand_from_cdf_arrays( array_t<T> const & xvals, array_t<T> cvals, T_gen & gen = rng )
{
	assert(xvals.size()>=2 and xvals.size()==cvals.size());

	// Quietly normalize the cvals
	cvals -= cvals[0];

	T cmax = cvals[cvals.size()-1];

	if(cmax<=0)
	{
		throw std::runtime_error("Invalid values used for generating random value: Final CDF value is <= 0.");
	}

	cvals /= cmax;

	// Generate a random value
	T const r = drand(0.,1.,rng);

	// Get the index on the cdf where this lies
	array_t<T> diffs = (cvals-r).abs();
	int_t i,j;
	diffs.minCoeff(&i,&j);

	// If the value at the index is below r, or the index is zero, move up one index
	while(((cvals[i]<r) and (i<cvals.size())) or (i==0))
	{
		++i;
	}

	// Interpolate to estimate the value
	T const clow = cvals[i - 1];
	T const chi = cvals[i];
	T const xlow = xvals[i - 1];
	T const xhi = xvals[i];

	T res = xlow + (xhi - xlow) / (chi - clow) * (r - clow);

	// Check for edge cases
	if( res < xvals[0] ) res = xvals[0];
	if( res > xvals[xvals.size()-1] ) res = xvals[xvals.size()-1];

	return res;

}

template< typename Tf, typename T=flt_t, typename T_gen=decltype(rng) >
T rand_from_cdf( Tf const & f, int_t const & N_samples=1000,
		T const & xlow=-5., T const & xhigh=5., T_gen & gen = rng )
{
	// Get an array of x points and cdf values at those points
	array_t<T> xvals = array_t<T>::LinSpaced(N_samples, xlow, xhigh);
	array_t<T> cvals = xvals.unaryExpr(f);

	T res = rand_from_cdf_arrays( xvals, cvals, gen );

	return res;
}

template< typename Tf, typename T=flt_t, typename T_gen=decltype(rng) >
T rand_from_pdf( Tf const & f, int_t const & N_samples=1000,
		T const & xlow=-5., T const & xhigh=5., T_gen & gen = rng )
{
	// Get an array of x points and pdf values at those points
	array_t<T> xvals = array_t<T>::LinSpaced(N_samples, xlow, xhigh);
	array_t<T> pvals = xvals.unaryExpr(f);

	// Get (unnormalized) cdf values
	array_t<T> cvals(pvals.size());

	std::partial_sum(pvals.data(), pvals.data()+pvals.size(), cvals.data(), std::plus<T>());

	T res = rand_from_cdf_arrays( xvals, cvals, gen );

	return res;
}

#endif // End global function declarations

}

#endif
