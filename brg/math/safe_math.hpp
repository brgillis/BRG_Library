/**       @file safe_math.hpp
 *
 *     Project: brg
 *        Path: /brg/math/safe_math.hpp
 *
 *  Created on: 29 Aug 2014
 *      Author: brg
 */

#ifndef _BRG_SAFE_MATH_HPP_INCLUDED_
#define _BRG_SAFE_MATH_HPP_INCLUDED_

#include <cstdlib>
#include <iostream>
#include <limits>

#include "brg/global.h"

namespace brgastro {

// "Safe" functions - perform the operation specified, but will
// take necessary actions to ensure it won't crash the program
// if the argument is invalid (eg. taking square-root of a negative
// number).
template< class T >
const T safe_sqrt( const T a )
{

#ifdef _BRG_WARN_FOR_SAFE_FUNCTIONS_TRIGGERED_
	if(a < 0)
	{
		std::cerr << "WARNING: safe_sqrt() prevented error from negative input.\n";
	}
#endif
	return std::sqrt( std::fabs( a ) );
}
inline const double safe_sqrt( const int a ) // Special case for integers due to -INT_MIN > INT_MAX
{
	using std::sqrt;

	int res;

#ifdef _BRG_WARN_FOR_SAFE_FUNCTIONS_TRIGGERED_
	if(a < 0)
	{
		std::cerr << "WARNING: safe_sqrt() prevented error from negative input.\n";
	}
#endif

	if ( a == std::numeric_limits<int>::min() )
	{
		res = std::numeric_limits<int>::max();
	}
	else
	{
		res = a < 0 ? -a : a;
	}
	return sqrt( res );
}
template< class Ta, class Tx >
inline const Ta safe_pow( const Ta a, const Tx x )
{
	Ta res;
	double ipart;

	std::modf( a, &ipart );

	if ( ( a < 0 ) && ( ipart != a ) )
	{
#ifdef _BRG_WARN_FOR_SAFE_FUNCTIONS_TRIGGERED_
		std::cerr << "WARNING: safe_pow() prevented error from negative input.\n";
#endif
		res = -a;
	}
	else
	{
		res = a;
	}

	return std::pow( res, x );
}

// Subtract one value from another in quadrature
template < typename T1, typename T2 >
inline const T1 safe_quad_sub( const T1 v1, const T2 v2 )
{
	return safe_sqrt( v1 * v1 - v2 * v2 );
}

// Safe_d is used a bit differently: Put it around the denominator to make
// sure you don't hit a divide-by-zero error.
template< class T >
const T safe_d( const T a )
{

#ifdef _BRG_WARN_FOR_SAFE_FUNCTIONS_TRIGGERED_
	if( (a == 0) || isbad(a) )
	{
		std::cerr << "WARNING: safe_d() prevented error from zero input or bad input.\n";
	}
#endif
	T min_d = a; // So it'll have the right units if we're using units here.
	min_d = MIN_DIVISOR;

	if ( isnan( a ) )
		return min_d;
	if ( isinf( a ) )
		return 1. / min_d;

	if (std::fabs(a)>min_d)
		return a;

	if(min_d == 0) return 1; // In case of integers

	return min_d;

}

} // namespace brgastro

#endif /* _BRG_SAFE_MATH_HPP_INCLUDED_ */
