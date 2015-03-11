/**********************************************************************\
  @file safe_math.hpp

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

#ifndef _BRG_SAFE_MATH_HPP_INCLUDED_
#define _BRG_SAFE_MATH_HPP_INCLUDED_

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <type_traits>

#include "brg/container/is_container.hpp"

#include "brg/global.h"
#include "../container/is_eigen_container.hpp"

namespace brgastro {

// "Safe" functions - perform the operation specified, but will
// take necessary actions to ensure it won't crash the program
// if the argument is invalid (eg. taking square-root of a negative
// number).
template< class T >
const T safe_sqrt( const T & a )
{

#ifdef _BRG_WARN_FOR_SAFE_FUNCTIONS_TRIGGERED_
	if(a < 0)
	{
		std::cerr << "WARNING: safe_sqrt() prevented error from negative input.\n";
	}
#endif
	return std::sqrt( std::fabs( a ) );
}
inline double safe_sqrt( const int a ) // Special case for integers due to -INT_MIN > INT_MAX
{
	using std::sqrt;

	double res;

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
inline const Ta safe_pow( const Ta & a, const Tx & x )
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
inline const T1 safe_quad_sub( const T1 & v1, const T2 & v2 )
{
	return safe_sqrt( v1 * v1 - v2 * v2 );
}

// Safe_d is used a bit differently: Put it around the denominator to make
// sure you don't hit a divide-by-zero error.
template< class T,
typename std::enable_if<!brgastro::is_const_container<T>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T>::value,char>::type = 0>
const T safe_d( const T & a )
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

// Container specialization
template< class T,
typename std::enable_if<brgastro::is_const_container<T>::value,char>::type = 0>
T safe_d( T array )
{
	for( auto & a : array)
	{

		#ifdef _BRG_WARN_FOR_SAFE_FUNCTIONS_TRIGGERED_
		if( (a == 0) || isbad(a) )
		{
			std::cerr << "WARNING: safe_d() prevented error from zero input or bad input.\n";
		}
		#endif

		#ifdef _BRG_USE_UNITS_
		decltype(a) min_d = a; // So it'll have the right units
		#else
		decltype(a) min_d;
		#endif

		if(std::is_integral<decltype(a)>::value)
			min_d = 1;
		else
			min_d = MIN_DIVISOR;

		if ( isnan( a ) )
		{
			a = min_d;
			continue;
		}
		if ( isinf( a ) )
		{
			a = 1. / min_d;
			continue;
		}

		if (std::fabs(a)<min_d)
		{
			a = min_d;
			continue;
		}

	}

	return array;
}

// Eigen-like container specialization
template< class T,
typename std::enable_if<!brgastro::is_const_container<T>::value,char>::type = 0,
typename std::enable_if<brgastro::is_eigen_container<T>::value,char>::type = 0>
T safe_d( T array )
{
	auto p = array.data();
	for( decltype(array.size()) i = 0; i < array.size(); ++i)
	{

		#ifdef _BRG_WARN_FOR_SAFE_FUNCTIONS_TRIGGERED_
		if( (*p == 0) || isbad(*p) )
		{
			std::cerr << "WARNING: safe_d() prevented error from zero input or bad input.\n";
		}
		#endif

		typename std::decay<decltype(*p)>::type min_d;

		if(std::is_integral<decltype(min_d)>::value)
			min_d = 1;
		else
			min_d = MIN_DIVISOR;

		if ( isnan( *p ) )
		{
			*p = min_d;
			continue;
		}
		if ( isinf( *p ) )
		{
			*p = 1. / min_d;
			continue;
		}

		if (std::fabs(*p)<min_d)
		{
			*p = min_d;
			continue;
		}

		++p;
	}

	return array;
}

} // namespace brgastro

#endif /* _BRG_SAFE_MATH_HPP_INCLUDED_ */
