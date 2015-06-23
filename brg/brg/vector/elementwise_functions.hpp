/**********************************************************************\
  @file elementwise_functions.hpp

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

#ifndef _BRG_ELEMENTWISE_FUNCTIONS_HPP_INCLUDED_
#define _BRG_ELEMENTWISE_FUNCTIONS_HPP_INCLUDED_

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <type_traits>
#include <vector>

#include "brg/common.h"

#include "brg/container/is_container.hpp"
#include "brg/container/is_eigen_container.hpp"
#include "brg/math/misc_math.hpp"
#include "brg/math/safe_math.hpp"
#include "brg/utility.hpp"

namespace brgastro {

// These apply a standard function to each element of a vector and return a vector of results

// Element-wise generic function
#if (1)

template< typename f, typename T >
T apply( const f & func, T v1)
{
	std::transform(v1.begin(), v1.end(), v1.begin(), func);

	return v1;
}

template< typename f, typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
T1 apply( const f & func,  T1 v1, const T2 & v2)
{
	std::transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), func);

	return v1;
}

template< typename f, typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0 >
T2 apply( const f & func, const T1 & v1, T2 v2)
{
	for(auto & v : v2) v = func(v1,v);

	return v2;
}

template< typename f, typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
T1 apply( const f & func, T1 v1, const T2 &v2 )
{
	for(auto & v : v1) v = func(v,v2);

	return v1;
}

template< typename f, typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0 >
T1 apply( const f & func, const T1 & v1, const T2 &v2 )
{
	T1 result = func(v1,v2);

	return result;
}

#endif // Element-wise generic function

// Random values
#if (1)

template<typename T, typename f>
T rand_vector_of_size(const f func, const int_type size)
{
	T result(size);

	for(auto & v : result) v = func();

	return result;
}

template<typename f, typename T1>
T1 rand_vector_of_size(const f func, const T1 & v1, const int_type size)
{
	T1 result(size);

	for(auto & v : result) v = func(v1);

	return result;
}

template<typename f, typename T1, typename T2>
T1 rand_vector_of_size(const f func, const T1 & v1, const T2 & v2, const int_type size)
{
	T1 result(size);

	for(auto & v : result) v = func(v1,v2);

	return result;
}

template<typename f, typename T1,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0>
T1 rand_vector(const f func, T1 v1 )
{
	for(typename T1::size_type i = 0; i < v1.size(); i++) v1[i] = func(v1[i]);

	return v1;
}

template<typename f, typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
T1 rand_vector(const f func, T1 v1, const T2 & v2)
{
	assert(v1.size()==v2.size());

	for(typename T1::size_type i = 0; i < v1.size(); i++) v1[i] = (func)(v1[i],v2[i]);

	return v1;
}

template<typename f, typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0 >
T1 rand_vector(const f func, T1 v1, const T2 & v2)
{
	for(auto & v : v1) v = func(v,v2);

	return v1;
}

template<typename f, typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
T2 rand_vector(const f func, const T1 & v1, T2 v2)
{
	for(auto & v : v2) v = func(v1,v);

	return v2;
}

#endif

// Element-wise math
#if (1)

// Element-wise addition
#if (1)

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto add( const T1 & v1, const T2 &v2 ) -> decltype(v1+v2)
{
	return v1+v2;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto add( const T1 & v1, const T2 &v2 ) -> std::vector<decltype(v1[0]+v2)>
{
	typedef decltype(v1[0]+v2) Tvout;
	typedef std::vector<Tvout> Tout;

	Tout vout(v1.size());

	for(ssize_t i=0; i<ssize(vout);++i)
	{
		vout[i] = add(v1[i],v2);
	}

	return vout;
}

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto add( const T1 & v1, const T2 &v2 ) -> std::vector<decltype(v1+v2[0])>
{
	return add(v2,v1);
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto add( const T1 & v1, const T2 &v2 ) -> std::vector<decltype(v1[0]+v2[0])>
{
	assert(v1.size()==v2.size());

	typedef decltype(v1[0]+v2[0]) Tvout;
	typedef std::vector<Tvout> Tout;

	Tout vout(v1.size());

	for(ssize_t i=0; i<ssize(vout);++i)
	{
		vout[i] = add(v1[i],v2[i]);
	}

	return vout;
}

#endif // Element-wise addition

// Element-wise subtraction
#if (1)

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto subtract( const T1 & v1, const T2 &v2 ) -> decltype(v1-v2)
{
	return v1-v2;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto subtract( const T1 & v1, const T2 &v2 ) -> std::vector<decltype(v1[0]-v2)>
{
	typedef decltype(v1[0]-v2) Tvout;
	typedef std::vector<Tvout> Tout;

	Tout vout(v1.size());

	for(ssize_t i=0; i<ssize(vout);++i)
	{
		vout[i] = subtract(v1[i],v2);
	}

	return vout;
}

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto subtract( const T1 & v1, const T2 &v2 ) -> std::vector<decltype(v1-v2[0])>
{
	typedef decltype(v1-v2[0]) Tvout;
	typedef std::vector<Tvout> Tout;

	Tout vout(v1.size());

	for(ssize_t i=0; i<ssize(vout);++i)
	{
		vout[i] = subtract(v1,v2[i]);
	}

	return vout;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto subtract( const T1 & v1, const T2 &v2 ) -> std::vector<decltype(v1[0]-v2[0])>
{
	assert(v1.size()==v2.size());

	typedef decltype(v1[0]-v2[0]) Tvout;
	typedef std::vector<Tvout> Tout;

	Tout vout(v1.size());

	for(ssize_t i=0; i<ssize(vout);++i)
	{
		vout[i] = subtract(v1[i],v2[i]);
	}

	return vout;
}

#endif // Element-wise subtraction

// Element-wise multiplication
#if (1)

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto multiply( const T1 & v1, const T2 &v2 ) -> decltype(v1*v2)
{
	return v1*v2;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto multiply( const T1 & v1, const T2 &v2 ) -> std::vector<decltype(v1[0]*v2)>
{
	typedef decltype(v1[0]*v2) Tvout;
	typedef std::vector<Tvout> Tout;

	Tout vout(v1.size());

	for(ssize_t i=0; i<ssize(vout);++i)
	{
		vout[i] = multiply(v1[i],v2);
	}

	return vout;
}

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto multiply( const T1 & v1, const T2 &v2 ) -> std::vector<decltype(v1*v2[0])>
{
	return multiply(v2,v1);
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto multiply( const T1 & v1, const T2 &v2 ) -> std::vector<decltype(v1[0]*v2[0])>
{
	assert(v1.size()==v2.size());

	typedef decltype(v1[0]*v2[0]) Tvout;
	typedef std::vector<Tvout> Tout;

	Tout vout(v1.size());

	for(ssize_t i=0; i<ssize(vout);++i)
	{
		vout[i] = multiply(v1[i],v2[i]);
	}

	return vout;
}

#endif // Element-wise multiplication

// Element-wise division
#if (1)

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto divide( const T1 & v1, const T2 &v2 ) -> decltype(v1/v2)
{
	return v1/v2;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto divide( const T1 & v1, const T2 &v2 ) -> std::vector<decltype(v1[0]/v2)>
{
	typedef decltype(v1[0]/v2) Tvout;
	typedef std::vector<Tvout> Tout;

	Tout vout(v1.size());

	for(ssize_t i=0; i<ssize(vout);++i)
	{
		vout[i] = divide(v1[i],v2);
	}

	return vout;
}

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto divide( const T1 & v1, const T2 &v2 ) -> std::vector<decltype(v1/v2[0])>
{
	typedef decltype(v1/v2[0]) Tvout;
	typedef std::vector<Tvout> Tout;

	Tout vout(v1.size());

	for(ssize_t i=0; i<ssize(vout);++i)
	{
		vout[i] = divide(v1,v2[i]);
	}

	return vout;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
auto divide( const T1 & v1, const T2 &v2 ) -> std::vector<decltype(v1[0]/v2[0])>
{
	assert(v1.size()==v2.size());

	typedef decltype(v1[0]/v2[0]) Tvout;
	typedef std::vector<Tvout> Tout;

	Tout vout(v1.size());

	for(ssize_t i=0; i<ssize(vout);++i)
	{
		vout[i] = divide(v1[i],v2[i]);
	}

	return vout;
}

#endif // Element-wise divide

// Element-wise power
#if (1)

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0 >
T1 pow( const T1 & v1, const T2 & v2 )
{
	return std::pow(v1, v2);
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0 >
T1 pow( T1 v1, const T2 &v2 )
{
	for(auto & v : v1)
	{
		v = pow(v, v2);
	}

	return v1;
}

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
T2 pow( const T1 & v1, T2 v2 )
{
	using std::pow;
	for(auto & v : v2)
	{
		v = pow(v1, v);
	}

	return v2;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
T1 pow( T1 v1, const T2 &v2 )
{
	assert(v1.size()==v2.size());
	for(typename T1::size_type i = 0; i < v1.size(); i++)
	{
		v1[i] = pow(v1[i], v2[i]);
	}

	return v1;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0 >
T1 runtime_ipow( T1 v1, const T2 &v2 )
{
	for(auto & v : v1)
	{
		v = runtime_ipow(v, v2);
	}

	return v1;
}

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
T1 runtime_ipow( const T1 & v1, const T2 & v2 )
{
	T1 result(v2.size());
	for(typename T2::size_type i = 0; i < v2.size(); i++)
	{
		result[i] = runtime_ipow(v1, v2[i]);
	}

	return v2;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
T1 runtime_ipow( T1 v1, const T2 &v2 )
{
	assert(v1.size()==v2.size());
	for(typename T1::size_type i = 0; i < v1.size(); i++)
	{
		v1[i] = runtime_ipow(v1[i], v2[i]);
	}

	return v1;
}

#endif // Element-wise power

// Element-wise safe power
#if (1)

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0 >
T1 safe_pow( T1 v1, const T2 &v2 )
{

	for(auto & v : v1)
	{
		v = safe_pow(v, v2);
	}

	return v1;
}

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
T2 safe_pow( const T1 & v1, T2 v2 )
{
	for(auto & v : v2)
	{
		v = safe_pow(v1, v);
	}

	return v2;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
T1 safe_pow( T1 v1, const T2 &v2 )
{

	assert(v1.size()==v2.size());
	for(int_type i = 0; i < v1.size(); i++)
	{
		v1[i] = safe_pow(v1[i], v2[i]);
	}

	return v1;
}

#endif // Element-wise safe power

// Element-wise max
#if (1)

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value || brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value || brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T1 max( const T1 & v1, const T2 &v2 )
{

	assert(v1.size()==v2.size());

	auto vr(v1);

	for(decltype(v1.size()) i = 0; i < vr.size(); i++)
	{
		vr[i] = max(vr[i], v2[i]);
	}

	return vr;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value || brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value || brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T1 max( T1 && v1, const T2 &v2 )
{

	assert(v1.size()==v2.size());

	for(decltype(v1.size()) i = 0; i < v1.size(); i++)
	{
		v1[i] = max(v1[i], v2[i]);
	}

	return v1;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T1 max( const T1 & v1, const T2 &v2 )
{

	auto vr(v1);

	for(auto & v : vr)
	{
		v = max(v, v2);
	}

	return vr;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T1 max( T1 && v1, const T2 &v2 )
{

	for(auto & v : v1)
	{
		v = max(v, v2);
	}

	return v1;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T1 max( const T1 & v1, const T2 &v2 )
{
	auto vr(v1);

	for(decltype(v1.size()) i = 0; i < vr.size(); i++)
	{
		vr[i] = max(vr[i], v2);
	}

	return vr;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T1 max( T1 && v1, const T2 &v2 )
{

	for(decltype(v1.size()) i = 0; i < v1.size(); i++)
	{
		v1[i] = max(v1[i], v2);
	}

	return v1;
}

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
T2 max( const T1 & v1, const T2 & v2 )
{
	auto vr(v2);

	for(auto & v : vr)
	{
		v = max(v1, v);
	}

	return vr;
}

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
T2 max( const T1 & v1, T2 && v2 )
{
	for(auto & v : v2)
	{
		v = max(v1, v);
	}

	return v2;
}

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T2 max( const T1 & v1, const T2 & v2 )
{
	auto vr(v2);

	for(decltype(v2.size()) i = 0; i < vr.size(); i++)
	{
		vr[i] = max(v1, vr[i]);
	}

	return vr;
}

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T2 max( const T1 & v1, T2 && v2 )
{
	for(decltype(v2.size()) i = 0; i < v2.size(); i++)
	{
		v2[i] = max(v1, v2[i]);
	}

	return v2;
}

#endif // Element-wise max

// Element-wise min
#if (1)

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value || brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value || brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T1 min( const T1 & v1, const T2 &v2 )
{

	assert(v1.size()==v2.size());

	auto vr(v1);

	for(decltype(v1.size()) i = 0; i < vr.size(); i++)
	{
		vr[i] = min(vr[i], v2[i]);
	}

	return vr;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value || brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value || brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T1 min( T1 && v1, const T2 &v2 )
{

	assert(v1.size()==v2.size());

	for(decltype(v1.size()) i = 0; i < v1.size(); i++)
	{
		v1[i] = min(v1[i], v2[i]);
	}

	return v1;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T1 min( const T1 & v1, const T2 &v2 )
{

	auto vr(v1);

	for(auto & v : vr)
	{
		v = min(v, v2);
	}

	return vr;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T1 min( T1 && v1, const T2 &v2 )
{

	for(auto & v : v1)
	{
		v = min(v, v2);
	}

	return v1;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T1 min( const T1 & v1, const T2 &v2 )
{
	auto vr(v1);

	for(decltype(v1.size()) i = 0; i < vr.size(); i++)
	{
		vr[i] = min(vr[i], v2);
	}

	return vr;
}

template< typename T1, typename T2,
typename std::enable_if<brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_stl_container<T2>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T1 min( T1 && v1, const T2 &v2 )
{

	for(decltype(v1.size()) i = 0; i < v1.size(); i++)
	{
		v1[i] = min(v1[i], v2);
	}

	return v1;
}

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
T2 min( const T1 & v1, const T2 & v2 )
{
	auto vr(v2);

	for(auto & v : vr)
	{
		v = min(v1, v);
	}

	return vr;
}

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
T2 min( const T1 & v1, T2 && v2 )
{
	for(auto & v : v2)
	{
		v = min(v1, v);
	}

	return v2;
}

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T2 min( const T1 & v1, const T2 & v2 )
{
	auto vr(v2);

	for(decltype(v2.size()) i = 0; i < vr.size(); i++)
	{
		vr[i] = min(v1, vr[i]);
	}

	return vr;
}

template< typename T1, typename T2,
typename std::enable_if<!brgastro::is_stl_container<T1>::value,char>::type = 0,
typename std::enable_if<!brgastro::is_eigen_container<T1>::value,char>::type = 0,
typename std::enable_if<brgastro::is_eigen_container<T2>::value,char>::type = 0 >
T2 min( const T1 & v1, T2 && v2 )
{
	for(decltype(v2.size()) i = 0; i < v2.size(); i++)
	{
		v2[i] = min(v1, v2[i]);
	}

	return v2;
}

#endif // Element-wise min

// Element-wise bound
#if (1)

template< typename T1, typename T2, typename T3,
typename std::enable_if<brgastro::is_stl_container<T1>::value || brgastro::is_eigen_container<T1>::value ||
						brgastro::is_stl_container<T2>::value || brgastro::is_eigen_container<T2>::value ||
						brgastro::is_stl_container<T3>::value || brgastro::is_eigen_container<T3>::value ,char >::type = 0 >
T2 bound( T1 && lower_bound, T2 && a, T3 && upper_bound)
{
	return brgastro::min( brgastro::max( std::forward<T2>(a), std::forward<T1>(lower_bound) ) , std::forward<T3>(upper_bound) );
}

#endif // Element-wise bound

// Element-wise negate
#if (1)

template< typename T,
typename std::enable_if<brgastro::is_stl_container<T>::value,T>::type* = nullptr>
T negate( T v )
{
	for(auto & vv : v)
	{
		vv = negate(vv);
	}

	return v;
}

template< typename T,
typename std::enable_if<!brgastro::is_stl_container<T>::value,T>::type* = nullptr>
T negate( const T & v )
{
	return -v;
}

#endif // Element-wise negate

// Element-wise abs
#if (1)

template< typename T,
typename std::enable_if<brgastro::is_stl_container<T>::value,T>::type* = nullptr>
T abs( T v )
{
	using std::abs;

	for(auto & vv : v)
	{
		vv = abs(vv);
	}

	return v;
}

#endif // Element-wise abs

// Element-wise fabs
#if (1)

template< typename T,
typename std::enable_if<brgastro::is_stl_container<T>::value,T>::type* = nullptr>
T fabs( T v )
{
	using std::fabs;

	for(auto & vv : v)
	{
		vv = fabs(vv);
	}

	return v;
}

template< typename T,
typename std::enable_if<!brgastro::is_stl_container<T>::value,T>::type* = nullptr>
T fabs( T v )
{
	using std::fabs;

	return fabs(v);
}

#endif // Element-wise fabs

// Element-wise square root
#if (1)

template< typename T,
typename std::enable_if<brgastro::is_stl_container<T>::value,T>::type* = nullptr>
T sqrt( T v )
{
	using std::sqrt;

	for(auto & vv : v)
	{
		vv = sqrt(vv);
	}

	return v;
}

#endif // Element-wise square root

// Element-wise safe square root
#if (1)

template< typename T,
typename std::enable_if<brgastro::is_stl_container<T>::value,T>::type* = nullptr>
T safe_sqrt( T v )
{
	for(auto & vv : v)
	{
		vv = safe_sqrt(vv);
	}

	return v;
}

#endif // Element-wise safe square root

// Element-wise exponential
#if (1)

template< typename T,
typename std::enable_if<brgastro::is_stl_container<T>::value,T>::type* = nullptr>
T exp( T v )
{
	using std::exp;

	for(auto & vv : v)
		vv = exp(vv);

	return v;
}

#endif // Element-wise exponential

// Element-wise log
#if (1)

template< typename T,
typename std::enable_if<brgastro::is_stl_container<T>::value,T>::type* = nullptr>
T log( T v )
{
	using std::log;

	for(auto & vv : v)
		vv = log(vv);

	return v;
}

#endif // Element-wise log

// Element-wise square
#if (1)

template< typename T,
typename std::enable_if<brgastro::is_stl_container<T>::value,T>::type* = nullptr >
T square( T v )
{
	for(auto & vv : v)
	{
		vv = square(vv);
	}

	return v;
}

#endif // Element-wise square

// Element-wise cube
#if (1)

template< typename T,
typename std::enable_if<brgastro::is_stl_container<T>::value,T>::type* = nullptr >
T cube( T v )
{
	for(auto & vv : v)
	{
		vv = cube(vv);
	}

	return v;
}

#endif // Element-wise cube

// Element-wise quart
#if (1)

template< typename T,
typename std::enable_if<brgastro::is_stl_container<T>::value,T>::type* = nullptr >
T quart( T v )
{
	for(auto & vv : v)
	{
		vv = quart(vv);
	}

	return v;
}

#endif // Element-wise quart

// Element-wise inverse
#if (1)

template< typename T,
typename std::enable_if<brgastro::is_stl_container<T>::value,T>::type* = nullptr >
T inverse( T v )
{
	for(auto & vv : v)
	{
		vv = inverse(vv);
	}

	return v;
}

#endif // Element-wise inverse

// Element-wise inv_square
#if (1)

template< typename T,
typename std::enable_if<brgastro::is_stl_container<T>::value,T>::type* = nullptr >
T inv_square( T v )
{
	for(auto & vv : v)
	{
		vv = inv_square(vv);
	}

	return v;
}

#endif // Element-wise inv_square

// Element-wise inv_cube
#if (1)

template< typename T,
typename std::enable_if<brgastro::is_stl_container<T>::value,T>::type* = nullptr >
T inv_cube( T v )
{
	for(auto & vv : v)
	{
		vv = inv_cube(vv);
	}

	return v;
}

#endif // Element-wise inv_cube

// Element-wise inv_quart
#if (1)

template< typename T,
typename std::enable_if<brgastro::is_stl_container<T>::value,T>::type* = nullptr >
T inv_quart( T v )
{
	for(auto & vv : v)
	{
		vv = inv_quart(vv);
	}

	return v;
}

#endif // Element-wise inv_quart

// Element-wise safe_d
#if (1)

template< typename T,
typename std::enable_if<brgastro::is_stl_container<T>::value,T>::type* = nullptr >
T safe_d( T v )
{
	for(auto & vv : v)
		vv = safe_d(vv);

	return v;
}

#endif // Element-wise safe_d

#endif // Element-wise math

// Element-wise comparison (always return a vector of bools)
#if (1)

// Element-wise not
#if (1)

inline const std::vector<bool> v_not( std::vector<bool> v )
{
	for(std::vector<bool>::size_type i=0; i<v.size(); ++i)
	{
		v[i] = !v[i];
	}
	return v;
}

#endif // Element-wise not

// Element-wise equal
#if (1)
template< typename T1, typename T2 >
std::vector<bool> equal( const T1 & v1, const T2 &v2 )
{
	assert(v1.size()==v2.size());

	std::vector<bool> result(v1.size());

	for(typename T1::size_type i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] == v2[i]);
	}

	return result;
}

template< typename T1, typename T2, typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
std::vector<bool> equal( const T1 & v1, const T2 &v2 )
{
	std::vector<bool> result(v1.size());

	for(typename T1::size_type i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] == v2);
	}

	return result;
}

template< typename T1, typename T2 >
std::vector<bool> equal( const T2 & v1, const T1 &v2 )
{
	std::vector<bool> result(v2.size());

	for(typename T2::size_type i = 0; i < v2.size(); i++)
	{
		result[i] = (v2[i] == v1);
	}

	return result;
}

template< typename T1, typename T2 >
bool equal( const T1 & v1, const T2 & v2 )
{
	return (v1 == v2);
}

#endif // Element-wise equal

// Element-wise not_equal
#if (1)
template< typename T1, typename T2 >
std::vector<bool> not_equal( const T1 & v1, const T2 &v2 )
{
	assert(v1.size()==v2.size());

	std::vector<bool> result(v1.size());

	for(typename T1::size_type i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] != v2[i]);
	}

	return result;
}

template< typename T1, typename T2, typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
std::vector<bool> not_equal( const T1 & v1, const T2 &v2 )
{
	std::vector<bool> result(v1.size());

	for(typename T1::size_type i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] != v2);
	}

	return result;
}

template< typename T1, typename T2 >
std::vector<bool> not_equal( const T2 & v1, const T1 &v2 )
{
	std::vector<bool> result(v2.size());

	for(typename T2::size_type i = 0; i < v2.size(); i++)
	{
		result[i] = (v2[i] != v1);
	}

	return result;
}

template< typename T1, typename T2 >
bool not_equal( const T1 & v1, const T2 & v2 )
{
	return (v1 != v2);
}

#endif // Element-wise not_equal

// Element-wise less_than
#if (1)
template< typename T1, typename T2 >
std::vector<bool> less_than( const T1 & v1, const T2 &v2 )
{
	assert(v1.size()==v2.size());

	std::vector<bool> result(v1.size());

	for(typename T1::size_type i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] < v2[i]);
	}

	return result;
}

template< typename T1, typename T2, typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
std::vector<bool> less_than( const T1 & v1, const T2 &v2 )
{
	std::vector<bool> result(v1.size());

	for(typename T1::size_type i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] < v2);
	}

	return result;
}

template< typename T1, typename T2 >
std::vector<bool> less_than( const T2 & v1, const T1 &v2 )
{
	std::vector<bool> result(v2.size());

	for(typename T2::size_type i = 0; i < v2.size(); i++)
	{
		result[i] = (v1 < v2[i]);
	}

	return result;
}

template< typename T1, typename T2 >
bool less_than( const T1 & v1, const T2 & v2 )
{
	return (v1 < v2);
}

#endif // Element-wise less_than

// Element-wise greater_than
#if (1)
template< typename T1, typename T2 >
std::vector<bool> greater_than( const T1 & v1, const T2 &v2 )
{
	assert(v1.size()==v2.size());

	std::vector<bool> result(v1.size());

	for(typename T1::size_type i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] > v2[i]);
	}

	return result;
}

template< typename T1, typename T2, typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
std::vector<bool> greater_than( const T1 & v1, const T2 &v2 )
{
	std::vector<bool> result(v1.size());

	for(typename T1::size_type i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] > v2);
	}

	return result;
}

template< typename T1, typename T2 >
std::vector<bool> greater_than( const T2 & v1, const T1 &v2 )
{
	std::vector<bool> result(v2.size());

	for(typename T2::size_type i = 0; i < v2.size(); i++)
	{
		result[i] = (v1 > v2[i]);
	}

	return result;
}

template< typename T1, typename T2 >
bool greater_than( const T1 & v1, const T2 & v2 )
{
	return (v1 > v2);
}

#endif // Element-wise equal

// Element-wise less_than_or_equal
#if (1)
template< typename T1, typename T2 >
std::vector<bool> less_than_or_equal( const T1 & v1, const T2 &v2 )
{
	assert(v1.size()==v2.size());

	std::vector<bool> result(v1.size());

	for(typename T1::size_type i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] <= v2[i]);
	}

	return result;
}

template< typename T1, typename T2, typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
std::vector<bool> less_than_or_equal( const T1 & v1, const T2 &v2 )
{
	std::vector<bool> result(v1.size());

	for(typename T1::size_type i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] <= v2);
	}

	return result;
}

template< typename T1, typename T2 >
std::vector<bool> less_than_or_equal( const T2 & v1, const T1 &v2 )
{
	std::vector<bool> result(v2.size());

	for(typename T1::size_type i = 0; i < v2.size(); i++)
	{
		result[i] = (v1 <= v2[i]);
	}

	return result;
}

template< typename T1, typename T2 >
bool less_than_or_equal( const T1 & v1, const T2 & v2 )
{
	return (v1 <= v2);
}

#endif // Element-wise less_than_or_equal

// Element-wise greater_than_or_equal
#if (1)
template< typename T1, typename T2 >
std::vector<bool> greater_than_or_equal( const T1 & v1, const T2 &v2 )
{
	assert(v1.size()==v2.size());

	std::vector<bool> result(v1.size());

	for(typename T1::size_type i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] >= v2[i]);
	}

	return result;
}

template< typename T1, typename T2, typename std::enable_if<brgastro::is_stl_container<T2>::value,char>::type = 0 >
std::vector<bool> greater_than_or_equal( const T1 & v1, const T2 &v2 )
{
	std::vector<bool> result(v1.size());

	for(typename T1::size_type i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] >= v2);
	}

	return result;
}

template< typename T1, typename T2 >
std::vector<bool> greater_than_or_equal( const T2 & v1, const T1 &v2 )
{
	std::vector<bool> result(v2.size());

	for(typename T2::size_type i = 0; i < v2.size(); i++)
	{
		result[i] = (v1 >= v2[i]);
	}

	return result;
}

template< typename T1, typename T2 >
bool greater_than_or_equal( const T1 & v1, const T2 & v2 )
{
	return (v1 >= v2);
}

#endif // Element-wise greater_than_or_equal

#endif // Element-wise comparison

// Functions on bool vectors
#if (1)

// not
#if (1)
template<typename T>
inline const T v_not(const T v)
{
	for(typename T::size_type i=0; i < v.size(); i++)
	{
		v[i] = v_not(v[i]);
	}
	return v;
}

inline const bool v_not(const bool v)
{
	return !v;
}
#endif // not

#endif // Functions on bool vectors


} // namespace brgastro

#endif // _BRG_ELEMENTWISE_FUNCTIONS_HPP_INCLUDED_
