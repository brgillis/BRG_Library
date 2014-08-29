/**       @file summary_functions.hpp
 *
 *     Project: brg
 *        Path: /brg/vector/summary_functions.hpp
 *
 *  Created on: 29 Aug 2014
 *      Author: brg
 */

#ifndef SUMMARY_FUNCTIONS_HPP_
#define SUMMARY_FUNCTIONS_HPP_

#include <vector>

#include "brg/global.h"

#include "brg/vector/elementwise_functions.hpp"

namespace brgastro {

// Sum
#if (1)

template< typename T >
const T sum(const std::vector<T> &v)
{
	T result = 0;
	for(unsigned int i = 0; i < v.size(); i++)
	{
		result += v[i];
	}
	return result;
}

template< typename T >
const T sum(const T v)
{
	return v;
}

#endif // Sum

// Product
#if (1)

template< typename T >
const T product(const std::vector<T> &v)
{
	T result = 1;
	for(unsigned int i = 0; i < v.size(); i++)
	{
		result *= v[i];
	}
	return result;
}

template< typename T >
const T product(const T v)
{
	return v;
}

#endif // Product

// Mean
#if (1)

template< typename T >
const T mean(const std::vector<T> &v)
{
	if(v.size()==0) return 0;
	return sum(v)/v.size();
}

template< typename T >
const T mean(const T v)
{
	return v;
}

#endif // Mean

// Standard Deviation
#if (1)

template< typename T >
const T std(const std::vector<T> &v)
{
	if(v.size()<=1) return 0;

	return std::sqrt( divide(subtract(sum( square(v) ), square(sum(v)) ), v.size() ) );
}

template< typename T >
const T stddev(const std::vector<T> &v)
{
	return std(v);
}

template< typename T >
const T std(const T v)
{
	return 0;
}

template< typename T >
const T stddev(const T v)
{
	return 0;
}

#endif // Standard Deviation

// all_true
#if (1)
inline const bool all_true(const std::vector<bool> v)
{
	for(unsigned int i=0; i < v.size(); i++)
	{
		if(!(v[i])) return false;
	}
	return true;
}

inline const bool all_true(const bool v)
{
	return v;
}
#endif // all_true

// all_false
#if (1)
inline const bool all_false(const std::vector<bool> v)
{
	for(unsigned int i=0; i < v.size(); i++)
	{
		if(v[i]) return false;
	}
	return true;
}

inline const bool all_false(const bool v)
{
	return (!v);
}
#endif // all_false

// not_all_true
#if (1)
inline const bool not_all_true(const std::vector<bool> v)
{
	return !all_true(v);
}

inline const bool not_all_true(const bool v)
{
	return !v;
}
#endif // not_all_true

// not_all_false
#if (1)
inline const bool not_all_false(const std::vector<bool> v)
{
	return !all_false(v);
}

inline const bool not_all_false(const bool v)
{
	return v;
}
#endif // all_false

// some_true_some_false
#if (1)
inline const bool some_true_some_false(const std::vector<bool> v)
{
	return ( (!all_true(v)) && (!all_false(v)) );
}

inline const bool some_true_some_false(const bool v)
{
	return false;
}
#endif // some_true_some_false

} // namespace brgastro



#endif /* SUMMARY_FUNCTIONS_HPP_ */
