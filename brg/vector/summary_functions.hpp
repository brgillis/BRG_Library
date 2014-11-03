/**********************************************************************\
  @file summary_functions.hpp

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

	return std::sqrt( subtract(sum( square(v) )/v.size(), square(sum(v)/v.size()) ) );
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

// Sample standard error
#if(1)

template< typename T >
const T stderr(const std::vector<T> &v)
{
	if(v.size()<=2) return std::numeric_limits<T>::max();

	return std(v)/std::sqrt(v.size()-1);
}

template< typename T >
const T stderr(const T &v)
{
	return std::numeric_limits<T>::max();
}

#endif

// Monotonically increasing
#if(1)

template<typename InputIterator>
bool is_monotonically_increasing(const InputIterator &first, const InputIterator &last)
{
	typedef decltype(*first) T;

	// Check if container is empty
	if(first==last) return false;

	T last_value = *first;
	InputIterator it=first;

	// Test the first step, to make sure the container is large enough
	++it;
	if(it==last)
	{
		// If we're here, the container only has one element
		return false;
	}
	if(*it<=last_value) return false;
	last_value = *it;

	// Now do remainder of steps
	for(++it; it!=last; ++it)
	{
		if(*it<=last_value) return false;
		last_value = *it;
	}
	return true;
}

template<typename T>
bool is_monotonically_increasing(const std::vector<T> &v)
{
	return is_monotonically_increasing(v.begin(),v.end());
}

template<typename T>
bool is_monotonically_increasing(T v)
{
	return false;
}

#endif // Monotonically increasing

// Monotonically decreasing
#if(1)

template<typename InputIterator>
bool is_monotonically_decreasing(const InputIterator &first, const InputIterator &last)
{
	typedef decltype(*first) T;

	// Check if container is empty
	if(first==last) return false;

	T last_value = *first;
	InputIterator it=first;

	// Test the first step, to make sure the container is large enough
	++it;
	if(it==last)
	{
		// If we're here, the container only has one element
		return false;
	}
	if(*it>=last_value) return false;
	last_value = *it;

	// Now do remainder of steps
	for(++it; it!=last; ++it)
	{
		if(*it>=last_value) return false;
		last_value = *it;
	}
	return true;
}

template<typename T>
bool is_monotonically_decreasing(const std::vector<T> &v)
{
	return is_monotonically_decreasing(v.begin(),v.end());
}

template<typename T>
bool is_monotonically_decreasing(T v)
{
	return false;
}

#endif // Monotonically decreasing

// all_true
#if (1)
inline bool all_true(const std::vector<bool> v)
{
	for(unsigned int i=0; i < v.size(); i++)
	{
		if(!(v[i])) return false;
	}
	return true;
}

inline bool all_true(const bool v)
{
	return v;
}
#endif // all_true

// all_false
#if (1)
inline bool all_false(const std::vector<bool> v)
{
	for(unsigned int i=0; i < v.size(); i++)
	{
		if(v[i]) return false;
	}
	return true;
}

inline bool all_false(const bool v)
{
	return (!v);
}
#endif // all_false

// not_all_true
#if (1)
inline bool not_all_true(const std::vector<bool> v)
{
	return !all_true(v);
}

inline bool not_all_true(const bool v)
{
	return !v;
}
#endif // not_all_true

// not_all_false
#if (1)
inline bool not_all_false(const std::vector<bool> v)
{
	return !all_false(v);
}

inline bool not_all_false(const bool v)
{
	return v;
}
#endif // all_false

// some_true_some_false
#if (1)
inline bool some_true_some_false(const std::vector<bool> v)
{
	return ( (!all_true(v)) && (!all_false(v)) );
}

inline bool some_true_some_false(const bool v)
{
	return false;
}
#endif // some_true_some_false

} // namespace brgastro



#endif /* SUMMARY_FUNCTIONS_HPP_ */
