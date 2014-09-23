/**********************************************************************\
 @file limit_vector.hpp
 ------------------

 TODO <Insert file description here>

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

#ifndef _BRG_LIMIT_VECTOR_HPP_INCLUDED_
#define _BRG_LIMIT_VECTOR_HPP_INCLUDED_

#include <cassert>
#include <stdexcept>
#include <vector>

#include "brg/global.h"

#include "brg/math/misc_math.hpp"

namespace brgastro {

template<typename T>
std::vector<T> make_limit_vector(T min, T max, T step)
{
	assert(max>min);
	assert(step>0);

	std::vector<T> result(1,min);

	if(!isinf(step))
	{
		if( isinf(min) || isinf(max) )
			throw std::logic_error("Cannot generate limit vector with finite step and infinite limits.");

		size_t num_bins = (size_t)std::ceil((max-min)/step);

		for(size_t limit_num=1; limit_num<num_bins; ++limit_num)
		{
			// Recalculate at each step to minimize round-off error
			result.push_back( min + step*limit_num );
		}
	}
	result.push_back(max);

	return result;
}

template<typename T>
std::vector<T> make_log_limit_vector(T min, T max, size_t num_bins)
{
	assert(max>min);
	assert(min>0);
	assert(num_bins>0);

	std::vector<T> result(1,min);

	double log_step = std::pow((double)max/(double)min,1./(num_bins));

	if( isinf(max) )
	{
		if(num_bins>1)
			throw std::logic_error("Cannot generate limit vector with finite step and infinite limits.");
		result.push_back(max);
		return result;
	}

	for(size_t limit_num=1; limit_num<num_bins; ++limit_num)
	{
		// Recalculate at each step to minimize round-off error
		result.push_back( min * std::pow(log_step,limit_num) );
	}

	result.push_back(max);


	return result;
}

template<typename T>
size_t checked_get_bin_index(T val, std::vector<T> vec)
{
	if(!is_monotonically_increasing(vec))
		throw std::logic_error("Invalid limit vector passed to checked_get_bin_index.");
	return get_bin_index(val,vec);
}

template<typename T>
size_t get_bin_index(T val, std::vector<T> vec)
{
	if((val<vec.front())||(val>vec.back()))
			throw std::runtime_error("Value outside limits of vector in get_bin_index.");
	for(size_t i=1; i<vec.size(); ++i)
	{
		if(vec[i]>=val) return i-1;
	}
	throw std::logic_error("Invalid limit vector passed to get_bin_index.");
}

} // namespace brgastro



#endif // _BRG_LIMIT_VECTOR_HPP_INCLUDED_
