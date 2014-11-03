/**********************************************************************\
 @file limit_vector_operations.hpp
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

#ifndef _BRG_LIMIT_VECTOR_OPERATIONS_HPP_INCLUDED_
#define _BRG_LIMIT_VECTOR_OPERATIONS_HPP_INCLUDED_

#include <cassert>
#include <stdexcept>
#include <vector>

#include "brg/global.h"

#include "brg/math/misc_math.hpp"
#include "brg/vector/summary_functions.hpp"

namespace brgastro {

template<typename T>
std::vector<T> make_limit_vector(const T & min, const T & max, const T & step)
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
std::vector<T> make_log_limit_vector(const T & min, const T & max, const size_t & num_bins)
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
		result.push_back( min * brgastro::ipow(log_step,limit_num) );
	}

	result.push_back(max);


	return result;
}

template<typename T>
bool above_limits(const T & val, const std::vector<T> & vec)
{
	assert(is_monotonically_increasing(vec));
	return val>vec.back();
}

template<typename T>
bool checked_above_limits(const T & val, const std::vector<T> & vec)
{
	if(!is_monotonically_increasing(vec))
		throw std::logic_error("Invalid limit vector passed to above_limits.");
	return above_limits(val,vec);
}

template<typename T>
bool under_limits(const T & val, const std::vector<T> & vec)
{
	assert(is_monotonically_increasing(vec));
	return val<vec.front();
}

template<typename T>
bool checked_under_limits(const T & val, const std::vector<T> & vec)
{
	if(!is_monotonically_increasing(vec))
		throw std::logic_error("Invalid limit vector passed to under_limits.");
	return under_limits(val,vec);
}

template<typename T>
bool outside_limits(const T & val, const std::vector<T> & vec)
{
	return (above_limits(val,vec) or under_limits(val,vec));
}

template<typename T>
bool checked_outside_limits(const T & val, const std::vector<T> & vec)
{
	return (checked_above_limits(val,vec) or checked_under_limits(val,vec));
}

template<typename T>
bool inside_limits(const T & val, const std::vector<T> & vec)
{
	return !outside_limits(val,vec);
}

template<typename T>
bool checked_inside_limits(const T & val, const std::vector<T> & vec)
{
	return !checked_outside_limits(val,vec);
}

template<typename T>
size_t checked_get_bin_index(const T & val, const std::vector<T> & vec)
{
	if(!is_monotonically_increasing(vec))
		throw std::logic_error("Invalid limit vector passed to checked_get_bin_index.");
	if(outside_limits(val,vec))
		throw std::runtime_error("Value outside limits of vector in checked_get_bin_index.");
	return get_bin_index(val,vec);
}

template<typename T>
size_t get_bin_index(const T & val, const std::vector<T> & vec)
{
	assert(is_monotonically_increasing(vec));
	for(size_t i=1; i<vec.size(); ++i)
	{
		if(vec[i]>=val) return i-1;
	}
	if(above_limits(val,vec)) return vec.size()-2;
	return 0;
}

template<typename T1, typename T2>
T2 interpolate_bins(const T2 & val, const std::vector<T1> & lim_vec, const std::vector<T2> & val_vec)
{
	assert(val_vec.size()>1);
	assert(lim_vec.size()==val_vec.size()+1);
	assert(is_monotonically_increasing(lim_vec));
	assert(lim_vec.size()>2);

	size_t bin_i=lim_vec.size()-3;

	for(size_t i=0; i<lim_vec.size()-2; ++i)
	{
		if((lim_vec[i]+lim_vec[i+1])/2>=val)
		{
			if(i==0)
				bin_i=0;
			else
				bin_i=i-1;
			break;
		}
	}

	T1 xlo = (lim_vec[bin_i]+lim_vec[bin_i+1])/2;
	T1 xhi = (lim_vec[bin_i+1]+lim_vec[bin_i+2])/2;
	const T2 & ylo = val_vec[bin_i];
	const T2 & yhi = val_vec[bin_i+1];

	return ylo + (yhi-ylo)/(xhi-xlo) * (val-xlo);
}

template<typename T1>
std::vector<T1> get_bin_mids_from_limits(std::vector<T1> vec)
{
	assert(is_monotonically_increasing(vec));

	for(size_t i=0; i<vec.size()-1; ++i)
	{
		vec[i] += (vec[i+1]-vec[i])/2;
	}

	vec.pop_back();
	return vec;
}

template<typename T1>
std::vector<T1> get_bin_limits_from_mids(std::vector<T1> vec)
{
	assert(is_monotonically_increasing(vec));

	size_t i;

	for(i=0; i<vec.size()-1; ++i)
	{
		vec[i] -= (vec[i+1]-vec[i])/2;
	}

	// Special handling for the final element
	T1 d_last = (vec[i]-vec[i-1])/3;
	vec.push_back(vec[i]+d_last);
	vec[i]-=d_last;

	return vec;
}


} // namespace brgastro



#endif // _BRG_LIMIT_VECTOR_OPERATIONS_HPP_INCLUDED_
