/**********************************************************************\
 @file correlation_function_estimator.cpp
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

#include <functional>
#include <utility>
#include <vector>

#include "brg/global.h"

#include "brg/math/misc_math.hpp"
#include "brg/math/safe_math.hpp"
#include "brg/vector/limit_vector.hpp"

#include "correlation_function_estimator.h"

namespace brgastro {

bool correlation_function_estimator::_set_up()
{
	return !( (_D1_pos_list_.empty()) || (_D2_pos_list_.empty())
			|| (_R1_pos_list_.empty()) || (_R2_pos_list_.empty()));
}


// Calculation functions
std::valarray<double> correlation_function_estimator::calculate_weighted(
		const std::function<double(double)> & weight_function)
{
	if(!_set_up())
		throw std::logic_error("Cannot calculate correlation function without data set up.\n");

	// We'll calculate four valarrays here, for the combinations of data and random lists
	std::valarray<double> D1D2_counts(_r_bin_limits_.num_bins());
	std::valarray<double> D1R2_counts(_r_bin_limits_.num_bins());
	std::valarray<double> R1D2_counts(_r_bin_limits_.num_bins());
	std::valarray<double> R1R2_counts(_r_bin_limits_.num_bins());

	const double & max_r = _r_bin_limits_.max();

	// Set up a function to add to the correct bin of a valarray
	auto increment_bin = [&] (std::valarray<double> & array, const std::pair<double,double> & p1,
			const std::pair<double,double> & p2)
	{
		double dx = p1.first-p2.first;
		if(std::fabs(dx)>max_r) return;
		double dy = p1.second-p2.second;
		if(std::fabs(dy)>max_r) return;

		double r = brgastro::quad_add(dx,dy);

		if(_r_bin_limits_.outside_limits(r)) return;

		double theta = std::atan2(dy,dx);

		size_t i = _r_bin_limits_.get_bin_index(r);

		array[i] += weight_function(theta);

		return;
	};

	// Loop over all combinations
	for(const auto & p1 : _D1_pos_list_)
	{
		for(const auto & p2 : _D2_pos_list_)
		{
			increment_bin(D1D2_counts,p1,p2);
		}
		for(const auto & p2 : _R2_pos_list_)
		{
			increment_bin(D1R2_counts,p1,p2);
		}
	}
	for(const auto & p1 : _R1_pos_list_)
	{
		for(const auto & p2 : _D2_pos_list_)
		{
			increment_bin(R1D2_counts,p1,p2);
		}
		for(const auto & p2 : _R2_pos_list_)
		{
			increment_bin(R1R2_counts,p1,p2);
		}
	}

	// And calculate the correlation function
	std::valarray<double> result = (D1D2_counts*R1R2_counts)/
			brgastro::safe_d(static_cast< std::valarray<double> >(D1R2_counts*R1D2_counts))-1;
	return result;
}
std::valarray<double> correlation_function_estimator::calculate()
{
	if(!_set_up())
		throw std::logic_error("Cannot calculate correlation function without data set up.\n");

	if(_unweighted_cached_value_.size()>0) return _unweighted_cached_value_;

	// We'll calculate four valarrays here, for the combinations of data and random lists
	std::valarray<double> D1D2_counts(_r_bin_limits_.num_bins());
	std::valarray<double> D1R2_counts(_r_bin_limits_.num_bins());
	std::valarray<double> R1D2_counts(_r_bin_limits_.num_bins());
	std::valarray<double> R1R2_counts(_r_bin_limits_.num_bins());

	const double & max_r = _r_bin_limits_.max();

	// Set up a function to add to the correct bin of a valarray
	auto increment_bin = [&] (std::valarray<double> & array, const std::pair<double,double> & p1,
			const std::pair<double,double> & p2)
	{
		double dx = std::fabs(p1.first-p2.first);
		if(dx>max_r) return;
		double dy = std::fabs(p1.second-p2.second);
		if(dy>max_r) return;

		double r = brgastro::quad_add(dx,dy);

		if(_r_bin_limits_.outside_limits(r)) return;

		size_t i = _r_bin_limits_.get_bin_index(r);

		array[i] += 1;

		return;
	};

	// Loop over all combinations
	for(const auto & p1 : _D1_pos_list_)
	{
		for(const auto & p2 : _D2_pos_list_)
		{
			increment_bin(D1D2_counts,p1,p2);
		}
		for(const auto & p2 : _R2_pos_list_)
		{
			increment_bin(D1R2_counts,p1,p2);
		}
	}
	for(const auto & p1 : _R1_pos_list_)
	{
		for(const auto & p2 : _D2_pos_list_)
		{
			increment_bin(R1D2_counts,p1,p2);
		}
		for(const auto & p2 : _R2_pos_list_)
		{
			increment_bin(R1R2_counts,p1,p2);
		}
	}

	// And calculate the correlation function
	_unweighted_cached_value_.resize(_r_bin_limits_.num_bins());
	_unweighted_cached_value_ = (D1D2_counts*R1R2_counts)/
			brgastro::safe_d(static_cast< std::valarray<double> >(D1R2_counts*R1D2_counts)) - 1;
	return _unweighted_cached_value_;
}

std::valarray<double> correlation_function_estimator::calculate_dipole(const double & offset)
{
	auto weight_function = [&offset] (const double & theta)
	{
		return 1+std::sin(theta + 2*pi*offset);
	};
	return calculate_weighted(weight_function)-calculate();
}
std::valarray<double> correlation_function_estimator::calculate_quadrupole(const double & offset)
{
	auto weight_function = [&offset] (const double & theta)
	{
		return 1+std::sin((theta + 2*pi*offset)/2);
	};
	return calculate_weighted(weight_function)-calculate();
}
std::valarray<double> correlation_function_estimator::calculate_octopole(const double & offset)
{
	auto weight_function = [&offset] (const double & theta)
	{
		return 1+std::sin((theta + 2*pi*offset)/4);
	};
	return calculate_weighted(weight_function)-calculate();
}

} // end namespace brgastro
