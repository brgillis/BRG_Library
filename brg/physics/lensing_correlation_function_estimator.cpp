/**********************************************************************\
 @file lensing_lensing_correlation_function_estimator.cpp
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
#include <tuple>
#include <vector>

#include "brg/global.h"

#include "brg/math/misc_math.hpp"
#include "brg/math/safe_math.hpp"
#include "brg/vector/limit_vector.hpp"

#include "lensing_correlation_function_estimator.h"

namespace brgastro {

bool lensing_correlation_function_estimator::_set_up()
{
	return !( (_D1_pos_list_.empty()) || (_D2_pos_list_.empty())
			|| (_R1_pos_list_.empty()) || (_R2_pos_list_.empty()));
}


// Calculation functions
std::valarray<double> lensing_correlation_function_estimator::calculate_weighted(
		const std::function<double(double)> & weight_function)
{
	if(!_set_up())
		throw std::logic_error("Cannot calculate correlation function without data set up.\n");

	// We'll calculate four valarrays here, for the combinations of data and random lists
	std::valarray<double> D1D2_counts(_r_bin_limits_.num_bins());
	std::valarray<double> D1R2_counts(_r_bin_limits_.num_bins());
	std::valarray<double> R1D2_counts(_r_bin_limits_.num_bins());
	std::valarray<double> R1R2_counts(_r_bin_limits_.num_bins());

	unsigned D1D2_pairs = 0;
	unsigned D1R2_pairs = 0;
	unsigned R1D2_pairs = 0;
	unsigned R1R2_pairs = 0;

	const double & max_r = _r_bin_limits_.max();

	// Set up a function to add to the correct bin of a valarray
	auto increment_bin = [&] (std::valarray<double> & array, unsigned & pair_counter,
			const std::tuple<double,double,double> & p1,
			const std::tuple<double,double,double> & p2)
	{
		double dz = std::get<2>(p2)-std::get<2>(p1);
		if(dz<_z_buffer_) return;

		++pair_counter;

		double dx = std::get<0>(p1)-std::get<0>(p2);
		if(std::fabs(dx)>max_r) return;
		double dy = std::get<1>(p1)-std::get<1>(p2);
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
			increment_bin(D1D2_counts,D1D2_pairs,p1,p2);
		}
		for(const auto & p2 : _R2_pos_list_)
		{
			increment_bin(D1R2_counts,D1R2_pairs,p1,p2);
		}
	}
	for(const auto & p1 : _R1_pos_list_)
	{
		for(const auto & p2 : _D2_pos_list_)
		{
			increment_bin(R1D2_counts,R1D2_pairs,p1,p2);
		}
		for(const auto & p2 : _R2_pos_list_)
		{
			increment_bin(R1R2_counts,R1R2_pairs,p1,p2);
		}
	}

	// Normalize by total number of pairs
	D1D2_counts /= D1D2_pairs;
	D1R2_counts /= D1R2_pairs;
	R1D2_counts /= R1D2_pairs;
	R1R2_counts /= R1R2_pairs;

	// And calculate the correlation function
	std::valarray<double> result = (D1D2_counts*R1R2_counts)/
			brgastro::safe_d(static_cast< std::valarray<double> >(D1R2_counts*R1D2_counts))-1;
	return result;
}
std::valarray<double> lensing_correlation_function_estimator::calculate()
{
	if(!_set_up())
		throw std::logic_error("Cannot calculate correlation function without data set up.\n");

	if(_unweighted_cached_value_.size()>0) return _unweighted_cached_value_;

	// We'll calculate four valarrays here, for the combinations of data and random lists
	std::valarray<double> D1D2_counts(_r_bin_limits_.num_bins());
	std::valarray<double> D1R2_counts(_r_bin_limits_.num_bins());
	std::valarray<double> R1D2_counts(_r_bin_limits_.num_bins());
	std::valarray<double> R1R2_counts(_r_bin_limits_.num_bins());

	unsigned D1D2_pairs = 0;
	unsigned D1R2_pairs = 0;
	unsigned R1D2_pairs = 0;
	unsigned R1R2_pairs = 0;

	const double & max_r = _r_bin_limits_.max();

	// Set up a function to add to the correct bin of a valarray
	auto increment_bin = [&] (std::valarray<double> & array, unsigned & pair_counter,
			const std::tuple<double,double,double> & p1,
			const std::tuple<double,double,double> & p2)
	{
		double dz = std::get<2>(p2)-std::get<2>(p1);
		if(dz<_z_buffer_) return;

		++pair_counter;

		double dx = std::get<0>(p1)-std::get<0>(p2);
		if(std::fabs(dx)>max_r) return;
		double dy = std::get<1>(p1)-std::get<1>(p2);
		if(std::fabs(dy)>max_r) return;

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
			increment_bin(D1D2_counts,D1D2_pairs,p1,p2);
		}
		for(const auto & p2 : _R2_pos_list_)
		{
			increment_bin(D1R2_counts,D1R2_pairs,p1,p2);
		}
	}
	for(const auto & p1 : _R1_pos_list_)
	{
		for(const auto & p2 : _D2_pos_list_)
		{
			increment_bin(R1D2_counts,R1D2_pairs,p1,p2);
		}
		for(const auto & p2 : _R2_pos_list_)
		{
			increment_bin(R1R2_counts,R1R2_pairs,p1,p2);
		}
	}

	// Normalize by total number of pairs
	D1D2_counts /= D1D2_pairs;
	D1R2_counts /= D1R2_pairs;
	R1D2_counts /= R1D2_pairs;
	R1R2_counts /= R1R2_pairs;

	// And calculate the correlation function
	_unweighted_cached_value_.resize(_r_bin_limits_.num_bins());
	_unweighted_cached_value_ = (D1D2_counts*R1R2_counts)/
			brgastro::safe_d(static_cast< std::valarray<double> >(D1R2_counts*R1D2_counts)) - 1;
	return _unweighted_cached_value_;
}

std::valarray<double> lensing_correlation_function_estimator::calculate_dipole(const double & offset)
{
	auto weight_function = [&offset] (const double & theta)
	{
		return 1+std::sin(theta + 2*pi*offset);
	};
	return calculate_weighted(weight_function)-calculate();
}
std::valarray<double> lensing_correlation_function_estimator::calculate_quadrupole(const double & offset)
{
	auto weight_function = [&offset] (const double & theta)
	{
		return 1+std::sin((theta + 2*pi*offset)/2);
	};
	return calculate_weighted(weight_function)-calculate();
}
std::valarray<double> lensing_correlation_function_estimator::calculate_octopole(const double & offset)
{
	auto weight_function = [&offset] (const double & theta)
	{
		return 1+std::sin((theta + 2*pi*offset)/4);
	};
	return calculate_weighted(weight_function)-calculate();
}

} // end namespace brgastro
