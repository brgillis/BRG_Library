/**********************************************************************\
 @file mag_correlation_function_estimator.cpp
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

#include "IceBRG_main/common.hpp"

#include "IceBRG_main/Eigen.hpp"

#include "IceBRG_main/math/misc_math.hpp"
#include "IceBRG_main/math/safe_math.hpp"
#include "IceBRG_main/utility.hpp"
#include "IceBRG_main/units/units.hpp"
#include "IceBRG_main/vector/limit_vector.hpp"

#include "SHE_GST_IceBRG_lensing/magnification/mag_correlation_function_estimator.hpp"
#include "SHE_GST_IceBRG_lensing/magnification/magnification_alpha.hpp"


namespace IceBRG {

bool mag_correlation_function_estimator::_set_up() const
{
	return !( (_D1_pos_list_.empty()) || (_D2_pos_list_.empty())
			|| (_R1_pos_list_.empty()) || (_R2_pos_list_.empty()));
}


// Calculation functions
Eigen::ArrayXd mag_correlation_function_estimator::calculate_weighted(
		const std::function<flt_t(angle_type)> & weight_function) const
{
	if(!_set_up())
		throw std::logic_error("Cannot calculate correlation function without data set up.\n");

	// We'll calculate four arrays here, for the combinations of data and random lists
	Eigen::ArrayXd D1D2_counts = Eigen::ArrayXd::Zero(_r_bin_limits_.num_bins());
	Eigen::ArrayXd D1R2_counts = Eigen::ArrayXd::Zero(_r_bin_limits_.num_bins());
	Eigen::ArrayXd R1D2_counts = Eigen::ArrayXd::Zero(_r_bin_limits_.num_bins());
	Eigen::ArrayXd R1R2_counts = Eigen::ArrayXd::Zero(_r_bin_limits_.num_bins());

	int_t D1D2_pairs = 0;
	int_t D1R2_pairs = 0;
	int_t R1D2_pairs = 0;
	int_t R1R2_pairs = 0;

	const angle_type & max_r = _r_bin_limits_.max();

	// Set up a function to add to the correct bin of an array
	auto increment_bin = [&] (Eigen::ArrayXd & array, int_t & pair_counter,
			const position & p1,
			const position & p2,
			const flt_t & scale )
	{
		flt_t dz = std::get<2>(p2)-std::get<2>(p1);
		if(dz<_z_buffer_) return;

		++pair_counter;

		angle_type dx = std::get<0>(p1)-std::get<0>(p2);
		if(abs(dx)>max_r) return;
		angle_type dy = std::get<1>(p1)-std::get<1>(p2);
		if(abs(dy)>max_r) return;

		angle_type r = IceBRG::quad_add(dx,dy);

		if(_r_bin_limits_.outside_limits(r)) return;

		angle_type theta = atan2(dy,dx);

		size_t i = _r_bin_limits_.get_bin_index(r);

		array[i] += scale*weight_function(theta);

		return;
	};

	flt_t am1m = 0;
	int_t p2_count = 0;

	// Loop over all combinations
	for(const auto & p1 : _D1_pos_list_)
	{
		for(const auto & p2 : _D2_pos_list_)
		{
			flt_t am1 = IceBRG::magnification_alpha(std::get<3>(p2),std::get<2>(p1)+_z_buffer_)-1;

			increment_bin(D1D2_counts,D1D2_pairs,p1,p2,am1);
		}
		for(const auto & p2 : _R2_pos_list_)
		{
			increment_bin(D1R2_counts,D1R2_pairs,p1,p2,1.);
		}
	}
	for(const auto & p1 : _R1_pos_list_)
	{
		for(const auto & p2 : _D2_pos_list_)
		{
			flt_t am1 = IceBRG::magnification_alpha(std::get<3>(p2),std::get<2>(p1)+_z_buffer_)-1;
			am1m += am1;
			++p2_count;

			increment_bin(R1D2_counts,R1D2_pairs,p1,p2,am1);
		}
		for(const auto & p2 : _R2_pos_list_)
		{
			increment_bin(R1R2_counts,R1R2_pairs,p1,p2,1.);
		}
	}

	// Normalize by total number of pairs

	am1m /= p2_count;

	D1D2_counts /= D1D2_pairs;
	D1R2_counts /= D1R2_pairs;
	R1D2_counts /= R1D2_pairs;
	R1R2_counts /= R1R2_pairs;

	// And calculate the correlation function
//	Eigen::ArrayXd result = D1D2_counts/IceBRG::safe_d(D1R2_counts)-1;
	Eigen::ArrayXd result = (D1D2_counts-am1m*D1R2_counts-R1D2_counts) / IceBRG::safe_d(R1R2_counts) + am1m;
	return result;
}
Eigen::ArrayXd mag_correlation_function_estimator::calculate() const
{
	if(!_set_up())
		throw std::logic_error("Cannot calculate correlation function without data set up.\n");

	if(_unweighted_cached_value_.size()>0) return _unweighted_cached_value_;

	// We'll calculate four arrays here, for the combinations of data and random lists
	Eigen::ArrayXd D1D2_unnorm_counts = Eigen::ArrayXd::Zero(_r_bin_limits_.num_bins());
	Eigen::ArrayXd D1D2_counts = Eigen::ArrayXd::Zero(_r_bin_limits_.num_bins());
	Eigen::ArrayXd D1R2_unnorm_counts = Eigen::ArrayXd::Zero(_r_bin_limits_.num_bins());
	Eigen::ArrayXd D1R2_counts = Eigen::ArrayXd::Zero(_r_bin_limits_.num_bins());
	Eigen::ArrayXd R1D2_unnorm_counts = Eigen::ArrayXd::Zero(_r_bin_limits_.num_bins());
	Eigen::ArrayXd R1D2_counts = Eigen::ArrayXd::Zero(_r_bin_limits_.num_bins());
	Eigen::ArrayXd R1R2_unnorm_counts = Eigen::ArrayXd::Zero(_r_bin_limits_.num_bins());
	Eigen::ArrayXd R1R2_counts = Eigen::ArrayXd::Zero(_r_bin_limits_.num_bins());

	long_int_t D1D2_unnorm_pairs = 0;
	long_int_t D1D2_pairs = 0;
	long_int_t D1R2_unnorm_pairs = 0;
	long_int_t D1R2_pairs = 0;
	long_int_t R1D2_unnorm_pairs = 0;
	long_int_t R1D2_pairs = 0;
	long_int_t R1R2_unnorm_pairs = 0;
	long_int_t R1R2_pairs = 0;

	const angle_type & max_r = _r_bin_limits_.max();

	// Set up a function to add to the correct bin of an array
	auto increment_bin = [&] (Eigen::ArrayXd & array, long_int_t & pair_counter,
			const position & p1,
			const position & p2,
			const flt_t & scale )
	{
		++pair_counter;

		angle_type dx = std::get<0>(p1)-std::get<0>(p2);
		if(abs(dx)>max_r) return;
		angle_type dy = std::get<1>(p1)-std::get<1>(p2);
		if(abs(dy)>max_r) return;

		angle_type r = IceBRG::quad_add(dx,dy);

		if(_r_bin_limits_.outside_limits(r)) return;

		size_t i = _r_bin_limits_.get_bin_index(r);

		array[i] += scale;

		return;
	};

	flt_t am1m = 0;
	int_t p2_count = 0;

	// Loop over all combinations
	for(const auto & p1 : _D1_pos_list_)
	{

		for(const auto & p2 : _D2_pos_list_)
		{
			flt_t dz = std::get<2>(p2)-std::get<2>(p1);
			if(dz<_z_buffer_) continue;

			flt_t am1 = IceBRG::magnification_alpha(std::get<3>(p2),std::get<2>(p1)+_z_buffer_)-1;
			am1m += am1;
			++p2_count;

			increment_bin(D1D2_unnorm_counts,D1D2_unnorm_pairs,p1,p2,1.);
			increment_bin(D1D2_counts,D1D2_pairs,p1,p2,am1);
		}
		for(const auto & p2 : _R2_pos_list_)
		{
			increment_bin(D1R2_unnorm_counts,D1R2_unnorm_pairs,p1,p2,1.);
			increment_bin(D1R2_counts,D1R2_pairs,p1,p2,1.);
		}
	}
	for(const auto & p1 : _R1_pos_list_)
	{
		for(const auto & p2 : _D2_pos_list_)
		{
			flt_t dz = std::get<2>(p2)-std::get<2>(p1);
			if(dz<_z_buffer_) continue;

			flt_t am1 = IceBRG::magnification_alpha(std::get<3>(p2),std::get<2>(p1)+_z_buffer_)-1;

			increment_bin(R1D2_unnorm_counts,R1D2_unnorm_pairs,p1,p2,1.);
			increment_bin(R1D2_counts,R1D2_pairs,p1,p2,am1);
		}
		for(const auto & p2 : _R2_pos_list_)
		{
			increment_bin(R1R2_unnorm_counts,R1R2_unnorm_pairs,p1,p2,1.);
			increment_bin(R1R2_counts,R1R2_pairs,p1,p2,1);
		}
	}

	am1m /= p2_count;

	// Normalize by total number of pairs

	D1D2_counts /= D1D2_pairs;
	D1R2_counts /= D1R2_pairs;
	R1D2_counts /= R1D2_pairs;
	R1R2_counts /= R1R2_pairs;

	// And calculate the correlation function

	_unweighted_cached_value_ = (D1D2_counts-am1m*D1R2_counts-R1D2_counts) / IceBRG::safe_d(R1R2_counts) + am1m; // LS

	auto percent_err_D1D2_counts = 1. /  D1D2_unnorm_counts.sqrt();
	auto percent_err_D1R2_counts = 1. /  D1R2_unnorm_counts.sqrt();
	auto percent_err_R1D2_counts = 1. /  R1D2_unnorm_counts.sqrt();
	auto percent_err_R1R2_counts = 1. /  R1R2_unnorm_counts.sqrt();
	auto percent_err_am1m = 1./std::sqrt(p2_count);

	_unweighted_cached_error_ = D1D2_counts/IceBRG::safe_d(R1R2_counts) * (percent_err_D1D2_counts.square()+percent_err_R1R2_counts.square()).sqrt()
		+ std::abs(am1m)*D1R2_counts/IceBRG::safe_d(R1R2_counts) * (percent_err_D1R2_counts.square()+percent_err_R1R2_counts.square()+
			square(percent_err_am1m)).sqrt()
		+ R1D2_counts/IceBRG::safe_d(R1R2_counts) * (percent_err_R1D2_counts.square()+percent_err_R1R2_counts.square()).sqrt()
		+ am1m*percent_err_am1m;

//	_unweighted_cached_error_ = D1D2_counts/IceBRG::safe_d<Eigen::ArrayXd>(D1D2_unnorm_counts.sqrt()*R1R2_counts)
//		+ std::abs(am1m)*D1R2_counts/IceBRG::safe_d<Eigen::ArrayXd>(D1R2_unnorm_counts.sqrt()*R1R2_counts)
//		+ R1D2_counts/IceBRG::safe_d<Eigen::ArrayXd>(R1D2_unnorm_counts.sqrt()*R1R2_counts);

	return _unweighted_cached_value_;
}

Eigen::ArrayXd mag_correlation_function_estimator::weights() const
{
	if(!_set_up())
		throw std::logic_error("Cannot get weights without data set up.\n");

	if(_unweighted_cached_error_.size()==0) calculate();

	//return 1./_unweighted_cached_error_.square();

	return Eigen::ArrayXd::Constant(_unweighted_cached_error_.size(),1);
}

Eigen::ArrayXd mag_correlation_function_estimator::errors() const
{
	if(!_set_up())
		throw std::logic_error("Cannot get weights without data set up.\n");

	if(_unweighted_cached_error_.size()==0) calculate();

	return _unweighted_cached_error_;
}

Eigen::ArrayXd mag_correlation_function_estimator::calculate_dipole(const flt_t & offset) const
{
	auto weight_function = [&offset] (const angle_type & theta)
	{
		return 1.+sin(theta + 2*pi*rad*offset);
	};
	return calculate_weighted(weight_function)-calculate();
}
Eigen::ArrayXd mag_correlation_function_estimator::calculate_quadrupole(const flt_t & offset) const
{
	auto weight_function = [&offset] (const angle_type & theta)
	{
		return 1.+sin((theta + 2*pi*rad*offset)/2.);
	};
	return calculate_weighted(weight_function)-calculate();
}
Eigen::ArrayXd mag_correlation_function_estimator::calculate_octopole(const flt_t & offset) const
{
	auto weight_function = [&offset] (const angle_type & theta)
	{
		return 1.+sin((theta + 2*pi*rad*offset)/4.);
	};
	return calculate_weighted(weight_function)-calculate();
}

} // end namespace IceBRG
