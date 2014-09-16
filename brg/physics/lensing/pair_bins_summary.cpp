/**********************************************************************\
 @file pair_binner.cpp
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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/iterator/zip_iterator.hpp>

#include "brg/global.h"

#include "brg/file_access/ascii_table.h"
#include "brg/physics/units/unit_obj.h"
#include "brg/utility.hpp"
#include "brg/vector/limit_vector.hpp"
#include "brg/vector/summary_functions.hpp"

#include "pair_binner.h"

namespace brgastro {

	// Private methods
#if(1)

void pair_binner::_check_limits()
{
	// Now check they're all monotonically increasing
	// Note that this function returns false if they're too small as well
	if(!is_monotonically_increasing(_R_bin_limits_))
	{
		_valid_limits_ = false;
		return;
	}
	if(!is_monotonically_increasing(_m_bin_limits_))
	{
		_valid_limits_ = false;
		return;
	}
	if(!is_monotonically_increasing(_z_bin_limits_))
	{
		_valid_limits_ = false;
		return;
	}
	if(!is_monotonically_increasing(_mag_bin_limits_))
	{
		_valid_limits_ = false;
		return;
	}

	_valid_limits_ = true;
}

void pair_binner::_sort() const
{
	if(_sorted_) return;

	if(!_valid_limits_) throw std::logic_error("Pairs can't be sorted without valid bin limits.");

	// Check if the pair_bins vector has been generated and is the proper size
	if(_pair_bins_.empty())
	{
		make_array4d(_pair_bins_,_R_bin_limits_.size()-1,_m_bin_limits_.size()-1,
				_z_bin_limits_.size()-1,_mag_bin_limits_.size()-1);
		for(size_t R_i=0; R_i<_R_bin_limits_.size()-1;++R_i)
		{
			for(size_t m_i=0; m_i<_m_bin_limits_.size()-1;++m_i)
			{
				for(size_t z_i=0; z_i<_z_bin_limits_.size()-1;++z_i)
				{
					for(size_t mag_i=0; mag_i<_mag_bin_limits_.size()-1;++mag_i)
					{
						_pair_bins_[R_i][m_i][z_i][mag_i] = pair_bin(
								_R_bin_limits_[R_i],_R_bin_limits_[R_i+1],
								_m_bin_limits_[m_i],_m_bin_limits_[m_i+1],
								_z_bin_limits_[z_i],_z_bin_limits_[z_i+1],
								_mag_bin_limits_[mag_i],_mag_bin_limits_[mag_i+1]);
					}
				}
			}
		}
	}
	else
	{
		if(_pair_bins_.size() != _R_bin_limits_.size()-1)
		{
			_resort();
			return;
		}
		if(_pair_bins_[0].size() != _m_bin_limits_.size()-1)
		{
			_resort();
			return;
		}
		if(_pair_bins_[0][0].size() != _z_bin_limits_.size()-1)
		{
			_resort();
			return;
		}
		if(_pair_bins_[0][0][0].size() != _mag_bin_limits_.size()-1)
		{
			_resort();
			return;
		}
	}

	// Now loop through and sort each pair into the proper bin
	// Note that we don't zero the sorting index here. This is intentional, so
	// we don't have to resort previously sorted pairs when new ones are added.
	for(; _sorting_index_<_pairs_.size(); ++_sorting_index_)
	{
		// Check bounds first
		BRG_DISTANCE R_proj = _pairs_[_sorting_index_].R_proj();
		if((R_proj < _R_bin_limits_.front()) || (R_proj > _R_bin_limits_.back()))
			continue;
		BRG_MASS m = _pairs_[_sorting_index_].m_lens();
		if((m < _m_bin_limits_.front()) || (m > _m_bin_limits_.back()))
			continue;
		double z = _pairs_[_sorting_index_].z_lens();
		if((z < _z_bin_limits_.front()) || (z > _z_bin_limits_.back()))
			continue;
		double mag = _pairs_[_sorting_index_].mag_lens();
		if((mag < _mag_bin_limits_.front()) || (mag > _mag_bin_limits_.back()))
			continue;

		// Find specific position by using a zip iterator to go over limits vectors
		// and the pair bins vector at the same time
		auto R_zip_begin = boost::make_zip_iterator(boost::make_tuple(
				_R_bin_limits_.begin()+1,_pair_bins_.begin()));
		auto R_zip_end = boost::make_zip_iterator(boost::make_tuple(
				_R_bin_limits_.end(),_pair_bins_.end()));

		auto R_func = &t1first_lt_v2<decltype(*R_zip_begin),BRG_DISTANCE>;

		auto & pair_R_vec = std::lower_bound(R_zip_begin,R_zip_end,R_proj,R_func)->get<1>();



		auto m_zip_begin = boost::make_zip_iterator(boost::make_tuple(
				_m_bin_limits_.begin()+1,pair_R_vec.begin()));
		auto m_zip_end = boost::make_zip_iterator(boost::make_tuple(
				_m_bin_limits_.end(),pair_R_vec.end()));

		auto m_func = &t1first_lt_v2<decltype(*m_zip_begin),BRG_MASS>;

		auto & pair_Rm_vec = std::lower_bound(m_zip_begin,m_zip_end,m,m_func)->get<1>();



		auto z_zip_begin = boost::make_zip_iterator(boost::make_tuple(
				_z_bin_limits_.begin()+1,pair_Rm_vec.begin()));
		auto z_zip_end = boost::make_zip_iterator(boost::make_tuple(
				_z_bin_limits_.end(),pair_Rm_vec.end()));

		auto z_func = &t1first_lt_v2<decltype(*z_zip_begin),double>;

		auto & pair_Rmz_vec = std::lower_bound(z_zip_begin,z_zip_end,z,z_func)->get<1>();



		auto mag_zip_begin = boost::make_zip_iterator(boost::make_tuple(
				_mag_bin_limits_.begin()+1,pair_Rmz_vec.begin()));
		auto mag_zip_end = boost::make_zip_iterator(boost::make_tuple(
				_mag_bin_limits_.end(),pair_Rmz_vec.end()));

		auto mag_func = &t1first_lt_v2<decltype(*mag_zip_begin),double>;

		auto & pair_Rmzmag = std::lower_bound(mag_zip_begin,mag_zip_end,mag,mag_func)->get<1>();

		// At this point, pair_Rmzmag is a reference to the pair bin we want to add this pair to
		pair_Rmzmag.add_pair(_pairs_[_sorting_index_]);
	}
}

void pair_binner::_resort() const
{
	_sorting_index_ = 0;
	_sorted_ = false;
	_pair_bins_.clear();
	_sort();
}

	// Private implementations of set/clear methods
#if(1)

// Set specific limits through a limits vector
void pair_binner::_set_R_limits(std::vector< BRG_DISTANCE > R_bin_limits)
{
	if(R_bin_limits.empty())
	{
		_clear_R_limits();
	}
	else
	{
		_R_bin_limits_ = R_bin_limits;
	}
}
void pair_binner::_set_m_limits(std::vector< BRG_MASS > m_bin_limits)
{
	if(m_bin_limits.empty())
	{
		_clear_m_limits();
	}
	else
	{
		_m_bin_limits_ = m_bin_limits;
	}
}
void pair_binner::_set_z_limits(std::vector< double > z_bin_limits)
{
	if(z_bin_limits.empty())
	{
		_clear_z_limits();
	}
	else
	{
		_z_bin_limits_ = z_bin_limits;
	}

}
void pair_binner::_set_mag_limits(std::vector< double > mag_bin_limits)
{
	if(mag_bin_limits.empty())
	{
		_clear_mag_limits();
	}
	else
	{
		_mag_bin_limits_ = mag_bin_limits;
	}
}

// Set specific limits through a linear spacing
void pair_binner::_set_linear_R_limits(CONST_BRG_DISTANCE_REF R_min,
		CONST_BRG_DISTANCE_REF R_max,
		CONST_BRG_DISTANCE_REF R_step)
{
	_R_bin_limits_ = make_limit_vector(R_min, R_max, R_step);
}
void pair_binner::_set_linear_m_limits(CONST_BRG_MASS_REF m_min,
		CONST_BRG_MASS_REF m_max,
		CONST_BRG_MASS_REF m_step)
{
	_m_bin_limits_ = make_limit_vector(m_min, m_max, m_step);
}
void pair_binner::_set_linear_z_limits(double z_min,
		double z_max,
		double z_step)
{
	_z_bin_limits_ = make_limit_vector(z_min, z_max, z_step);
}
void pair_binner::_set_linear_mag_limits(double mag_min,
		double mag_max,
		double mag_step)
{
	_mag_bin_limits_ = make_limit_vector(mag_min, mag_max, mag_step);
}

// Set specific limits through a log spacing
void pair_binner::_set_log_R_limits(CONST_BRG_DISTANCE_REF R_min,
		CONST_BRG_DISTANCE_REF R_max,
		size_t R_num_bins)
{
	_R_bin_limits_ = make_log_limit_vector(R_min, R_max, R_num_bins);
}
void pair_binner::_set_log_m_limits(CONST_BRG_MASS_REF m_min,
		CONST_BRG_MASS_REF m_max,
		size_t m_num_bins)
{
	_m_bin_limits_ = make_log_limit_vector(m_min, m_max, m_num_bins);
}
void pair_binner::_set_log_z_limits(double z_min,
		double z_max,
		size_t z_num_bins)
{
	_z_bin_limits_ = make_log_limit_vector(z_min, z_max, z_num_bins);
}
void pair_binner::_set_log_mag_limits(double mag_min,
		double mag_max,
		size_t mag_num_bins)
{
	_mag_bin_limits_ = make_log_limit_vector(mag_min, mag_max, mag_num_bins);
}

// Clear limits. That is, make them unbound - one bin from neg infinity to pos infinity
void pair_binner::_clear_R_limits()
{
	_R_bin_limits_.clear();
	_R_bin_limits_.push_back(-std::numeric_limits<double>::infinity());
	_R_bin_limits_.push_back(std::numeric_limits<double>::infinity());
}
void pair_binner::_clear_m_limits()
{
	_m_bin_limits_.clear();
	_m_bin_limits_.push_back(-std::numeric_limits<double>::infinity());
	_m_bin_limits_.push_back(std::numeric_limits<double>::infinity());
}
void pair_binner::_clear_z_limits()
{
	_z_bin_limits_.clear();
	_z_bin_limits_.push_back(-std::numeric_limits<double>::infinity());
	_z_bin_limits_.push_back(std::numeric_limits<double>::infinity());
}
void pair_binner::_clear_mag_limits()
{
	_mag_bin_limits_.clear();
	_mag_bin_limits_.push_back(-std::numeric_limits<double>::infinity());
	_mag_bin_limits_.push_back(std::numeric_limits<double>::infinity());
}

#endif // Private implementations of set/clear methods

#endif // Private methods

	// Constructors
#if(1)

	// Set limits by vectors
pair_binner::pair_binner(std::vector< BRG_DISTANCE > R_bin_limits,
				std::vector< BRG_MASS > m_bin_limits,
				std::vector< double > z_bin_limits,
				std::vector< double > mag_bin_limits)
:	_R_bin_limits_(R_bin_limits),
 	_m_bin_limits_(m_bin_limits),
 	_z_bin_limits_(z_bin_limits),
 	_mag_bin_limits_(mag_bin_limits),
 	_valid_limits_(false),
 	_sorting_index_(0),
	_sorted_(false)
{
	_check_limits();
}

	// Set limits by min, max, and step
pair_binner::pair_binner(CONST_BRG_DISTANCE_REF R_min,
				CONST_BRG_DISTANCE_REF R_max,
				CONST_BRG_DISTANCE_REF R_step,
				CONST_BRG_MASS_REF m_min,
				CONST_BRG_MASS_REF m_max,
				CONST_BRG_MASS_REF m_step,
				double z_min,
				double z_max,
				double z_step,
				double mag_min,
				double mag_max,
				double mag_step)
:	_R_bin_limits_(make_limit_vector(R_min,R_max,R_step)),
 	_m_bin_limits_(make_limit_vector(m_min,m_max,m_step)),
 	_z_bin_limits_(make_limit_vector(z_min,z_max,z_step)),
 	_mag_bin_limits_(make_limit_vector(mag_min,mag_max,mag_step)),
 	_valid_limits_(false),
 	_sorting_index_(0),
	_sorted_(false)
{
	_check_limits();
}

#endif // Constructors

// Set/change limits
#if(1)

// Set specific limits through a limits vector
void pair_binner::set_R_limits(std::vector< BRG_DISTANCE > R_bin_limits)
{
	_set_R_limits(R_bin_limits);
	_check_limits();
}
void pair_binner::set_m_limits(std::vector< BRG_MASS > m_bin_limits)
{
	_set_m_limits(m_bin_limits);
	_check_limits();
}
void pair_binner::set_z_limits(std::vector< double > z_bin_limits)
{
	_set_z_limits(z_bin_limits);
	_check_limits();
}
void pair_binner::set_mag_limits(std::vector< double > mag_bin_limits)
{
	_set_mag_limits(mag_bin_limits);
	_check_limits();
}

// Set specific limits through a linear spacing
void pair_binner::set_linear_R_limits(CONST_BRG_DISTANCE_REF R_min,
		CONST_BRG_DISTANCE_REF R_max,
		CONST_BRG_DISTANCE_REF R_step)
{
	_set_linear_R_limits(R_min,R_max,R_step);
	_check_limits();
}
void pair_binner::set_linear_m_limits(CONST_BRG_MASS_REF m_min,
		CONST_BRG_MASS_REF m_max,
		CONST_BRG_MASS_REF m_step)
{
	_set_linear_m_limits(m_min,m_max,m_step);
	_check_limits();
}
void pair_binner::set_linear_z_limits(double z_min,
		double z_max,
		double z_step)
{
	_set_linear_z_limits(z_min,z_max,z_step);
	_check_limits();
}
void pair_binner::set_linear_mag_limits(double mag_min,
		double mag_max,
		double mag_step)
{
	_set_linear_mag_limits(mag_min,mag_max,mag_step);
	_check_limits();
}

// Set specific limits through a log spacing
void pair_binner::set_log_R_limits(CONST_BRG_DISTANCE_REF R_min,
		CONST_BRG_DISTANCE_REF R_max,
		size_t R_num_bins)
{
	_set_log_R_limits(R_min,R_max,R_num_bins);
	_check_limits();
}
void pair_binner::set_log_m_limits(CONST_BRG_MASS_REF m_min,
		CONST_BRG_MASS_REF m_max,
		size_t m_num_bins)
{
	_set_log_m_limits(m_min,m_max,m_num_bins);
	_check_limits();
}
void pair_binner::set_log_z_limits(double z_min,
		double z_max,
		size_t z_num_bins)
{
	_set_log_z_limits(z_min,z_max,z_num_bins);
	_check_limits();
}
void pair_binner::set_log_mag_limits(double mag_min,
		double mag_max,
		size_t mag_num_bins)
{
	_set_log_mag_limits(mag_min,mag_max,mag_num_bins);
	_check_limits();
}

// Clear limits. That is, make them unbound - one bin from neg infinity to pos infinity
void pair_binner::clear_R_limits()
{
	_clear_R_limits();
	_check_limits();
}
void pair_binner::clear_m_limits()
{
	_clear_m_limits();
	_check_limits();
}
void pair_binner::clear_z_limits()
{
	_clear_z_limits();
	_check_limits();
}
void pair_binner::clear_mag_limits()
{
	_clear_mag_limits();
	_check_limits();
}

void pair_binner::set_limits(std::vector< BRG_DISTANCE > R_bin_limits,
			std::vector< BRG_MASS > m_bin_limits,
			std::vector< double > z_bin_limits,
			std::vector< double > mag_bin_limits)
{
	_set_R_limits(R_bin_limits);
	_set_m_limits(m_bin_limits);
	_set_z_limits(z_bin_limits);
	_set_mag_limits(mag_bin_limits);
	_check_limits();
}

void pair_binner::set_linear_limits(CONST_BRG_DISTANCE_REF R_min,
			CONST_BRG_DISTANCE_REF R_max,
			CONST_BRG_DISTANCE_REF R_step,
			CONST_BRG_MASS_REF m_min,
			CONST_BRG_MASS_REF m_max,
			CONST_BRG_MASS_REF m_step,
			double z_min,
			double z_max,
			double z_step,
			double mag_min,
			double mag_max,
			double mag_step)
{
	_set_linear_R_limits(R_min,R_max,R_step);
	_set_linear_m_limits(m_min,m_max,m_step);
	_set_linear_z_limits(z_min,z_max,z_step);
	_set_linear_mag_limits(mag_min,mag_max,mag_step);
	_check_limits();
}

void pair_binner::set_log_limits(CONST_BRG_DISTANCE_REF R_min,
			CONST_BRG_DISTANCE_REF R_max,
			size_t R_num_bins,
			CONST_BRG_MASS_REF m_min,
			CONST_BRG_MASS_REF m_max,
			size_t m_num_bins,
			double z_min,
			double z_max,
			size_t z_num_bins,
			double mag_min,
			double mag_max,
			size_t mag_num_bins)
{
	_set_log_R_limits(R_min,R_max,R_num_bins);
	_set_log_m_limits(m_min,m_max,m_num_bins);
	_set_log_z_limits(z_min,z_max,z_num_bins);
	_set_log_mag_limits(mag_min,mag_max,mag_num_bins);
	_check_limits();
}

#endif

// Adding and clearing data
#if(1)

void pair_binner::add_pair( const lens_source_pair & new_pair)
{
	_pairs_.push_back(new_pair);
	_sorted_ = false;
}
void pair_binner::clear_pairs()
{
	_pair_bins_.clear();
	_pairs_.clear();
}
void pair_binner::sort()
{
	_sort();
}

#endif // Adding and clearing data

// Accessing summary data for bins
#if(1)

// Access by index (will throw if out of bounds)
#if(1)
BRG_UNITS pair_binner::delta_Sigma_t_mean_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	_sort();
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_mean();
}
BRG_UNITS pair_binner::delta_Sigma_x_mean_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	_sort();
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_x_mean();
}

BRG_UNITS pair_binner::delta_Sigma_t_std_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	_sort();
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_std();
}
BRG_UNITS pair_binner::delta_Sigma_x_std_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	_sort();
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_x_std();
}

BRG_UNITS pair_binner::delta_Sigma_t_stderr_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	_sort();
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_stderr();
}
BRG_UNITS pair_binner::delta_Sigma_x_stderr_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	_sort();
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_stderr();
}
#endif // Access by index

// Access by position
#if(1)
BRG_UNITS pair_binner::delta_Sigma_t_mean_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	_sort();
	size_t R_i = get_bin_index(R,_R_bin_limits_);
	size_t m_i = get_bin_index(m,_m_bin_limits_);
	size_t z_i = get_bin_index(z,_z_bin_limits_);
	size_t mag_i = get_bin_index(mag,_mag_bin_limits_);
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_mean();
}
BRG_UNITS pair_binner::delta_Sigma_x_mean_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	_sort();
	size_t R_i = get_bin_index(R,_R_bin_limits_);
	size_t m_i = get_bin_index(m,_m_bin_limits_);
	size_t z_i = get_bin_index(z,_z_bin_limits_);
	size_t mag_i = get_bin_index(mag,_mag_bin_limits_);
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_x_mean();
}

BRG_UNITS pair_binner::delta_Sigma_t_std_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	_sort();
	size_t R_i = get_bin_index(R,_R_bin_limits_);
	size_t m_i = get_bin_index(m,_m_bin_limits_);
	size_t z_i = get_bin_index(z,_z_bin_limits_);
	size_t mag_i = get_bin_index(mag,_mag_bin_limits_);
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_std();
}
BRG_UNITS pair_binner::delta_Sigma_x_std_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	_sort();
	size_t R_i = get_bin_index(R,_R_bin_limits_);
	size_t m_i = get_bin_index(m,_m_bin_limits_);
	size_t z_i = get_bin_index(z,_z_bin_limits_);
	size_t mag_i = get_bin_index(mag,_mag_bin_limits_);
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_x_std();
}

BRG_UNITS pair_binner::delta_Sigma_t_stderr_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	_sort();
	size_t R_i = get_bin_index(R,_R_bin_limits_);
	size_t m_i = get_bin_index(m,_m_bin_limits_);
	size_t z_i = get_bin_index(z,_z_bin_limits_);
	size_t mag_i = get_bin_index(mag,_mag_bin_limits_);
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_stderr();
}
BRG_UNITS pair_binner::delta_Sigma_x_stderr_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	_sort();
	size_t R_i = get_bin_index(R,_R_bin_limits_);
	size_t m_i = get_bin_index(m,_m_bin_limits_);
	size_t z_i = get_bin_index(z,_z_bin_limits_);
	size_t mag_i = get_bin_index(mag,_mag_bin_limits_);
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_x_stderr();
}
#endif // Access by index

#endif // Accessing summary data for bins

// Print data for all bins
void pair_binner::print_bin_data(std::ostream &out)
{
	_sort();

	// Set up the data and header to be printed
	table_t<double> data;
	header_t header;

	size_t num_columns = 19;

	header.resize(num_columns);
	header[0] = "R_min";
	header[1] = "R_max";
	header[2] = "R_mean";
	header[3] = "m_min";
	header[4] = "m_max";
	header[5] = "m_mean";
	header[6] = "z_min";
	header[7] = "z_max";
	header[8] = "z_mean";
	header[9] = "mag_min";
	header[10]= "mag_max";
	header[11]= "mag_mean";
	header[12]= "N";
	header[13]= "dS_t_mean";
	header[14]= "dS_x_mean";
	header[15]= "dS_t_stddev";
	header[16]= "dS_x_stddev";
	header[17]= "dS_t_stderr";
	header[18]= "dS_x_stderr";

	data.resize(num_columns);
	for(size_t R_i=0; R_i<_pair_bins_.size(); ++R_i)
	{
		for(size_t m_i=0; m_i<_pair_bins_[R_i].size(); ++m_i)
		{
			for(size_t z_i=0; z_i<_pair_bins_[R_i][m_i].size(); ++z_i)
			{
				for(size_t mag_i=0; mag_i<_pair_bins_[R_i][m_i][z_i].size(); ++mag_i)
				{
					pair_bin & bin = _pair_bins_[R_i][m_i][z_i][mag_i];
					data[0].push_back(bin.R_min());
					data[1].push_back(bin.R_max());
					data[2].push_back(bin.R_mean());
					data[3].push_back(bin.m_min());
					data[4].push_back(bin.m_max());
					data[5].push_back(bin.m_mean());
					data[6].push_back(bin.z_min());
					data[7].push_back(bin.z_max());
					data[8].push_back(bin.z_mean());
					data[9].push_back(bin.mag_min());
					data[10].push_back(bin.mag_max());
					data[11].push_back(bin.mag_mean());
					data[12].push_back(bin.count());
					data[13].push_back(bin.delta_Sigma_t_mean());
					data[14].push_back(bin.delta_Sigma_x_mean());
					data[15].push_back(bin.delta_Sigma_t_std());
					data[16].push_back(bin.delta_Sigma_x_std());
					data[17].push_back(bin.delta_Sigma_t_stderr());
					data[18].push_back(bin.delta_Sigma_x_stderr());
				}
			}
		}
	}

	// And now print it out
	print_table<double>(out,data,header);

}
void pair_binner::print_bin_data(std::string file_name)
{
	std::ofstream fo;
	open_file_output(fo,file_name);

	print_bin_data(fo);
}

} // end namespace brgastro
