/**********************************************************************\
 @file pair_bins_summary.cpp
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
#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "brg/global.h"

#include "brg/file_access/ascii_table.h"
#include "brg/physics/lensing/pair_binner.h"
#include "brg/physics/units/unit_obj.h"
#include "brg/utility.hpp"
#include "brg/vector/limit_vector.hpp"
#include "brg/vector/summary_functions.hpp"

#include "pair_bins_summary.h"

namespace brgastro {

	// Private methods
#if(1)

void pair_bins_summary::_check_limits()
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

	// Private implementations of set/clear methods
#if(1)

// Set specific limits through a limits vector
void pair_bins_summary::_set_R_limits(std::vector< BRG_DISTANCE > R_bin_limits)
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
void pair_bins_summary::_set_m_limits(std::vector< BRG_MASS > m_bin_limits)
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
void pair_bins_summary::_set_z_limits(std::vector< double > z_bin_limits)
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
void pair_bins_summary::_set_mag_limits(std::vector< double > mag_bin_limits)
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
void pair_bins_summary::_set_linear_R_limits(CONST_BRG_DISTANCE_REF R_min,
		CONST_BRG_DISTANCE_REF R_max,
		CONST_BRG_DISTANCE_REF R_step)
{
	_R_bin_limits_ = make_limit_vector(R_min, R_max, R_step);
}
void pair_bins_summary::_set_linear_m_limits(CONST_BRG_MASS_REF m_min,
		CONST_BRG_MASS_REF m_max,
		CONST_BRG_MASS_REF m_step)
{
	_m_bin_limits_ = make_limit_vector(m_min, m_max, m_step);
}
void pair_bins_summary::_set_linear_z_limits(double z_min,
		double z_max,
		double z_step)
{
	_z_bin_limits_ = make_limit_vector(z_min, z_max, z_step);
}
void pair_bins_summary::_set_linear_mag_limits(double mag_min,
		double mag_max,
		double mag_step)
{
	_mag_bin_limits_ = make_limit_vector(mag_min, mag_max, mag_step);
}

// Set specific limits through a log spacing
void pair_bins_summary::_set_log_R_limits(CONST_BRG_DISTANCE_REF R_min,
		CONST_BRG_DISTANCE_REF R_max,
		size_t R_num_bins)
{
	_R_bin_limits_ = make_log_limit_vector(R_min, R_max, R_num_bins);
}
void pair_bins_summary::_set_log_m_limits(CONST_BRG_MASS_REF m_min,
		CONST_BRG_MASS_REF m_max,
		size_t m_num_bins)
{
	_m_bin_limits_ = make_log_limit_vector(m_min, m_max, m_num_bins);
}
void pair_bins_summary::_set_log_z_limits(double z_min,
		double z_max,
		size_t z_num_bins)
{
	_z_bin_limits_ = make_log_limit_vector(z_min, z_max, z_num_bins);
}
void pair_bins_summary::_set_log_mag_limits(double mag_min,
		double mag_max,
		size_t mag_num_bins)
{
	_mag_bin_limits_ = make_log_limit_vector(mag_min, mag_max, mag_num_bins);
}

// Clear limits. That is, make them unbound - one bin from neg infinity to pos infinity
void pair_bins_summary::_clear_R_limits()
{
	_R_bin_limits_.clear();
	_R_bin_limits_.push_back(-std::numeric_limits<double>::infinity());
	_R_bin_limits_.push_back(std::numeric_limits<double>::infinity());
}
void pair_bins_summary::_clear_m_limits()
{
	_m_bin_limits_.clear();
	_m_bin_limits_.push_back(-std::numeric_limits<double>::infinity());
	_m_bin_limits_.push_back(std::numeric_limits<double>::infinity());
}
void pair_bins_summary::_clear_z_limits()
{
	_z_bin_limits_.clear();
	_z_bin_limits_.push_back(-std::numeric_limits<double>::infinity());
	_z_bin_limits_.push_back(std::numeric_limits<double>::infinity());
}
void pair_bins_summary::_clear_mag_limits()
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
pair_bins_summary::pair_bins_summary(std::vector< BRG_DISTANCE > R_bin_limits,
				std::vector< BRG_MASS > m_bin_limits,
				std::vector< double > z_bin_limits,
				std::vector< double > mag_bin_limits)
:	_R_bin_limits_(R_bin_limits),
 	_m_bin_limits_(m_bin_limits),
 	_z_bin_limits_(z_bin_limits),
 	_mag_bin_limits_(mag_bin_limits),
 	_valid_limits_(false)
{
	_check_limits();
}

// Set limits by min, max, and step
pair_bins_summary::pair_bins_summary(CONST_BRG_DISTANCE_REF R_min,
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
 	_valid_limits_(false)
{
	_check_limits();
}

// Construct from pair_bins
pair_bins_summary::pair_bins_summary( const pair_binner & bins)
:	_R_bin_limits_(bins.R_limits()),
 	_m_bin_limits_(bins.m_limits()),
 	_z_bin_limits_(bins.z_limits()),
 	_mag_bin_limits_(bins.mag_limits()),
 	_valid_limits_(false)
{
	_check_limits();
	if(_valid_limits_)
	{
		bins.sort();
		make_array4d(_pair_bin_summaries_,bins.pair_bins());

		for(size_t R_i=0; R_i<_pair_bin_summaries_.size(); ++R_i)
		{
			for(size_t m_i=0; m_i<_pair_bin_summaries_[R_i].size(); ++m_i)
			{
				for(size_t z_i=0; z_i<_pair_bin_summaries_[R_i][m_i].size(); ++z_i)
				{
					for(size_t mag_i=0; mag_i<_pair_bin_summaries_[R_i][m_i][z_i].size(); ++mag_i)
					{
						_pair_bin_summaries_[R_i][m_i][z_i][mag_i] =
								bins.pair_bins()[R_i][m_i][z_i][mag_i];
					}
				}
			}
		}
	}
}

#endif // Constructors

// Set/change limits
#if(1)

// Set specific limits through a limits vector
void pair_bins_summary::set_R_limits(std::vector< BRG_DISTANCE > R_bin_limits)
{
	_set_R_limits(R_bin_limits);
	_check_limits();
}
void pair_bins_summary::set_m_limits(std::vector< BRG_MASS > m_bin_limits)
{
	_set_m_limits(m_bin_limits);
	_check_limits();
}
void pair_bins_summary::set_z_limits(std::vector< double > z_bin_limits)
{
	_set_z_limits(z_bin_limits);
	_check_limits();
}
void pair_bins_summary::set_mag_limits(std::vector< double > mag_bin_limits)
{
	_set_mag_limits(mag_bin_limits);
	_check_limits();
}

// Set specific limits through a linear spacing
void pair_bins_summary::set_linear_R_limits(CONST_BRG_DISTANCE_REF R_min,
		CONST_BRG_DISTANCE_REF R_max,
		CONST_BRG_DISTANCE_REF R_step)
{
	_set_linear_R_limits(R_min,R_max,R_step);
	_check_limits();
}
void pair_bins_summary::set_linear_m_limits(CONST_BRG_MASS_REF m_min,
		CONST_BRG_MASS_REF m_max,
		CONST_BRG_MASS_REF m_step)
{
	_set_linear_m_limits(m_min,m_max,m_step);
	_check_limits();
}
void pair_bins_summary::set_linear_z_limits(double z_min,
		double z_max,
		double z_step)
{
	_set_linear_z_limits(z_min,z_max,z_step);
	_check_limits();
}
void pair_bins_summary::set_linear_mag_limits(double mag_min,
		double mag_max,
		double mag_step)
{
	_set_linear_mag_limits(mag_min,mag_max,mag_step);
	_check_limits();
}

// Set specific limits through a log spacing
void pair_bins_summary::set_log_R_limits(CONST_BRG_DISTANCE_REF R_min,
		CONST_BRG_DISTANCE_REF R_max,
		size_t R_num_bins)
{
	_set_log_R_limits(R_min,R_max,R_num_bins);
	_check_limits();
}
void pair_bins_summary::set_log_m_limits(CONST_BRG_MASS_REF m_min,
		CONST_BRG_MASS_REF m_max,
		size_t m_num_bins)
{
	_set_log_m_limits(m_min,m_max,m_num_bins);
	_check_limits();
}
void pair_bins_summary::set_log_z_limits(double z_min,
		double z_max,
		size_t z_num_bins)
{
	_set_log_z_limits(z_min,z_max,z_num_bins);
	_check_limits();
}
void pair_bins_summary::set_log_mag_limits(double mag_min,
		double mag_max,
		size_t mag_num_bins)
{
	_set_log_mag_limits(mag_min,mag_max,mag_num_bins);
	_check_limits();
}

// Clear limits. That is, make them unbound - one bin from neg infinity to pos infinity
void pair_bins_summary::clear_R_limits()
{
	_clear_R_limits();
	_check_limits();
}
void pair_bins_summary::clear_m_limits()
{
	_clear_m_limits();
	_check_limits();
}
void pair_bins_summary::clear_z_limits()
{
	_clear_z_limits();
	_check_limits();
}
void pair_bins_summary::clear_mag_limits()
{
	_clear_mag_limits();
	_check_limits();
}

void pair_bins_summary::set_limits(std::vector< BRG_DISTANCE > R_bin_limits,
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

void pair_bins_summary::set_linear_limits(CONST_BRG_DISTANCE_REF R_min,
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

void pair_bins_summary::set_log_limits(CONST_BRG_DISTANCE_REF R_min,
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

void pair_bins_summary::clear()
{
	_pair_bin_summaries_.clear();
}

#endif // Adding and clearing data

// Accessing summary data for bins
#if(1)

// Access by index (will throw if out of bounds)
#if(1)
BRG_UNITS pair_bins_summary::delta_Sigma_t_mean_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	return _pair_bin_summaries_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_mean();
}
BRG_UNITS pair_bins_summary::delta_Sigma_x_mean_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	return _pair_bin_summaries_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_x_mean();
}

BRG_UNITS pair_bins_summary::delta_Sigma_t_std_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	return _pair_bin_summaries_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_std();
}
BRG_UNITS pair_bins_summary::delta_Sigma_x_std_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	return _pair_bin_summaries_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_x_std();
}

BRG_UNITS pair_bins_summary::delta_Sigma_t_stderr_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	return _pair_bin_summaries_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_stderr();
}
BRG_UNITS pair_bins_summary::delta_Sigma_x_stderr_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	return _pair_bin_summaries_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_stderr();
}
#endif // Access by index

// Access by position
#if(1)
BRG_UNITS pair_bins_summary::delta_Sigma_t_mean_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	size_t R_i = get_bin_index(R,_R_bin_limits_);
	size_t m_i = get_bin_index(m,_m_bin_limits_);
	size_t z_i = get_bin_index(z,_z_bin_limits_);
	size_t mag_i = get_bin_index(mag,_mag_bin_limits_);
	return _pair_bin_summaries_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_mean();
}
BRG_UNITS pair_bins_summary::delta_Sigma_x_mean_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	size_t R_i = get_bin_index(R,_R_bin_limits_);
	size_t m_i = get_bin_index(m,_m_bin_limits_);
	size_t z_i = get_bin_index(z,_z_bin_limits_);
	size_t mag_i = get_bin_index(mag,_mag_bin_limits_);
	return _pair_bin_summaries_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_x_mean();
}

BRG_UNITS pair_bins_summary::delta_Sigma_t_std_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	size_t R_i = get_bin_index(R,_R_bin_limits_);
	size_t m_i = get_bin_index(m,_m_bin_limits_);
	size_t z_i = get_bin_index(z,_z_bin_limits_);
	size_t mag_i = get_bin_index(mag,_mag_bin_limits_);
	return _pair_bin_summaries_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_std();
}
BRG_UNITS pair_bins_summary::delta_Sigma_x_std_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	size_t R_i = get_bin_index(R,_R_bin_limits_);
	size_t m_i = get_bin_index(m,_m_bin_limits_);
	size_t z_i = get_bin_index(z,_z_bin_limits_);
	size_t mag_i = get_bin_index(mag,_mag_bin_limits_);
	return _pair_bin_summaries_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_x_std();
}

BRG_UNITS pair_bins_summary::delta_Sigma_t_stderr_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	size_t R_i = get_bin_index(R,_R_bin_limits_);
	size_t m_i = get_bin_index(m,_m_bin_limits_);
	size_t z_i = get_bin_index(z,_z_bin_limits_);
	size_t mag_i = get_bin_index(mag,_mag_bin_limits_);
	return _pair_bin_summaries_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_stderr();
}
BRG_UNITS pair_bins_summary::delta_Sigma_x_stderr_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	size_t R_i = get_bin_index(R,_R_bin_limits_);
	size_t m_i = get_bin_index(m,_m_bin_limits_);
	size_t z_i = get_bin_index(z,_z_bin_limits_);
	size_t mag_i = get_bin_index(mag,_mag_bin_limits_);
	return _pair_bin_summaries_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_x_stderr();
}
#endif // Access by index

#endif // Accessing summary data for bins

// Print data for all bins
void pair_bins_summary::print_bin_data(std::ostream &out)
{
	// Set up the data and header to be printed
	table_t<double> data;
	header_t header;

	size_t num_columns = 22;

	header.reserve(num_columns);
	header.push_back("R_min");
	header.push_back("R_max");
	header.push_back("R_mean");
	header.push_back("m_min");
	header.push_back("m_max");
	header.push_back("m_mean");
	header.push_back("z_min");
	header.push_back("z_max");
	header.push_back("z_mean");
	header.push_back("mag_min");
	header.push_back("mag_max");
	header.push_back("mag_mean");
	header.push_back("N");
	header.push_back("N_eff");
	header.push_back("dS_t_mean");
	header.push_back("dS_x_mean");
	header.push_back("dS_t_stddev");
	header.push_back("dS_x_stddev");
	header.push_back("dS_t_stderr");
	header.push_back("dS_x_stderr");
	header.push_back("mu");
	header.push_back("mu_stderr");

	assert(header.size()==num_columns);

	data.resize(num_columns);
	for(size_t R_i=0; R_i<_pair_bin_summaries_.size(); ++R_i)
	{
		for(size_t m_i=0; m_i<_pair_bin_summaries_[R_i].size(); ++m_i)
		{
			for(size_t z_i=0; z_i<_pair_bin_summaries_[R_i][m_i].size(); ++z_i)
			{
				for(size_t mag_i=0; mag_i<_pair_bin_summaries_[R_i][m_i][z_i].size(); ++mag_i)
				{
					pair_bin_summary & bin = _pair_bin_summaries_[R_i][m_i][z_i][mag_i];

					size_t col_i = 0;
					data[col_i].push_back(bin.R_min());
					data[++col_i].push_back(bin.R_max());
					data[++col_i].push_back(bin.R_mean());
					data[++col_i].push_back(bin.m_min());
					data[++col_i].push_back(bin.m_max());
					data[++col_i].push_back(bin.m_mean());
					data[++col_i].push_back(bin.z_min());
					data[++col_i].push_back(bin.z_max());
					data[++col_i].push_back(bin.z_mean());
					data[++col_i].push_back(bin.mag_min());
					data[++col_i].push_back(bin.mag_max());
					data[++col_i].push_back(bin.mag_mean());
					data[++col_i].push_back(bin.count());
					data[++col_i].push_back(bin.effective_count());
					data[++col_i].push_back(bin.delta_Sigma_t_mean());
					data[++col_i].push_back(bin.delta_Sigma_x_mean());
					data[++col_i].push_back(bin.delta_Sigma_t_std());
					data[++col_i].push_back(bin.delta_Sigma_x_std());
					data[++col_i].push_back(bin.delta_Sigma_t_stderr());
					data[++col_i].push_back(bin.delta_Sigma_x_stderr());
					data[++col_i].push_back(bin.mu_hat());
					data[++col_i].push_back(1/bin.mu_W());
				}
			}
		}
	}

	// And now print it out
	print_table<double>(out,data,header);

}
void pair_bins_summary::print_bin_data(std::string file_name)
{
	std::ofstream fo;
	open_file_output(fo,file_name);

	print_bin_data(fo);
}

// Operators to combine data
#if (1)

pair_bins_summary & pair_bins_summary::operator+=( const pair_bins_summary & other )
{
	// Check for same size
	assert(_R_bin_limits_==other.R_limits());
	assert(_m_bin_limits_==other.m_limits());
	assert(_z_bin_limits_==other.z_limits());
	assert(_mag_bin_limits_==other.mag_limits());

	// Make sure the other is sorted
	other.sort();

	// And add the bins together

	for(size_t R_i=0; R_i<_pair_bin_summaries_.size(); ++R_i)
	{
		for(size_t m_i=0; m_i<_pair_bin_summaries_[R_i].size(); ++m_i)
		{
			for(size_t z_i=0; z_i<_pair_bin_summaries_[R_i][m_i].size(); ++z_i)
			{
				for(size_t mag_i=0; mag_i<_pair_bin_summaries_[R_i][m_i][z_i].size(); ++mag_i)
				{
					_pair_bin_summaries_[R_i][m_i][z_i][mag_i] +=
							other.pair_bin_summaries()[R_i][m_i][z_i][mag_i];
				}
			}
		}
	}

	return *this;
}

#endif

} // end namespace brgastro
