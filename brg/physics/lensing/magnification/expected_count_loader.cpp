/**********************************************************************\
 @file expected_count_loader.cpp
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

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/lexical_cast.hpp>

#include "brg/global.h"

#include "brg/file_access/ascii_table_map.hpp"
#include "brg/math/misc_math.hpp"
#include "brg/vector/limit_vector.hpp"

#include "brg/physics/lensing/magnification/expected_count_loader.h"
#include "brg/physics/lensing/magnification/mag_global_values.h"

namespace brgastro {

// Initialisation of static vars
#if (1)
bool expected_count_loader::_loaded_(false);

std::vector<double> expected_count_loader::_z_limits_(brgastro::make_limit_vector<double>(mag_z_min,
		mag_z_max,mag_z_step));
std::vector<std::vector<double>> expected_count_loader::_mag_limits_;
std::vector<std::vector<double>> expected_count_loader::_smoothed_count_;
std::vector<std::vector<double>> expected_count_loader::_smoothed_count_derivative_;

std::string expected_count_loader::_filename_base_("/disk2/brg/git/CFHTLenS_cat/Data/magnitude_hist_gz");
std::string expected_count_loader::_filename_tail_(".dat");

#endif

void expected_count_loader::_load()
{
	if(_loaded_) return;

	#ifdef _OPENMP
	#pragma omp critical(brg_expected_count_loader_load)
	#endif
	{

	if(!_loaded_)
	{

		// Resize arrays as necessary
		size_t num_z_bins = _z_limits_.size()-1;

		_mag_limits_.resize(num_z_bins);
		_smoothed_count_.resize(num_z_bins);
		_smoothed_count_derivative_.resize(num_z_bins);

		// Loop over z, loading in all files
		for(size_t z_i=0; z_i<num_z_bins; ++z_i)
		{
			const unsigned short z100(brgastro::round_int(100*_z_limits_[z_i]));
			const std::string filename(_filename_base_ + boost::lexical_cast<std::string>(z100)
					+ _filename_tail_);

			std::vector<double> count;

			try
			{
				auto table_map = brgastro::load_table_map<double>(filename);
				_mag_limits_[z_i] = table_map.at("mag_bin_lower");
				count = table_map.at("count");
				_smoothed_count_[z_i] = table_map.at("smoothed_count");
				_smoothed_count_derivative_[z_i] = table_map.at("smoothed_alpha");

				// Add a final value to the _mag_limits_ vector
				_mag_limits_[z_i].push_back(2*_mag_limits_[z_i].back()-
						_mag_limits_[z_i].at(_mag_limits_[z_i].size()-2));
			}
			catch(const std::runtime_error &e)
			{
				std::cerr << "ERROR: Cannot load data for expected_count_loader. Check that the filename\n"
						<< "base and root match the data (including path), and the z limits are set up to match.\n"
						<< "Original exception:" << e.what();
				throw;
			}

			// Normalize by magnitude bin size
			double bin_size = _mag_limits_[z_i].at(1)-_mag_limits_[z_i].at(0);
			_smoothed_count_[z_i] = brgastro::divide(_smoothed_count_[z_i],bin_size);

			// Correct spurious smoothed count values in this vector
			for(size_t m_i=0; m_i<_smoothed_count_.size(); ++m_i)
			{
				if((count[m_i]==0)&&(_smoothed_count_[z_i][m_i]==brgastro::mag_min_count))
				{
					_smoothed_count_[z_i][m_i] = 0;
				}
			}
		} // Loop over z, loading in all files

		_loaded_ = true;

	}

	}
} // void expected_count_loader::_load()


double expected_count_loader::_get_interp(const double & mag, const double & z,
		const std::vector<std::vector<double>> & table, const double & def)
{
	// Load if necessary
	_load();

	size_t z_i = brgastro::get_bin_index(z,_z_limits_);

	if(z_i==_z_limits_.size()-2) --z_i;
	const double & z_lo = _z_limits_[z_i];
	const double & z_hi = _z_limits_[z_i+1];

	double tot_weight;
	double temp_result;

	// Get the interpolated value at both the lower redshift and the higher

	// At the lower redshift first
	double lo_result;
	if(brgastro::outside_limits(mag,_mag_limits_[z_i]))
	{
		lo_result = def;
	}
	else
	{
		lo_result = interpolate_bins(mag,_mag_limits_[z_i],table[z_i]);
	}

	// At the higher redshift now
	double hi_result;

	if(brgastro::outside_limits(mag,_mag_limits_[z_i+1]))
	{
		hi_result = def;
	}
	else
	{
		hi_result = interpolate_bins(mag,_mag_limits_[z_i+1],table[z_i+1]);
	}

	// And now interpolate between these results

	tot_weight = z_hi-z_lo;

	temp_result = 0;
	temp_result += lo_result * (z_hi-z);
	temp_result += hi_result * (z-z_lo);
	return temp_result / tot_weight;

} // double expected_count_loader::_get_interp(...)

// Setting parameters for where the data can be loaded from
#if(1)
void expected_count_loader::set_z_limits(const std::vector<double> & new_limits_vector)
{
	if(!is_monotonically_increasing(new_limits_vector))
	{
		throw std::logic_error("New z_limits vector must be monotonically increasing.");
	}
	_z_limits_ = new_limits_vector;
}
void expected_count_loader::set_z_limits(std::vector<double> && new_limits_vector)
{
	if(!is_monotonically_increasing(new_limits_vector))
	{
		throw std::logic_error("New z_limits vector must be monotonically increasing.");
	}
	_z_limits_ = std::move(new_limits_vector);
}
void expected_count_loader::set_filename_base(const std::string & new_filename_base)
{
	_filename_base_ = new_filename_base;
}
void expected_count_loader::set_filename_base(std::string && new_filename_base)
{
	_filename_base_ = std::move(new_filename_base);
}
void expected_count_loader::set_filename_tail(const std::string & new_filename_tail)
{
	_filename_tail_ = new_filename_tail;
}
void expected_count_loader::set_filename_tail(std::string && new_filename_tail)
{
	_filename_tail_ = std::move(new_filename_tail);
}
#endif // Setting parameters for where the data is

// Access data
#if(1)

double expected_count_loader::get_count(const double & mag, const double & z)
{
	return _get_interp(mag,z,_smoothed_count_,0.);
}
double expected_count_loader::get_derivative(const double & mag, const double & z)
{
	return _get_interp(mag,z,_smoothed_count_derivative_,1.);
}

#endif // Access data

} // end namespace brgastro
