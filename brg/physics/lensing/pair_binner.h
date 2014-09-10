/**********************************************************************\
 @file pair_binner.h
 ------------------

 A class which stores data on lens-source pairs and sorts them into
 bins.

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

// body file: pair_binner.cpp

#ifndef _BRG_PAIR_BINNER_H_INCLUDED_
#define _BRG_PAIR_BINNER_H_INCLUDED_

#include <limits>
#include <vector>

#include <boost/tuple/tuple.hpp>

#include "brg/global.h"

#include "brg/physics/lensing/pair_bin.h"
#include "brg/physics/units/unit_obj.h"

namespace brgastro {

/**
 *
 */
class pair_binner {
private:

	// The bins, contained in a 4D vector, for one dimension for each parameter
	mutable std::vector< std::vector< std::vector< std::vector<pair_bin> > > > _pair_bins_;

	// Limits for the bins
#if(1)

	std::vector< BRG_DISTANCE > _R_bin_limits_;
	std::vector< BRG_MASS > _m_bin_limits_;
	std::vector< double > _z_bin_limits_;
	std::vector< double > _mag_bin_limits_;

	bool _valid_limits_;

#endif // Limits for the bins

	// Data to be sorted into bins
#if(1)

	std::vector<lens_source_pair> _pairs_;
	mutable size_t _sorting_index_;
	mutable bool _sorted_;

#endif

	// Private methods
#if(1)

	void _check_limits();
	void _sort() const;
	void _resort() const;

	// Private implementations of setting/clearing limits
	// Unlike the public versions, these don't check for valid
	// limits at the end.
#if(1)

	// Set specific limits through a limits vector
	void _set_R_limits(std::vector< BRG_DISTANCE > R_bin_limits);
	void _set_m_limits(std::vector< BRG_MASS > m_bin_limits);
	void _set_z_limits(std::vector< double > z_bin_limits);
	void _set_mag_limits(std::vector< double > mag_bin_limits);

	// Set specific limits through a linear spacing
	void _set_linear_R_limits(CONST_BRG_DISTANCE_REF R_min,
			CONST_BRG_DISTANCE_REF R_max,
			CONST_BRG_DISTANCE_REF R_step);
	void _set_linear_m_limits(CONST_BRG_MASS_REF m_min,
			CONST_BRG_MASS_REF m_max,
			CONST_BRG_MASS_REF m_step);
	void _set_linear_z_limits(double z_min,
			double z_max,
			double z_step);
	void _set_linear_mag_limits(double mag_min,
			double mag_max,
			double mag_step);

	// Set specific limits through a log spacing
	void _set_log_R_limits(CONST_BRG_DISTANCE_REF R_min,
			CONST_BRG_DISTANCE_REF R_max,
			size_t R_num_bins=1);
	void _set_log_m_limits(CONST_BRG_MASS_REF m_min,
			CONST_BRG_MASS_REF m_max,
			size_t m_num_bins=1);
	void _set_log_z_limits(double z_min,
			double z_max,
			size_t z_num_bins=1);
	void _set_log_mag_limits(double mag_min,
			double mag_max,
			size_t mag_num_bins=1);

	// Clear limits. That is, make them unbound - one bin from neg infinity to pos infinity
	void _clear_R_limits();
	void _clear_m_limits();
	void _clear_z_limits();
	void _clear_mag_limits();

#endif

#endif

public:

	// Constructors and destructor
#if(1)
	// Default constructor
	pair_binner()
	: _valid_limits_(false),
	  _sorting_index_(0),
	  _sorted_(false)
	{
	}

	// Set limits by vectors
	pair_binner(std::vector< BRG_DISTANCE > R_bin_limits,
				std::vector< BRG_MASS > m_bin_limits=std::vector<double>(),
				std::vector< double > z_bin_limits=std::vector<double>(),
				std::vector< double > mag_bin_limits=std::vector<double>());

	// Set limits by min, max, and step
	pair_binner(CONST_BRG_DISTANCE_REF R_min,
				CONST_BRG_DISTANCE_REF R_max,
				CONST_BRG_DISTANCE_REF R_step,
				CONST_BRG_MASS_REF m_min=-std::numeric_limits<double>::infinity(),
				CONST_BRG_MASS_REF m_max=std::numeric_limits<double>::infinity(),
				CONST_BRG_MASS_REF m_step=std::numeric_limits<double>::infinity(),
				double z_min=-std::numeric_limits<double>::infinity(),
				double z_max=std::numeric_limits<double>::infinity(),
				double z_step=std::numeric_limits<double>::infinity(),
				double mag_min=-std::numeric_limits<double>::infinity(),
				double mag_max=std::numeric_limits<double>::infinity(),
				double mag_step=std::numeric_limits<double>::infinity());


	virtual ~pair_binner()
	{
	}
#endif

	// Set/change limits
#if(1)

	// Set specific limits through a limits vector
	void set_R_limits(std::vector< BRG_DISTANCE > R_bin_limits);
	void set_m_limits(std::vector< BRG_MASS > m_bin_limits);
	void set_z_limits(std::vector< double > z_bin_limits);
	void set_mag_limits(std::vector< double > mag_bin_limits);

	// Set specific limits through a linear spacing
	void set_linear_R_limits(CONST_BRG_DISTANCE_REF R_min,
			CONST_BRG_DISTANCE_REF R_max,
			CONST_BRG_DISTANCE_REF R_step);
	void set_linear_m_limits(CONST_BRG_MASS_REF m_min,
			CONST_BRG_MASS_REF m_max,
			CONST_BRG_MASS_REF m_step);
	void set_linear_z_limits(double z_min,
			double z_max,
			double z_step);
	void set_linear_mag_limits(double mag_min,
			double mag_max,
			double mag_step);

	// Set specific limits through a log spacing
	void set_log_R_limits(CONST_BRG_DISTANCE_REF R_min,
			CONST_BRG_DISTANCE_REF R_max,
			size_t R_num_bins=1);
	void set_log_m_limits(CONST_BRG_MASS_REF m_min,
			CONST_BRG_MASS_REF m_max,
			size_t m_num_bins=1);
	void set_log_z_limits(double z_min,
			double z_max,
			size_t z_num_bins=1);
	void set_log_mag_limits(double mag_min,
			double mag_max,
			size_t mag_num_bins=1);

	// Clear limits. That is, make them unbound - one bin from neg infinity to pos infinity
	void clear_R_limits();
	void clear_m_limits();
	void clear_z_limits();
	void clear_mag_limits();

	void set_limits(std::vector< BRG_DISTANCE > R_bin_limits,
				std::vector< BRG_MASS > m_bin_limits=std::vector<double>(),
				std::vector< double > z_bin_limits=std::vector<double>(),
				std::vector< double > mag_bin_limits=std::vector<double>());

	void set_linear_limits(CONST_BRG_DISTANCE_REF R_min,
				CONST_BRG_DISTANCE_REF R_max,
				CONST_BRG_DISTANCE_REF R_step,
				CONST_BRG_MASS_REF m_min=-std::numeric_limits<double>::infinity(),
				CONST_BRG_MASS_REF m_max=std::numeric_limits<double>::infinity(),
				CONST_BRG_MASS_REF m_step=std::numeric_limits<double>::infinity(),
				double z_min=-std::numeric_limits<double>::infinity(),
				double z_max=std::numeric_limits<double>::infinity(),
				double z_step=std::numeric_limits<double>::infinity(),
				double mag_min=-std::numeric_limits<double>::infinity(),
				double mag_max=std::numeric_limits<double>::infinity(),
				double mag_step=std::numeric_limits<double>::infinity());

	void set_log_limits(CONST_BRG_DISTANCE_REF R_min,
				CONST_BRG_DISTANCE_REF R_max,
				size_t R_num_bins=1,
				CONST_BRG_MASS_REF m_min=-std::numeric_limits<double>::infinity(),
				CONST_BRG_MASS_REF m_max=std::numeric_limits<double>::infinity(),
				size_t m_num_bins=1,
				double z_min=-std::numeric_limits<double>::infinity(),
				double z_max=std::numeric_limits<double>::infinity(),
				size_t z_num_bins=1,
				double mag_min=-std::numeric_limits<double>::infinity(),
				double mag_max=std::numeric_limits<double>::infinity(),
				size_t mag_num_bins=1);

#endif // Set/change limits

	// Adding and clearing data
#if(1)

	void add_pair( const lens_source_pair & new_pair);
	void clear_pairs();

#endif // Adding and clearing data

	// Accessing summary data for bins
#if(1)

	// Access by index (will throw if out of bounds)
#if(1)
	BRG_UNITS delta_Sigma_t_mean_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i);
	BRG_UNITS delta_Sigma_x_mean_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i);

	BRG_UNITS delta_Sigma_t_std_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i);
	BRG_UNITS delta_Sigma_x_std_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i);

	BRG_UNITS delta_Sigma_t_stderr_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i);
	BRG_UNITS delta_Sigma_x_stderr_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i);
#endif // Access by index

	// Access by position
#if(1)
	BRG_UNITS delta_Sigma_t_mean_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
			double z, double mag);
	BRG_UNITS delta_Sigma_x_mean_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
			double z, double mag);

	BRG_UNITS delta_Sigma_t_std_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
			double z, double mag);
	BRG_UNITS delta_Sigma_x_std_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
			double z, double mag);

	BRG_UNITS delta_Sigma_t_stderr_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
			double z, double mag);
	BRG_UNITS delta_Sigma_x_stderr_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
			double z, double mag);
#endif // Access by index

#endif // Accessing summary data for bins
};

// Function template implementations
template<typename tup_type, typename val_type>
bool t1first_lt_v2(tup_type tup, val_type val)
{
	return boost::get<0>(tup) < val;
}

} // end namespace brgastro

#endif // _BRG_PAIR_BINNER_H_INCLUDED_
