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

#include <functional>
#include <limits>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <boost/container/flat_set.hpp>
#include <boost/optional.hpp>

#include "brg/global.h"

#include "brg_lensing/pair_bin.h"
#include "brg_lensing/pair_bins_summary.h"

#include "brg_physics/sky_obj/galaxy.h"
#include "brg_physics/units/unit_obj.h"

namespace brgastro {

struct lens_id_lt
{
	bool operator()(const lens_id & lens1, const lens_id & lens2) const
	{
		return lens1.id<lens2.id;
	}
};

/**
 *
 */
class pair_binner: public pair_bins_summary {
private:

	double _z_buffer_;
	static constexpr double _default_z_buffer_ = 0.1;

	// Data to be sorted into bins
#if(1)

	// The ID of the lens we're currently buffering (and its associated data), if we're indeed currently buffering one
	mutable boost::optional<lens_id> _buffering_lens_id_;

	// A buffer for the pair bins (It's recreated each time we sort, but we store the structure since it's rather large,
	// to save on allocation/deallocation calls).
	mutable std::vector<pair_bin> _pair_bins_buffer_;

	// Pairs which have yet to be sorted into bins
	mutable std::vector<lens_source_pair> _pairs_;
	mutable bool _sorted_;

#endif

	mutable boost::container::flat_set<lens_id,lens_id_lt> _past_lens_ids_;

	// Data on the unmasked fraction of annuli
#if(1)

	std::vector<BRG_DISTANCE> _unmasked_frac_bin_limits_;
	std::vector<double> _unmasked_fracs_;

#endif

	// Private methods
#if(1)
	void _sort() const;
	void _empty_buffer() const;
#endif

public:

	// Constructors and destructor
#if(1)
	// Default constructor
	pair_binner()
	: _z_buffer_(_default_z_buffer_),
	  _sorted_(false)
	{
	}

	// Set limits by vectors
	pair_binner(brgastro::limit_vector< BRG_DISTANCE > R_bin_limits,
			brgastro::limit_vector< BRG_MASS > m_bin_limits=brgastro::limit_vector<double>(),
			brgastro::limit_vector< double > z_bin_limits=brgastro::limit_vector<double>(),
			brgastro::limit_vector< double > mag_bin_limits=brgastro::limit_vector<double>())
	:	pair_bins_summary(R_bin_limits,m_bin_limits,z_bin_limits,mag_bin_limits),
		_z_buffer_(_default_z_buffer_),
	 	_sorted_(false)
	{
	}

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
				double mag_step=std::numeric_limits<double>::infinity())
	:	pair_bins_summary(R_min,R_max,R_step,m_min,m_max,m_step,z_min,z_max,z_step,
			mag_min,mag_max,mag_step),
			_z_buffer_(_default_z_buffer_),
		 	_sorted_(false)
	{
	}

	virtual ~pair_binner()
	{
	}
#endif

	// Setting and accessing z_buffer
#if(1)
	void set_z_buffer(const double & new_z_buffer)
	{
		_z_buffer_ = new_z_buffer;
	}
	double z_buffer() const
	{
		return _z_buffer_;
	}
#endif

	// Adding, sorting, and clearing data
#if(1)

	bool binnable( const galaxy & lens) const;
	void add_pair( const lens_source_pair & new_pair);
	void add_lens_id( const size_t & new_lens_id, const BRG_MASS & m, const double & z,
			const double & mag, const double & weight=1);
	void clear();
	void empty();
	void sort() const;

	template<typename Tv1, typename Tv2>
	void set_unmasked_fractions( Tv1 && bin_limits,
			Tv2 && unmasked_fractions)
	{
		// Sort at this point, so any lenses already added get the right unmasked fraction
		_sort();

		_unmasked_frac_bin_limits_ = std::forward<Tv1>(bin_limits);
		_unmasked_fracs_ = std::forward<Tv2>(unmasked_fractions);
	}

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
										   const double & z, const double & mag);
	BRG_UNITS delta_Sigma_x_mean_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
										   const double & z, const double & mag);

	BRG_UNITS delta_Sigma_t_std_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
										   const double & z, const double & mag);
	BRG_UNITS delta_Sigma_x_std_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
										   const double & z, const double & mag);

	BRG_UNITS delta_Sigma_t_stderr_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
										   const double & z, const double & mag);
	BRG_UNITS delta_Sigma_x_stderr_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
										   const double & z, const double & mag);
#endif // Access by index

#endif // Accessing summary data for bins
};

} // end namespace brgastro

#endif // _BRG_PAIR_BINNER_H_INCLUDED_
