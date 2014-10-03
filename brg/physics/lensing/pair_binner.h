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

#include <boost/tuple/tuple.hpp>

#include "brg/global.h"

#include "brg/physics/sky_obj/galaxy.h"
#include "brg/physics/lensing/pair_bin.h"
#include "brg/physics/lensing/pair_bins_summary.h"
#include "brg/physics/units/unit_obj.h"

namespace brgastro {

// lens_id struct
#if(1)
struct lens_id
{
	size_t id;
	BRG_MASS m;
	double z;
	double mag;

	lens_id(const size_t & id, const BRG_MASS & m, const double & z, const double & mag)
	: id(id),
	  m(m),
	  z(z),
	  mag(mag)
	{
	}
};

inline bool lens_id_lt(const lens_id & lens1, const lens_id & lens2)
{
	return lens1.id<lens2.id;
}
#endif

/**
 *
 */
class pair_binner: public pair_bins_summary {
private:

	// The bins, contained in a 4D vector, for one dimension for each parameter
	mutable std::vector< std::vector< std::vector< std::vector<pair_bin> > > > _pair_bins_;

	// Data to be sorted into bins
#if(1)

	std::vector<lens_source_pair> _pairs_;
	std::set<lens_id,decltype(&lens_id_lt)> _lens_ids_;
	mutable size_t _sorting_index_;
	mutable bool _sorted_;

#endif

	// Private methods
#if(1)
	void _sort() const;
	void _resort() const;
#endif

public:

	// Constructors and destructor
#if(1)
	// Default constructor
	pair_binner()
	: _sorting_index_(0),
	  _sorted_(false)
	{
	}

	// Set limits by vectors
	pair_binner(std::vector< BRG_DISTANCE > R_bin_limits,
				std::vector< BRG_MASS > m_bin_limits=std::vector<double>(),
				std::vector< double > z_bin_limits=std::vector<double>(),
				std::vector< double > mag_bin_limits=std::vector<double>())
	:	pair_bins_summary(R_bin_limits,m_bin_limits,z_bin_limits,mag_bin_limits),
	 	_lens_ids_(&lens_id_lt),
	 	_sorting_index_(0),
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
		 	_lens_ids_(&lens_id_lt),
		 	_sorting_index_(0),
		 	_sorted_(false)
	{
	}

	virtual ~pair_binner()
	{
	}
#endif

	// pair_bins accessor
#if(1)

	const std::vector< std::vector< std::vector< std::vector<pair_bin> > > > & pair_bins() const
	{
		return _pair_bins_;
	}

#endif

	// Adding, sorting, and clearing data
#if(1)

	bool binnable( const galaxy & lens) const;
	void add_pair( const lens_source_pair & new_pair);
	void add_lens_id( const size_t & new_lens_id, const BRG_MASS & m, const double & z,
			const double & mag);
	void clear();
	void empty();
	void sort() const;

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
