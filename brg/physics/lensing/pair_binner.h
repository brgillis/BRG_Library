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

#include <vector>

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
	std::vector< std::vector< std::vector< std::vector<pair_bin> > > > _pair_bins_;

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
	size_t _sorting_index_;

#endif

	// Private methods
#if(1)

	void _check_limits();

#endif

public:

	// Constructors and destructor
#if(1)
	// Default constructor
	pair_binner()
	: _valid_limits_(false),
	  _sorting_index_(0)
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
};

} // end namespace brgastro

#endif // _BRG_PAIR_BINNER_H_INCLUDED_
