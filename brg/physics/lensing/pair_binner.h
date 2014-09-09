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

#include <stack>
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

#endif // Limits for the bins

	// Data to be sorted into bins
#if(1)

	std::stack<lens_source_pair> _sorted_pairs_;
	std::stack<lens_source_pair> _unsorted_pairs_;

#endif

public:
	pair_binner() {
		// TODO Auto-generated constructor stub

	}
	virtual ~pair_binner() {
		// TODO Auto-generated destructor stub
	}
};

} // end namespace brgastro

#endif // _BRG_PAIR_BINNER_H_INCLUDED_
