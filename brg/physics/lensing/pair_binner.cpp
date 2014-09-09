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

#include <vector>

#include "brg/global.h"

#include "brg/physics/units/unit_obj.h"
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

#endif

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
 	_sorting_index_(0)
{
	_check_limits();
}

	// Set limits by min, max, and step
pair_binner::pair_binner(CONST_BRG_DISTANCE_REF R_min,
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
:	_R_bin_limits_(make_limit_vector(R_min,R_max,R_step)),
 	_R_bin_limits_(make_limit_vector(m_min,m_max,m_step)),
 	_R_bin_limits_(make_limit_vector(z_min,z_max,z_step)),
 	_R_bin_limits_(make_limit_vector(mag_min,mag_max,mag_step)),
 	_valid_limits_(false),
 	_sorting_index_(0)
{
	_check_limits();
}

#endif // Constructors

} // end namespace brgastro
