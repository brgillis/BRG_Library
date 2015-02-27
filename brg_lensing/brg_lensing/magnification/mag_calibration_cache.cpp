/**********************************************************************\
 @file mag_calibration_cache.cpp
 ------------------

 TODO <Insert file description here>

 **********************************************************************

 Copyright (C) 2015  Bryan R. Gillis

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

#include <cstdlib>

#include "brg/global.h"

#include "brg/math/misc_math.hpp"

#include "brg_lensing/magnification/mag_calibration_loader.h"
#include "brg_lensing/magnification/mag_global_values.h"
#include "brg_lensing/magnification/magnification_alpha.h"

#include "mag_calibration_cache.h"

// Initialise the cache
DEFINE_BRG_CACHE_STATIC_VARS( mag_calibration_cache,
							  brgastro::mag_z_min,brgastro::mag_z_max,brgastro::mag_z_step);

// brgastro::mag_calibration_cache class methods
#if (1)
const double brgastro::mag_calibration_cache::_calculate( const double in_param_1 ) const
{
	return mag_calibration_loader::get(in_param_1);
}
#endif


