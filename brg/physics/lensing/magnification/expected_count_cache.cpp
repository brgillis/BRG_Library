/**********************************************************************\
 @file expected_count_cache.cpp
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

#include <cstdlib>

#include "brg/global.h"

#include "brg/math/calculus/integrate.hpp"
#include "brg/math/safe_math.hpp"
#include "brg/physics/lensing/magnification/expected_count_loader.h"
#include "brg/physics/lensing/magnification/mag_global_values.h"
#include "brg/physics/lensing/magnification/Schechter_like_functor.h"

#include "expected_count_cache.h"

// Initialise the cache
DEFINE_BRG_CACHE_2D_STATIC_VARS( expected_count_cache,
		brgastro::mag_m_min,brgastro::mag_m_max,0.01,
		brgastro::mag_z_min,brgastro::mag_z_max-0.01,0.01);

// brgastro::expected_count_cache class methods
#if (1)
const double brgastro::expected_count_cache::_calculate( const double in_param_1, const double in_param_2) const
{
	const double m = std::fabs(in_param_1);
	const long double z_min = std::fabs(in_param_2);

	// We have to integrate over z from the input value to the max here

	auto f = [&] (long double z, const bool silent)
		{
			return Schechter_like_functor(expected_count_loader::get(z))(m);
		};

	return max(integrate_Romberg(&f,z_min,static_cast<long double>(brgastro::mag_z_max)),0);
}
#endif
