/**********************************************************************\
 @file mag_signal_integral_cache.cpp
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
#include "brg/physics/lensing/magnification/expected_count_cache.h"
#include "brg/physics/lensing/magnification/mag_global_values.h"
#include "brg/physics/lensing/magnification/magnification_functors.h"

#include "mag_signal_integral_cache.h"

// Initialise the cache
DEFINE_BRG_CACHE_STATIC_VARS( mag_signal_integral_cache,
		0.2,brgastro::mag_z_max-0.01,0.01);

// brgastro::mag_signal_integral_cache class methods
#if (1)
const double brgastro::mag_signal_integral_cache::_calculate( const double in_param_1 ) const
{
	mu_signal_integration_functor func(in_param_1);

	return integrate_Romberg<decltype(func),long double>(&func,brgastro::mag_m_min,brgastro::mag_m_max);
}
void brgastro::mag_signal_integral_cache::_load_cache_dependencies() const
{
	expected_count_cache().get(0,0);
}
#endif