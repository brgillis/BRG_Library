/**********************************************************************\
  @file astro_caches.cpp

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
#include <iostream>
#include <exception>
#include <string>

#include "IceBRG_main/common.h"

#include "IceBRG_main/math/cache/cache.hpp"
#include "IceBRG_main/math/cache/cache_2d.hpp"
#include "IceBRG_main/math/calculus/integrate.hpp"
#include "IceBRG_main/units/units.hpp"

#include "IceBRG_physics/astro.h"
#include "astro_caches.h"

/** Static Class Initialisation **/
#if (1)

namespace IceBRG {

// Initialisation for IceBRG::dfa_cache
DEFINE_BRG_CACHE( dfa_cache, flt_t,
		decltype(custom_unit_type<1, 0, 0, -1, 0>()), 0, 5, 0.001
		,
			return integrate_add( 0, in_param )/radian;
		,

);

// Initialisation for IceBRG::add_cache
DEFINE_BRG_CACHE_2D( add_cache, flt_t, flt_t, distance_type,
		0, 5, 0.01, 0, 5, 0.01,
			return IceBRG::integrate_add(in_param_1,in_param_2);
		,

);

// Initialisation for IceBRG::tfa_cache
DEFINE_BRG_CACHE( tfa_cache, flt_t, time_type, 0.001, 1.02, 0.001
		,
			return -integrate_ltd( 0, zfa( in_param ) ) / c;
		,

);

// Initialisation for IceBRG::lum_func_integral_cache
DEFINE_BRG_CACHE_2D( lum_func_integral_cache, flt_t, flt_t, decltype(custom_unit_type<-3,0,0,0,0>()),
		-25, -19, 0.1,
		-25, -19, 0.1
		,
			flt_t res = integrate_Romberg(differential_luminosity_function,in_param_1,in_param_2);
			return res;
		,

);

}

#endif // end Static Class Initialisation
