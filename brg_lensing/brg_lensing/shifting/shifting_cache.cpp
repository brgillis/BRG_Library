/**********************************************************************\
 @file shifting_cache.cpp
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

#include "brg/common.h"

#include "brg/math/safe_math.hpp"

#include "brg_lensing/shifting/shifting_loader.h"
#include "brg_physics/units/unit_conversions.hpp"

#include "shifting_cache.h"

// Initialise the cache
DEFINE_BRG_CACHE_2D_STATIC_VARS( shifting_cache,
	0,3800*brgastro::unitconv::amintorad,1*brgastro::unitconv::amintorad,
	0.1,1.5,0.1);


// brgastro::tNFW_sig_cache class methods
#if (1)
const flt_type brgastro::shifting_cache::_calculate( const flt_type in_param_1, const flt_type in_param_2) const
{
	const flt_type t = std::fabs(in_param_1);
	const flt_type z = std::fabs(in_param_2);

	if(t==0) return 0;

	const flt_type zero_shift = shifting_loader().get(0,z);
	return safe_sqrt(shifting_loader().get(t,z)-zero_shift);
}

#endif // end brgastro::tNFW_sig_cache methods
