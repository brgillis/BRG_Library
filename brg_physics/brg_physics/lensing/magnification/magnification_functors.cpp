/**********************************************************************\
 @file magnification_functors.cpp
 ------------------

 Source file for the classes defined in magnification_functors.h.

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

#include "brg/global.h"

#include "brg/math/misc_math.hpp"
#include "brg/physics/astro.h"
#include "brg/physics/lensing/magnification/expected_count_cache.h"
#include "brg/physics/lensing/magnification/magnification_alpha.h"
#include "brg/physics/lensing/magnification/magnification_functors.h"

namespace brgastro {

constexpr double test_fudge_factor = 1;

BRG_UNITS mag_expected_count_functor::operator() (const long double & m, bool silent) const
{
	return expected_count_cache().get(m,_z_mean_);
}

BRG_UNITS mu_signal_integration_functor::operator() (const long double & m, bool silent) const
{
	long double alpha = magnification_alpha(m,_z_mean_);
	long double count = test_fudge_factor*expected_count_cache().get(m,_z_mean_);
	return count*(alpha-1)*(alpha-2);
}

BRG_UNITS mu_weight_integration_functor::operator() (const long double & m, bool silent) const
{
	long double alpha = magnification_alpha(m,_z_mean_);
	long double count = test_fudge_factor*expected_count_cache().get(m,_z_mean_);
	return count*square(alpha-1);
}

} // namespace brgastro

