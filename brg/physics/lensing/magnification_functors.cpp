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

#include "brg/physics/astro.h"
#include "brg/physics/lensing/magnification_alpha.h"
#include "brg/physics/lensing/magnification_functors.h"
#include "brg/physics/units/unit_conversions.hpp"
#include "brg/physics/units/unit_obj.h"

namespace brgastro {

double mag_expected_count_functor::operator() (double m)
{
	// TODO Fill-in proper, calibrated value
	return 1;
}

double mu_signal_integration_functor::operator() (double m)
{
	double alpha = magnification_alpha(m);
	return mag_expected_count_functor(_area_)(m)*(alpha-1)*(alpha-2);
}

double mu_weight_integration_functor::operator() (double m)
{
	return mag_expected_count_functor(_area_)(m)*(magnification_alpha(m)-1);
}

} // namespace brgastro

