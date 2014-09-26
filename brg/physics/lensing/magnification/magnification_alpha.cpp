/**********************************************************************\
 @file magnification_alpha.cpp
 ------------------

 File containing the implementation of magnification_alpha.

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

#include "brg/math/calculus/differentiate.hpp"
#include "brg/physics/lensing/magnification/magnification_alpha.h"
#include "brg/physics/lensing/magnification/magnification_functors.h"

namespace brgastro {

double magnification_alpha(double m, double z)
{
	mag_expected_count_functor func(z);
	if(func(m)<=0) return 1.; // 1 corresponds to a flat profile, so it's appropriate for zero counts

	auto ln = [&] (long double mp, bool silent=true)
			{return std::log(func(mp,silent));};
	// Area cancels out in calculation
	auto result = 2.5*differentiate(&ln,m);
	return result;
}

} // namespace brgastro
