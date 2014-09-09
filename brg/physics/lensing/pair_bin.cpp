/**********************************************************************\
 @file pair_bin.cpp
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
#include <brg/physics/lensing/pair_bin.h>

namespace brgastro {

// Calculations on stored values
#if (1)

BRG_UNITS pair_bin::delta_Sigma_t_mean()
{
	return mean(_delta_Sigma_t_values_);
}
BRG_UNITS pair_bin::delta_Sigma_x_mean()
{
	return mean(_delta_Sigma_x_values_);
}

BRG_UNITS pair_bin::delta_Sigma_t_std()
{
	return std(_delta_Sigma_t_values_);
}
BRG_UNITS pair_bin::delta_Sigma_x_std()
{
	return std(_delta_Sigma_x_values_);
}

BRG_UNITS pair_bin::delta_Sigma_t_stderr()
{
	return stderr(_delta_Sigma_t_values_);
}
BRG_UNITS pair_bin::delta_Sigma_x_stderr()
{
	return stderr(_delta_Sigma_x_values_);
}

#endif

} // end namespace brgastro
