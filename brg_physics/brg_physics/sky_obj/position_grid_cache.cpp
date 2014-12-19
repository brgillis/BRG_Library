/**********************************************************************\
  @file position_grid_cache.cpp

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

#include "brg_physics/units/unit_obj.h"

#include "position_grid_cache.h"

// Initialisation for brgastro::grid_cache
#if (1)

unsigned int brgastro::grid_cache::_ra_grid_change_num_ = 0;
unsigned int brgastro::grid_cache::_dec_grid_change_num_ = 0;
unsigned int brgastro::grid_cache::_z_grid_change_num_ = 0;
BRG_ANGLE brgastro::grid_cache::_ra_grid_min_ = -pi;
BRG_ANGLE brgastro::grid_cache::_ra_grid_max_ = pi;
BRG_ANGLE brgastro::grid_cache::_ra_grid_step_ = pi / 8;
BRG_ANGLE brgastro::grid_cache::_dec_grid_min_ = -pi / 2;
BRG_ANGLE brgastro::grid_cache::_dec_grid_max_ = pi / 2;
BRG_ANGLE brgastro::grid_cache::_dec_grid_step_ = pi / 8;
double brgastro::grid_cache::_z_grid_min_ = 0;
double brgastro::grid_cache::_z_grid_max_ = 2;
double brgastro::grid_cache::_z_grid_step_ = 0.1;

#endif // End initialisation for brgastro::grid_cache

