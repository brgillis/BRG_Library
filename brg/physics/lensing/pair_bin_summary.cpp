/**********************************************************************\
 @file pair_bin_summary.cpp
 ------------------

 Source file for the pair_bin class.

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

#include <cassert>

#include "brg/global.h"

#include "brg/physics/lensing/pair_bin.h"
#include "brg/physics/units/unit_obj.h"
#include "brg/vector/summary_functions.hpp"

#include "pair_bin_summary.h"

namespace brgastro {


pair_bin_summary::pair_bin_summary( CONST_BRG_DISTANCE_REF init_R_min, CONST_BRG_DISTANCE_REF init_R_max,
		CONST_BRG_MASS_REF init_m_min, CONST_BRG_MASS_REF init_m_max,
		double init_z_min, double init_z_max,
		double init_mag_min, double init_mag_max )
:	_R_min_(init_R_min),
	_R_max_(init_R_max),
	_m_min_(init_m_min),
	_m_max_(init_m_max),
	_z_min_(init_z_min),
	_z_max_(init_z_max),
	_mag_min_(init_mag_min),
	_mag_max_(init_mag_max),
	_count_(0),
	_R_mean_(0),
	_m_mean_(0),
	_z_mean_(0),
	_mag_mean_(0),
	_delta_Sigma_t_mean_(0),
	_delta_Sigma_x_mean_(0),
	_delta_Sigma_t_mean_square_(0),
	_delta_Sigma_x_mean_square_(0)
{
}

pair_bin_summary::pair_bin_summary( const pair_bin & bin)
: 	_R_min_(bin.R_min()),
	_R_max_(bin.R_max()),
	_m_min_(bin.m_min()),
	_m_max_(bin.m_max()),
	_z_min_(bin.z_min()),
	_z_max_(bin.z_max()),
	_mag_min_(bin.mag_min()),
	_mag_max_(bin.mag_min()),
	_count_(bin.count()),
	_R_mean_(bin.R_mean()),
	_m_mean_(bin.m_mean()),
	_z_mean_(bin.z_mean()),
	_mag_mean_(bin.mag_mean()),
	_delta_Sigma_t_mean_(bin.delta_Sigma_t_mean()),
	_delta_Sigma_x_mean_(bin.delta_Sigma_x_mean()),
	_delta_Sigma_t_mean_square_(bin.delta_Sigma_t_mean_square()),
	_delta_Sigma_x_mean_square_(bin.delta_Sigma_x_mean_square())
{
}

// Adding and clearing data
#if(1)

void pair_bin_summary::clear()
{
	*this = pair_bin_summary(_R_min_,_R_max_,_m_min_,_m_max_,_z_min_,_z_max_,_mag_min_,_mag_max_);
}

#endif // Adding and clearing data

// Combining summaries together
#if(1)

pair_bin_summary & pair_bin_summary::operator+=( const pair_bin_summary & other )
{
	// Check bin limits are all the same
	assert(_R_min_==other.R_min());
	assert(_R_max_==other.R_max());
	assert(_m_min_==other.m_min());
	assert(_m_max_==other.m_max());
	assert(_z_min_==other.z_min());
	assert(_z_max_==other.z_max());
	assert(_mag_min_==other.mag_min());
	assert(_mag_max_==other.mag_max());

	size_t new_count = _count_ + other.count();

	_R_mean_ = ( _R_mean_*count() + other.R_mean()*other.count())/new_count;
	_m_mean_ = ( _m_mean_*count() + other.m_mean()*other.count())/new_count;
	_z_mean_ = ( _z_mean_*count() + other.z_mean()*other.count())/new_count;
	_mag_mean_ = ( _mag_mean_*count() + other.mag_mean()*other.count())/new_count;

	_delta_Sigma_t_mean_ = ( _delta_Sigma_t_mean_*count() + other.delta_Sigma_t_mean()*other.count())
			/new_count;
	_delta_Sigma_t_mean_square_ = ( _delta_Sigma_t_mean_square_*count() +
			other.delta_Sigma_t_mean_square()*other.count())
			/new_count;
	_delta_Sigma_x_mean_ = ( _delta_Sigma_t_mean_square_*count() +
			other.delta_Sigma_x_mean_square()*other.count())
			/new_count;

	_count_ = new_count;

	return *this;
}

#endif // Combining summaries together

} // end namespace brgastro
