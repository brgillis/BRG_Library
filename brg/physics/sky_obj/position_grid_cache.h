/**********************************************************************\
  @file position_grid_cache.h

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

// body file: position_grid_cache.cpp

#ifndef _BRG_POSITION_GRID_CACHE_H_INCLUDED_
#define _BRG_POSITION_GRID_CACHE_H_INCLUDED_

#include <cmath>

#include "brg/global.h"

#include "brg/math/safe_math.hpp"
#include "brg/physics/units/unit_obj.h"

namespace brgastro
{

class grid_cache
{
private:
	static unsigned int _ra_grid_change_num_, _dec_grid_change_num_,
			_z_grid_change_num_;
	static BRG_ANGLE _ra_grid_min_, _ra_grid_max_, _ra_grid_step_;
	static BRG_ANGLE _dec_grid_min_, _dec_grid_max_, _dec_grid_step_;
	static double _z_grid_min_, _z_grid_max_, _z_grid_step_;
public:
	// Set functions
#if (1)
	static void set_ra_grid( CONST_BRG_ANGLE_REF new_ra_grid_min,
			CONST_BRG_ANGLE_REF new_ra_grid_max, CONST_BRG_ANGLE_REF new_ra_grid_step )
	{
		_ra_grid_min_ = new_ra_grid_min;
		_ra_grid_max_ = new_ra_grid_max;
		_ra_grid_step_ = new_ra_grid_step;
		_ra_grid_change_num_++;
	}

	static void set_dec_grid( CONST_BRG_ANGLE_REF new_dec_grid_min,
			CONST_BRG_ANGLE_REF new_dec_grid_max,
			CONST_BRG_ANGLE_REF new_dec_grid_step )
	{
		_dec_grid_min_ = new_dec_grid_min;
		_dec_grid_max_ = new_dec_grid_max;
		_dec_grid_step_ = new_dec_grid_step;
		_dec_grid_change_num_++;
	}

	static void set_z_grid( const double new_z_grid_min,
			const double new_z_grid_max, const double new_z_grid_step )
	{
		_z_grid_min_ = new_z_grid_min;
		_z_grid_max_ = new_z_grid_max;
		_z_grid_step_ = new_z_grid_step;
		_z_grid_change_num_++;
	}
#endif

	// Get functions
#if (1)

	static const int ra_grid_change_num()
	{
		return _ra_grid_change_num_;
	}
	static const int dec_grid_change_num()
	{
		return _dec_grid_change_num_;
	}
	static const int z_grid_change_num()
	{
		return _z_grid_change_num_;
	}
	static CONST_BRG_ANGLE_REF ra_grid_min()
	{
		return _ra_grid_min_;
	}
	static CONST_BRG_ANGLE_REF dec_grid_min()
	{
		return _dec_grid_min_;
	}
	static const double z_grid_min()
	{
		return _z_grid_min_;
	}
	static CONST_BRG_ANGLE_REF ra_grid_max()
	{
		return _ra_grid_max_;
	}
	static CONST_BRG_ANGLE_REF dec_grid_max()
	{
		return _dec_grid_max_;
	}
	static const double z_grid_max()
	{
		return _z_grid_max_;
	}
	static CONST_BRG_ANGLE_REF ra_grid_step()
	{
		return _ra_grid_step_;
	}
	static CONST_BRG_ANGLE_REF dec_grid_step()
	{
		return _dec_grid_step_;
	}
	static const double z_grid_step()
	{
		return _z_grid_step_;
	}
#endif

	// Functions to get grid integers or grid boundaries from integers
#if(1)
	static int get_ra_grid( CONST_BRG_ANGLE_REF ra )
	{
		return (int)floor(
				( ra - _ra_grid_min_ ) / safe_d( _ra_grid_step_ ) );
	}
	static int get_dec_grid( CONST_BRG_ANGLE_REF dec )
	{
		return (int)floor( ( dec - _dec_grid_min_ ) / _dec_grid_step_ );
	}
	static int get_z_grid( const double z )
	{
		return (int)floor( ( z - _z_grid_min_ ) / _z_grid_step_ );
	}

	static BRG_ANGLE get_ra_grid_lower( const int ra_grid )
	{
		return _ra_grid_min_ + _ra_grid_step_ * ra_grid;
	}
	static BRG_ANGLE get_dec_grid_lower( const int dec_grid )
	{
		return _dec_grid_min_ + _dec_grid_step_ * dec_grid;
	}
	static double get_z_grid_lower( const int z_grid )
	{
		return _z_grid_min_ + _z_grid_step_ * z_grid;
	}

	static BRG_ANGLE get_ra_grid_upper( const int ra_grid )
	{
		return _ra_grid_min_ + _ra_grid_step_ * ( ra_grid + 1 );
	}
	static BRG_ANGLE get_dec_grid_upper( const int dec_grid )
	{
		return _dec_grid_min_ + _dec_grid_step_ * ( dec_grid + 1 );
	}
	static double get_z_grid_upper( const int z_grid )
	{
		return _z_grid_min_ + _z_grid_step_ * ( z_grid + 1 );
	}

	static BRG_ANGLE get_ra_grid_mid( const int ra_grid )
	{
		return _ra_grid_min_ + _ra_grid_step_ * ( ra_grid + 0.5 );
	}
	static BRG_ANGLE get_dec_grid_mid( const int dec_grid )
	{
		return _dec_grid_min_ + _dec_grid_step_ * ( dec_grid + 0.5 );
	}
	static double get_z_grid_mid( const int z_grid )
	{
		return _z_grid_min_ + _z_grid_step_ * ( z_grid + 0.5 );
	}
#endif

};
// class grid_cache

} // end namespace brgastro

#endif // __BRG_POSITION_GRID_CACHE_H_INCLUDED__

