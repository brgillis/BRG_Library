/**********************************************************************\
 position_grid_cache.h
 -----------

 If this header is used, the source file position_grid_cache.h must be
 included and compiled with the project.

 \**********************************************************************/

#ifndef _BRG_POSITION_GRID_CACHE_H_INCLUDED_
#define _BRG_POSITION_GRID_CACHE_H_INCLUDED_

#include "brg/global.h"

#include "brg/physics/units/unit_obj.h"

namespace brgastro
{

class grid_cache
{
private:
	static int _ra_grid_change_num_, _dec_grid_change_num_,
			_z_grid_change_num_;
	static BRG_ANGLE _ra_grid_min_, _ra_grid_max_, _ra_grid_step_;
	static BRG_ANGLE _dec_grid_min_, dec_grid_max_val, _dec_grid_step_;
	static double _z_grid_min_, _z_grid_max_, _z_grid_step_;
public:
	// Set functions
#if (1)
	const int set_ra_grid( CONST_BRG_ANGLE_REF new_ra_grid_min,
			CONST_BRG_ANGLE_REF new_ra_grid_max, CONST_BRG_ANGLE_REF new_ra_grid_step )
	{
		_ra_grid_min_ = new_ra_grid_min;
		_ra_grid_max_ = new_ra_grid_max;
		_ra_grid_step_ = new_ra_grid_step;
		_ra_grid_change_num_++;
		return 0;
	}

	const int set_dec_grid( CONST_BRG_ANGLE_REF new_dec_grid_min,
			CONST_BRG_ANGLE_REF new_dec_grid_max,
			CONST_BRG_ANGLE_REF new_dec_grid_step )
	{
		_dec_grid_min_ = new_dec_grid_min;
		dec_grid_max_val = new_dec_grid_max;
		_dec_grid_step_ = new_dec_grid_step;
		_dec_grid_change_num_++;
		return 0;
	}

	const int set_z_grid( const double new_z_grid_min,
			const double new_z_grid_max, const double new_z_grid_step )
	{
		_z_grid_min_ = new_z_grid_min;
		_z_grid_max_ = new_z_grid_max;
		_z_grid_step_ = new_z_grid_step;
		_z_grid_change_num_++;
		return 0;
	}
#endif

	// Get functions
#if (1)

	const int ra_grid_change_num()
	{
		return _ra_grid_change_num_;
	}
	const int dec_grid_change_num()
	{
		return _dec_grid_change_num_;
	}
	const int z_grid_change_num()
	{
		return _z_grid_change_num_;
	}
	CONST_BRG_ANGLE_REF ra_grid_min()
	{
		return _ra_grid_min_;
	}
	CONST_BRG_ANGLE_REF dec_grid_min()
	{
		return _dec_grid_min_;
	}
	const double z_grid_min()
	{
		return _z_grid_min_;
	}
	CONST_BRG_ANGLE_REF ra_grid_max()
	{
		return _ra_grid_max_;
	}
	CONST_BRG_ANGLE_REF dec_grid_max()
	{
		return dec_grid_max_val;
	}
	const double z_grid_max()
	{
		return _z_grid_max_;
	}
	CONST_BRG_ANGLE_REF ra_grid_step()
	{
		return _ra_grid_step_;
	}
	CONST_BRG_ANGLE_REF dec_grid_step()
	{
		return _dec_grid_step_;
	}
	const double z_grid_step()
	{
		return _z_grid_step_;
	}
#endif

};
// class grid_cache

} // end namespace brgastro

#endif // __BRG_POSITION_GRID_CACHE_H_INCLUDED__

