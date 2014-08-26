/**********************************************************************\
brg_phase.hpp
 -----------



 Everything in this file is declared in the namespace brgastro.

 \**********************************************************************/

#ifndef __BRG_PHASE_HPP_INCLUDED__
#define __BRG_PHASE_HPP_INCLUDED__

#include "brg_global.h"

#ifdef _BRG_USE_UNITS_
#include "brg_units.h"
#endif

namespace brgastro
{

/** Class definitions **/
#if (1)
struct phase
/**********************************************
 phase
 -----

 A structure representing the full phase of an
 object (position and velocity) and also time,
 to help limit the number of variables that need
 to be passed to certain functions.

 All member variables are public, and may be
 accessed directly. Units should be in the
 default set: m for position, m/s for velocity,
 s for time.

 **********************************************/
{
	BRG_DISTANCE x, y, z;

	BRG_VELOCITY vx, vy, vz;

	BRG_TIME t;
	phase( CONST_BRG_DISTANCE_REF init_x = 0, CONST_BRG_DISTANCE_REF init_y = 0,
			CONST_BRG_DISTANCE_REF init_z = 0, CONST_BRG_VELOCITY_REF init_vx = 0,
			CONST_BRG_VELOCITY_REF init_vy = 0, CONST_BRG_VELOCITY_REF init_vz = 0,
	BRG_TIME init_t = 0 )
	{
		x = init_x;
		y = init_y;
		z = init_z;
		vx = init_vx;
		vy = init_vy;
		vz = init_vz;
		t = init_t;
	}
	const int set_phase( CONST_BRG_DISTANCE_REF init_x = 0, CONST_BRG_DISTANCE_REF init_y = 0,
			CONST_BRG_DISTANCE_REF init_z = 0, CONST_BRG_VELOCITY_REF init_vx = 0,
			CONST_BRG_VELOCITY_REF init_vy = 0, CONST_BRG_VELOCITY_REF init_vz = 0,
	BRG_TIME init_t = 0 )
	{
		x = init_x;
		y = init_y;
		z = init_z;
		vx = init_vx;
		vy = init_vy;
		vz = init_vz;
		t = init_t;

		return 0;
	}
};

#endif // end class declarations

}

#endif // __BRG_PHASE_HPP_INCLUDED__
