/**********************************************************************\
  @file sky_obj.h

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

// body file: sky_obj.cpp

#ifndef _BRG_SKY_OBJ_H_
#define _BRG_SKY_OBJ_H_

#include "brg/global.h"

#include "brg/physics/astro.h"
#include "brg/physics/sky_obj/position_grid_cache.h"

namespace brgastro {

class sky_obj: public virtual redshift_obj
{
	/**********************************
	 sky_obj class
	 -------------

	 An abstract class for anything in the sky.

	 Derived classes:
	 galaxy
	 group

	 **********************************/
private:
#if (1)
	std::string _ID_; // Name for it or ID number

	int _index_; // Position in an array

	double _weight_;

	BRG_ANGLE _ra_, _ra_err_, _dec_, _dec_err_;
#endif
public:
#if (1)
	// Public member functions

	// Constructor
	sky_obj( CONST_BRG_ANGLE_REF init_ra = 0, CONST_BRG_ANGLE_REF init_dec = 0, double init_z = 0,
			CONST_BRG_ANGLE_REF init_ra_err = 0, CONST_BRG_ANGLE_REF init_dec_err = 0, double init_z_err =
			0 ); // Normal constructor

	//sky_obj(const sky_obj & other_sky_obj); // Copy constructor (Default is fine for us)

	virtual ~sky_obj()
	{
	} // Virtual destructor, in case it needs to be overridden

	virtual void clear(); // Resets all variables to zero
	virtual void partial_clear(); // Resets all variables which can't be initialized

#if (1) // Set functions
	virtual void set_ra( CONST_BRG_ANGLE_REF new_ra );
	virtual void set_dec( CONST_BRG_ANGLE_REF new_dec );
	void set_ra_err( CONST_BRG_ANGLE_REF new_ra_err );
	void set_dec_err( CONST_BRG_ANGLE_REF new_dec_err );
	virtual void set_ra_dec( CONST_BRG_ANGLE_REF new_ra,
			CONST_BRG_ANGLE_REF new_dec ); // Sets ra and dec
	virtual void set_ra_dec_z( CONST_BRG_ANGLE_REF new_ra,
			CONST_BRG_ANGLE_REF new_dec, const double new_z ); // Sets all values
	virtual void set_ra_dec_z_err( CONST_BRG_ANGLE_REF new_ra,
			CONST_BRG_ANGLE_REF new_dec, const double new_z,
			CONST_BRG_ANGLE_REF new_ra_err, CONST_BRG_ANGLE_REF new_dec_err,
			const double new_z_err ); // Sets all values and error
	virtual void set_ra_dec_err( CONST_BRG_ANGLE_REF new_ra,
			CONST_BRG_ANGLE_REF new_dec, CONST_BRG_ANGLE_REF new_ra_err,
			CONST_BRG_ANGLE_REF new_dec_err ); // Sets ra and dec and error

	void set_weight( const double new_weight );
	void set_index( const int new_index );
	void set_ID( const std::string &new_ID );
#endif // end set functions

#if (1) //Get functions

	CONST_BRG_ANGLE_REF ra() const
	{
		return _ra_;
	}
	CONST_BRG_ANGLE_REF dec() const
	{
		return _dec_;
	}
	CONST_BRG_ANGLE_REF ra_err() const
	{
		return _ra_err_;
	}
	CONST_BRG_ANGLE_REF dec_err() const
	{
		return _dec_err_;
	}

	const double weight() const
	{
		return _weight_;
	}
	const int index() const
	{
		return _index_;
	}
	const std::string & ID() const
	{
		return _ID_;
	}

#endif // end get functions

	// Pure virtual get functions
#if (1)

	virtual BRG_MASS m() const = 0;
	virtual double mag() const = 0;

#endif

	// Clone function (to enable copies to be made of pointed-to objects)
	virtual redshift_obj *redshift_obj_clone() const =0;
	virtual sky_obj *sky_obj_clone() const =0;

#endif

};
// class sky_obj

BRG_DISTANCE dfa( const sky_obj *obj1, const sky_obj *obj2,
		const double z = -1 );

BRG_ANGLE skydist2d( const sky_obj *obj1, const sky_obj *obj2 );

} // end namespace brgastro

#endif /* _BRG_SKY_OBJ_H_ */
