/**********************************************************************\
  @file galaxy.h

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

// body file: galaxy.cpp

#ifndef _BRG_GALAXY_H_
#define _BRG_GALAXY_H_

#include "brg/global.h"

#include "brg_physics/sky_obj/sky_obj.h"

namespace brgastro {

// Forward declare galaxy_group class so we can point to it
class galaxy_group;

class galaxy: public virtual sky_obj
{
	/**********************************
	 galaxy class
	 -------------

	 A class for galaxies.

	 Parent class: sky_obj

	 **********************************/

private:
	// private member variables (none needed in current implementation. Should make all private for consistency in future though)

public:

	// Public member variables

	BRG_MASS stellar_mass;
	double umag, umag_err;
	double gmag, gmag_err;
	double rmag, rmag_err;
	double imag, imag_err;
	double zmag, zmag_err;

	double z_phot, z_phot_err;
	double odds;
	double phot_template;

	galaxy_group *host_group; // Pointer to group this galaxy resides within
	int host_group_index;

	// Public member functions

	// Constructors
	galaxy( CONST_BRG_ANGLE_REF init_ra = 0, CONST_BRG_ANGLE_REF init_dec = 0, double init_z = 0,
			CONST_BRG_ANGLE_REF init_ra_err = 0, CONST_BRG_ANGLE_REF init_dec_err = 0, double init_z_err =
			0, CONST_BRG_MASS_REF init_stellar_mass=0,
			double init_mag=0, double init_mag_err=0 );

	// Copy constructor
	//galaxy(const galaxy other_galaxy); // Implicit is fine for us

	// Virtual destructor
	virtual ~galaxy()
	{
	}

	// Implementations of pure virtual functions from sky_obj
#if (1)
	BRG_MASS m() const
	{
		return stellar_mass;
	}
	double mag() const
	{
		return imag;
	}
#endif

	// Clear function
	virtual void clear();

	// Clone functions
	virtual redshift_obj *redshift_obj_clone() const
	{
		return new galaxy( *this );
	}
	virtual sky_obj *sky_obj_clone() const
	{
		return new galaxy( *this );
	}
	virtual galaxy *galaxy_clone() const
	{
		return new galaxy( *this );
	}

};
// class galaxy

} // end namespace brgastro

#endif /* _BRG_GALAXY_H_ */
