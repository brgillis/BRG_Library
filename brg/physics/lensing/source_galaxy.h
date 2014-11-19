/**********************************************************************\
 @file source_galaxy.h
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

// body file: source_galaxy.cpp

#ifndef _BRG_SOURCE_GALAXY_H_INCLUDED_
#define _BRG_SOURCE_GALAXY_H_INCLUDED_

#include "brg/physics/sky_obj/galaxy.h"
#include "brg/physics/lensing/source_obj.hpp"

namespace brgastro {

/**
 *
 */
class source_galaxy: public source_obj, public galaxy {
public:

	// Constructors and destructors
#if(1)
	source_galaxy( CONST_BRG_ANGLE_REF init_ra = 0, CONST_BRG_ANGLE_REF init_dec = 0, const double & init_z = 0,
			const double & init_gamma_1 = 0, const double & init_gamma_2 = 0, const double & init_kappa = 0,
			CONST_BRG_MASS_REF init_mstar = 0, const double & init_mag = 0, const double & init_weight = 1)
	: sky_obj(init_ra,init_dec,init_z),
	  source_obj(init_ra, init_dec, init_z,
	  			init_gamma_1, init_gamma_2, init_kappa),
	  galaxy(init_ra,init_dec,init_z,0,0,0,init_mstar,init_mag)
	{
		set_weight(init_weight);
	}
	virtual ~source_galaxy() {
		// TODO Auto-generated destructor stub
	}
#endif

	// Implementations of pure virtual functions
#if(1)

	BRG_MASS m() const
	{
		return galaxy::m();
	}
	double mag() const
	{
		return galaxy::mag();
	}

	// Clone functions
	redshift_obj *redshift_obj_clone() const
	{
		return new source_galaxy( *this );
	}
	sky_obj *sky_obj_clone() const
	{
		return new source_galaxy( *this );
	}
	galaxy *galaxy_clone() const
	{
		return new source_galaxy( *this );
	}
	source_obj *source_obj_clone() const
	{
		return new source_galaxy( *this );
	}

#endif
};

} // end namespace brgastro

#endif // _BRG_SOURCE_GALAXY_H_INCLUDED_
