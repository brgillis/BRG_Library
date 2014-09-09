/**********************************************************************\
  @file galaxy_group.h

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

// body file: galaxy_group.cpp

#ifndef _BRG_GALAXY_GROUP_H_
#define _BRG_GALAXY_GROUP_H_

#include <vector>

#include "brg/global.h"

#include "brg/physics/sky_obj/galaxy.h"
#include "brg/physics/sky_obj/sky_obj.h"

namespace brgastro {


class galaxy_group: public sky_obj
{
	/**********************************
	 group class
	 -------------

	 A class for groups and clusters of galaxies.

	 Parent class: sky_obj

	 Derived classes:
	 (none)

	 **********************************/

private:

public:
	// Public member variables

	double z_phot, z_phot_err;
	double odds;

	size_t BCG_index;

	size_t num_members;
	std::vector< size_t > member_indices;
	std::vector< galaxy * > members;

	// Public member functions

	// Constructor
	galaxy_group( size_t init_num_members = 0 );
	galaxy_group( double init_mass, double init_z, double init_c = -1,
			double init_tau = -1 );

	// Copy constructor
	//group( const group &other ); // Implicit is fine

	// Destructor
	virtual ~galaxy_group();

	// Functions to set parameters
	virtual void clear();
	virtual void resize( size_t new_num_members );
	virtual void set_member( size_t index, galaxy * new_member,
			const bool silent = false );
	virtual void set_member_index( size_t index, size_t new_member_index,
			const bool silent = false );
	virtual void add_member( galaxy * new_member, const bool silent =
			false );
	virtual void remove_member( galaxy * rem_member, const bool silent =
			false );

	// Implementations of pure virtual functions from sky_obj
#if (1)
	BRG_MASS m() const
	{
		return 0; // TODO: Should be combine magnitude of member galaxies
	}
	double mag() const
	{
		return 0; // TODO: Should be combine magnitude of member galaxies
	}
#endif

	// Clone functions
	virtual redshift_obj *redshift_obj_clone() const
	{
		return new galaxy_group( *this );
	}
	virtual sky_obj *sky_obj_clone() const
	{
		return new galaxy_group( *this );
	}
	virtual galaxy_group *galaxy_group_clone() const
	{
		return new galaxy_group( *this );
	}

};
// class galaxy_group

} // end namespace brgastro

#endif /* _BRG_GALAXYGROUP_H_ */
