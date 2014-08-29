/**       @file galaxygroup.h
 *
 *     Project: brg
 *        Path: /brg/sky_obj/galaxygroup.h
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#ifndef _BRG_GALAXY_GROUP_H_
#define _BRG_GALAXY_GROUP_H_

#include <vector>

#include "brg/global.h"

#include "brg/physics/astro/sky_obj/galaxy.h"
#include "brg/physics/astro/sky_obj/sky_obj.h"

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

	// Clone functions
	virtual redshift_obj *redshift_obj_clone()
	{
		return new galaxy_group( *this );
	}
	virtual sky_obj *sky_obj_clone()
	{
		return new galaxy_group( *this );
	}
	virtual galaxy_group *galaxy_group_clone()
	{
		return new galaxy_group( *this );
	}

};
// class galaxy_group

} // end namespace brgastro

#endif /* _BRG_GALAXYGROUP_H_ */
