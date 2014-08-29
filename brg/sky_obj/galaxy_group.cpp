/**       @file galaxygalaxy_group.cpp
 *
 *     Project: brg
 *        Path: /brg/sky_obj/galaxygalaxy_group.cpp
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#include <cstdio>
#include <stdexcept>
#include <vector>

#include "brg/brg_global.h"

#include "brg/sky_obj/galaxy.h"

#include "galaxy_group.h"

// brgastro::galaxy_group class methods
#if (1)

// Constructor
brgastro::galaxy_group::galaxy_group( size_t init_num_members )
{
	if ( init_num_members < 0 )
		init_num_members = 0;

	z_phot = z_phot_err = 0;
	odds = 1;

	BCG_index = 0;

	num_members = init_num_members;
	if ( num_members > 0 )
	{
		member_indices.reserve( num_members );
		members.reserve( num_members );
	}
	for ( size_t i = 0; i < num_members; i++ )
	{
		member_indices[i] = -1;
		members[i] = 0;
	}
}

// Destructor
brgastro::galaxy_group::~galaxy_group()
{
}

// Other functions
void brgastro::galaxy_group::clear()
{
	z_phot = z_phot_err = 0;
	odds = 1;

	BCG_index = 0;

	num_members = 0;
	if ( num_members > 0 )
	{
		member_indices.reserve( num_members );
		members.reserve( num_members );
	}
	for ( size_t i = 0; i < num_members; i++ )
	{
		member_indices[i] = -1;
		members[i] = 0;
	}

	sky_obj::clear();
}

void brgastro::galaxy_group::resize( size_t new_num_members )
{
	member_indices.resize( new_num_members, -1 );
	members.resize( new_num_members, 0 );
	num_members = new_num_members;
}

void brgastro::galaxy_group::set_member( size_t member_index, galaxy * new_member,
		const bool silent )
{
	if ( member_index >= num_members )
	{
		throw std::out_of_range("Group member index out of range.");
	}
	members[member_index] = new_member;
	member_indices[member_index] = new_member->index();
}
void brgastro::galaxy_group::set_member_index( size_t member_index,
		size_t new_member_index, const bool silent )
{
	if ( member_index >= num_members )
	{
		throw std::out_of_range("ERROR: Member index out of bounds in galaxy_group::set_member_index\n");
	}
	member_indices[member_index] = new_member_index;
}
void brgastro::galaxy_group::add_member( galaxy * new_member, const bool silent )
{
	int new_num_members = num_members + 1;
	resize( new_num_members );
	members[new_num_members - 1] = new_member;
	member_indices[new_num_members - 1] = new_member->index();
}
void brgastro::galaxy_group::remove_member( galaxy * rem_member,
		const bool silent )
{
	size_t i;
	for ( i = 0; i < num_members; i++ )
	{
		if ( members[i]->ID() == rem_member->ID() )
			break;
	}
	if ( i == num_members ) // Then we didn't find it
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Could not find member to remove from galaxy_group.\n";
		return;
	}
	for ( size_t j = i; j < num_members - 1; j++ )
	{
		members[j] = members[j + 1];
		member_indices[j] = member_indices[j + 1];
	}
	resize( num_members - 1 );
}

#endif // End brgastro::galaxy_group class functions
