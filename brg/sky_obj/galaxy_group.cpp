/**       @file galaxygalaxy_group.cpp
 *
 *     Project: brg
 *        Path: /brg/sky_obj/galaxygalaxy_group.cpp
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#include <cstdio>
#include <vector>

#include "../brg_global.h"

#include "galaxy_group.h"
#include "galaxy.h"

// brgastro::galaxy_group class methods
#if (1)

// Constructor
brgastro::galaxy_group::galaxy_group( int init_num_members )
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
	for ( int i = 0; i < num_members; i++ )
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
const int brgastro::galaxy_group::clear()
{
	galaxy_group();
	return 0;
}

const int brgastro::galaxy_group::resize( int new_num_members )
{
	member_indices.resize( new_num_members, -1 );
	members.resize( new_num_members, 0 );
	num_members = new_num_members;
	return 0;
}

const int brgastro::galaxy_group::set_member( int member_index, galaxy * new_member,
		const bool silent )
{
	if ( member_index >= num_members )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Member index out of bounds in galaxy_group::set_member\n";
		return OUT_OF_BOUNDS_ERROR;
	}
	members[member_index] = new_member;
	member_indices[member_index] = new_member->index();
	return 0;
}
const int brgastro::galaxy_group::set_member_index( int member_index,
		int new_member_index, const bool silent )
{
	if ( member_index >= num_members )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Member index out of bounds in galaxy_group::set_member_index\n";
		return OUT_OF_BOUNDS_ERROR;
	}
	member_indices[member_index] = new_member_index;
	return 0;
}
const int brgastro::galaxy_group::add_member( galaxy * new_member, const bool silent )
{
	int new_num_members = num_members + 1;
	resize( new_num_members );
	members[new_num_members - 1] = new_member;
	member_indices[new_num_members - 1] = new_member->index();
	return 0;
}
const int brgastro::galaxy_group::remove_member( galaxy * rem_member,
		const bool silent )
{
	int i;
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
		return OUT_OF_BOUNDS_ERROR;
	}
	for ( int j = i; j < num_members - 1; j++ )
	{
		members[j] = members[j + 1];
		member_indices[j] = member_indices[j + 1];
	}
	resize( num_members - 1 );
	return 0;
}

#endif // End brgastro::galaxy_group class functions
