/**********************************************************************\
  @file sky_obj.cpp

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

#include "brg/common.h"

#include "brg/units/units.hpp"

#include "brg_physics/astro_caches.h"
#include "brg_physics/astro.h"

#include "sky_obj.h"

// brgastro::sky_obj class methods
#if (1)

brgastro::sky_obj::sky_obj( const angle_type & init_ra, const angle_type & init_dec,
		const flt_type & init_z, const angle_type & init_ra_err,
		const angle_type & init_dec_err, const flt_type & init_z_err )
{
	partial_clear();
	set_ra_dec_z_err( init_ra, init_dec, init_z, init_ra_err, init_dec_err,
			init_z_err );
}

void brgastro::sky_obj::clear()
{
	set_ra_dec_z_err( 0, 0, 0, 0, 0, 0 );
	return partial_clear();
}

void brgastro::sky_obj::partial_clear()
{
	_index_ = 0;
	_ID_ = "0";
	_weight_ = 1;
}

void brgastro::sky_obj::set_ra( const angle_type & new_ra )
{
	_ra_ = new_ra;
}
void brgastro::sky_obj::set_ra_err( const angle_type & new_ra_err )
{
	_ra_err_ = new_ra_err;
}
void brgastro::sky_obj::set_dec( const angle_type & new_dec )
{
	_dec_ = new_dec;
}
void brgastro::sky_obj::set_dec_err( const angle_type & new_dec_err )
{
	_dec_err_ = new_dec_err;
}
void brgastro::sky_obj::set_ra_dec( const angle_type & new_ra,
		const angle_type & new_dec )
{
	set_ra( new_ra );
	set_dec( new_dec );
}
void brgastro::sky_obj::set_ra_dec_z( const angle_type & new_ra,
		const angle_type & new_dec, const flt_type & new_z )
{
	set_ra( new_ra );
	set_dec( new_dec );
	set_z( new_z );
}
void brgastro::sky_obj::set_ra_dec_z_err( const angle_type & new_ra,
		const angle_type & new_dec, const flt_type & new_z,
		const angle_type & new_ra_err, const angle_type & new_dec_err,
		const flt_type & new_z_err )
{
	set_ra( new_ra );
	set_dec( new_dec );
	set_z( new_z );
	set_ra_err( new_ra_err );
	set_dec_err( new_dec_err );
	set_z_err( new_z_err );
}
void brgastro::sky_obj::set_ra_dec_err( const angle_type & new_ra,
		const angle_type & new_dec, const angle_type & new_ra_err,
		const angle_type & new_dec_err )
{
	set_ra( new_ra );
	set_dec( new_dec );
	set_ra_err( new_ra_err );
	set_dec_err( new_dec_err );
}
void brgastro::sky_obj::set_weight( const flt_type & new_weight )
{
	_weight_ = new_weight;
}
void brgastro::sky_obj::set_index( const int_type new_index )
{
	_index_ = new_index;
}
void brgastro::sky_obj::set_ID( const std::string &new_ID )
{
	_ID_ = new_ID;
}

#endif // end brgastro::sky_obj class methods

brgastro::distance_type brgastro::dfa( const brgastro::sky_obj *obj1,
		const brgastro::sky_obj *obj2, const flt_type & z )
{
	flt_type z_to_use;
	if ( z == -1 )
		z_to_use = obj1->z();
	else
		z_to_use = z;
	return brgastro::dfa(
			skydist2d( obj1->ra(), obj1->dec(), obj2->ra(), obj2->dec() ),
			z_to_use );
}

brgastro::angle_type brgastro::skydist2d( const sky_obj *obj1, const sky_obj *obj2 )
{
	return skydist2d(obj1->ra(),obj1->dec(),obj2->ra(),obj2->dec());
}
