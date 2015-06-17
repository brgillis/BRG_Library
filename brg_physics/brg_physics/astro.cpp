/**********************************************************************\
  @file astro.cpp

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

#include <cstdlib>
#include <iostream>
#include <string>

#include "brg/common.h"

#include "brg/math/cache/cache.hpp"
#include "brg/math/cache/cache_2d.hpp"

#include "brg_physics/astro_caches.h"
#include "brg_physics/astro.h"
#include "brg_physics/sky_obj/position_grid_cache.h"
#include "brg_physics/units/unit_obj.h"

// brgastro::redshift_obj class methods
#if (1)

BRG_UNITS brgastro::redshift_obj::H() const
{
	// If not cached, calculate and cache it
	if(!_H_cached_)
	{
		if(_z_==0)
		{
			_H_cache_ = H_0;
		}
		else
		{
			_H_cache_ = brgastro::H(_z_);
		}
		_H_cached_ = true;
	}
	return _H_cache_;
}

BRG_UNITS brgastro::redshift_obj::rho_crit() const
{
	return 3*square(H())/(8*pi*Gc);
}

#endif // brgastro::redshift_obj class methods

/** Global Function Definitions **/
#if (1)

BRG_UNITS brgastro::H( const flt_type test_z )
{
	// Friedmann equation, assuming omega = -1
	if(test_z==0) return H_0;
	flt_type zp1 = 1.+test_z;
	return H_0
			* std::sqrt( Omega_r * quart( zp1 )
							+ Omega_m * cube( zp1 )
							+ Omega_k * square( zp1 ) + Omega_l );
}

// grid functions
#if (1)

#endif // end grid functions

// dfa and afd functions
#if (1)

BRG_DISTANCE brgastro::dfa( const BRG_ANGLE & da, const flt_type z )
{
	return da * dfa_cache().get( z );
}
BRG_DISTANCE brgastro::dfa( const BRG_ANGLE & a1, const BRG_ANGLE & a2,
		const flt_type z )
{
	return brgastro::dfa( a2 - a1, z );
}
BRG_DISTANCE brgastro::dfa( const BRG_ANGLE & a1x, const BRG_ANGLE & a1y,
		const BRG_ANGLE & a2x, const BRG_ANGLE & a2y, const flt_type z )
{
	return brgastro::dfa( skydist2d( a1x, a1y, a2x, a2y ), z );
}
BRG_ANGLE brgastro::afd( const BRG_DISTANCE & dd, const flt_type z )
{
	return dd / safe_d( dfa_cache().get( z ) );
}
BRG_ANGLE brgastro::afd( const BRG_DISTANCE & d1, const BRG_DISTANCE & d2,
		const flt_type z )
{
	return brgastro::afd( fabs( d2 - d1 ), z );
}
BRG_ANGLE brgastro::afd( const BRG_DISTANCE & d1x,
		const BRG_DISTANCE & d1y, const BRG_DISTANCE & d2x,
		const BRG_DISTANCE & d2y, const flt_type z )
{
	return brgastro::afd( dist2d( d1x, d1y, d2x, d2y ), z );
}

flt_type brgastro::zfa( const flt_type a )
{
	return 1. / safe_d( a ) - 1.;
}
flt_type brgastro::afz( const flt_type z )
{
	return 1. / safe_d( 1 + z );
}

BRG_TIME brgastro::tfz( const flt_type z )
{
	return brgastro::tfa( afz( z ) );
}
BRG_TIME brgastro::tfa( const flt_type a )
{
	return tfa_cache().get( a );
}
flt_type brgastro::zft( const BRG_TIME & t )
{
	return brgastro::zfa( brgastro::aft( t ) );
}
flt_type brgastro::aft( const BRG_TIME & t )
{
	brgastro::tfa_cache cache;
	return cache.inverse_get( t );
}

#endif

// Functions to integrate out distances
#if (1)
flt_type brgastro::integrate_add( const flt_type z1, const flt_type z2 )
{
	return brgastro::integrate_distance( z1, z2, 0, 100000 );
}
flt_type brgastro::integrate_cmd( const flt_type z1, const flt_type z2 )
{
	return brgastro::integrate_distance( z1, z2, 1, 10000000 );
}
flt_type brgastro::integrate_Ld( const flt_type z1, const flt_type z2 )
{
	return brgastro::integrate_distance( z1, z2, 2, 10000000 );
}
flt_type brgastro::integrate_ltd( const flt_type z1, const flt_type z2 )
{
	return brgastro::integrate_distance( z1, z2, 3, 10000000 );
}
flt_type brgastro::integrate_add( const flt_type z )
{
	return brgastro::integrate_distance( 0, z, 0, 100000 );
}
flt_type brgastro::integrate_cmd( const flt_type z )
{
	return brgastro::integrate_distance( 0, z, 1, 10000000 );
}
flt_type brgastro::integrate_Ld( const flt_type z )
{
	return brgastro::integrate_distance( 0, z, 2, 10000000 );
}
flt_type brgastro::integrate_ltd( const flt_type z )
{
	return brgastro::integrate_distance( 0, z, 3, 10000000 );
}
flt_type brgastro::integrate_distance( const flt_type z1_init,
		const flt_type z2_init, const int_type mode, const int_type n )
{
	// Function that will help calculate cosmic distances, thanks to Richard Powell - http://www.atlasoftheuniverse.com/
	// NOTE: Will return a negative value if z2 < z1. This is by design.

	flt_type OK0;
	flt_type OM0 = Omega_m, OR0 = Omega_r, OL0 = Omega_l;
	flt_type OM, OR, OL, OK;
	flt_type z1 = z1_init, z2 = z2_init;
	flt_type HD; //Hubble distance in billions of lightyears
	flt_type z, a, a1, a2, adot, h1;
	flt_type DC = std::numeric_limits<flt_type>::max(), DCC = 0, DT = std::numeric_limits<flt_type>::max(), DTT = 0, DM;
	//flt_type age, size;
	int_type i;
	short_int_type sign = 1;

	if(z1==z2) return 0;

	OK0 = 1 - OM0 - OL0 - OR0;

	if ( z2 < z1 )
	{
		std::swap(z1,z2);
		sign = -1;
	}

	if ( z1 == 0 )
	{
		z = z2;
		h1 = H_0;
		OM = OM0;
		OR = OR0;
		OL = OL0;
		OK = OK0;
	}
	else
	{
		a1 = 1 / ( 1 + z1 );
		a2 = 1 / ( 1 + z2 );
		z = ( a1 / a2 ) - 1;
		h1 = H_0
				* std::sqrt( OR0 * inv_quart( a1 ) + OM0 * inv_cube( a1 )
								+ OK0 * inv_square( a1 ) + OL0 );
		OM = OM0 * square( H_0 / h1 ) * inv_cube( a1 );
		OR = OR0 * square( H_0 / h1 ) * inv_quart( a1 );
		OL = OL0 * square( H_0 / h1 );
		OK = 1 - OM - OR - OL;
	}

	HD = c / h1;

	for ( i = n; i >= 1; i-- )        // This loop is the numerical integration
	{
		a = ( i - 0.5 ) / n;              // Steadily decrease the scale factor
		// Comoving formula (See section 4 of Hogg, but I've added a radiation term too):
		adot = a * std::sqrt( OM * inv_cube(a) + OK * inv_square(a)
								+ OL + OR * inv_quart(a) ); // Note that "a" is equivalent to 1/(1+z)
		 // Collect DC and DT until the correct scale is reached
		DCC = DCC + 1 / ( a * adot ) / n; // Running total of the comoving distance
		DTT = DTT + 1 / adot / n; // Running total of the light travel time (see section 10 of Hogg)
		 // Collect DC and DT until the correct scale is reached
		DC = DCC;                 // Comoving distance DC
		DT = DTT;                 // Light travel time DT
		if ( a <= 1 / ( 1 + z ) ) // Collect DC and DT until the correct scale is reached
		{
			break;
		}
	}

	// Transverse comoving distance DM from section 5 of Hogg:
	if ( OK > 0.0001 )
		DM = ( 1 / std::sqrt( OK ) ) * sinh( std::sqrt( OK ) * DC );
	else if ( OK < -0.0001 )
		DM = ( 1 / std::sqrt( fabs( OK ) ) ) * sin( std::sqrt( fabs( OK ) ) * DC );
	else
		DM = DC;

	//age = HD*DTT;                 // Age of the universe (billions of years)
	//size = HD*DCC;                // Comoving radius of the observable universe

	switch ( mode )
	{
	case 0:
		return sign * HD * DM / ( 1 + z ); // Angular diameter distance (section 6 of Hogg)
		break;
	case 1:
		return sign * HD * DC;             // Comoving distance
		break;
	case 2:
		return sign * HD * DM * ( 1 + z );  // Luminosity distance (section 7 of Hogg)
		break;
	case 3:
		return sign * HD * DT;             // Light travel distance
		break;
	default:
		return sign * HD * DT;             // Light travel distance
	}
}
#endif

// Other functions
#if (1)

BRG_DISTANCE brgastro::ad_distance( flt_type z1, flt_type z2 )
{
	if ( z2 < z1 )
		std::swap( z1, z2 );
	return brgastro::add_cache().get( z1, z2 );
}

BRG_UNITS brgastro::sigma_crit( const flt_type z_lens,
		const flt_type z_source )
{
	return square( c ) / ( 4. * pi * Gc )
			* brgastro::ad_distance( 0, z_source )
			/ ( brgastro::ad_distance( 0, z_lens )
					* brgastro::ad_distance( z_lens, z_source ) );

}
#endif // end other functions

#endif // end Global function definitions
