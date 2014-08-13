#include <cstdlib>
#include <iostream>
#include <string>

#include "brg_global.h"

#include "brg_units.h"
#include "brg_astro.h"
#include "cache/brg_cache.hpp"
#include "cache/brg_cache_2d.hpp"

using namespace std;

/** Static Class Initialisation **/
#if (1)
// Initialisation for brgastro::grid_cache
#if (1)
int brgastro::grid_cache::_ra_grid_change_num_ = 0;
int brgastro::grid_cache::_dec_grid_change_num_ = 0;
int brgastro::grid_cache::_z_grid_change_num_ = 0;
BRG_ANGLE brgastro::grid_cache::_ra_grid_min_ = -pi;
BRG_ANGLE brgastro::grid_cache::_ra_grid_max_ = pi;
BRG_ANGLE brgastro::grid_cache::_ra_grid_step_ = pi / 8;
BRG_ANGLE brgastro::grid_cache::_dec_grid_min_ = -pi / 2;
BRG_ANGLE brgastro::grid_cache::dec_grid_max_val = pi / 2;
BRG_ANGLE brgastro::grid_cache::_dec_grid_step_ = pi / 8;
double brgastro::grid_cache::_z_grid_min_ = 0;
double brgastro::grid_cache::_z_grid_max_ = 2;
double brgastro::grid_cache::_z_grid_step_ = 0.1;
#endif // End initialisation for brgastro::grid_cache

// Initialisation for brgastro::dfa_cache
DEFINE_BRG_CACHE_STATIC_VARS( dfa_cache, 0, 5, 0.01 );

// Initialisation for brgastro::add_cache
DEFINE_BRG_CACHE_2D_STATIC_VARS( add_cache, 0, 5, 0.01,
		                                    0, 5, 0.01);

// Initialisation for brgastro::tfa_cache
DEFINE_BRG_CACHE_STATIC_VARS( tfa_cache, 0.001, 1.02, 0.001 );

#endif // end Static Class Initialisation

/** Class Method Definitions **/
#if (1)
// brgastro::redshift_obj class methods
#if (1)
const int brgastro::redshift_obj::z_grid() const
{
	if ( _z_grid_cached_ )
	{
		if ( _local_z_grid_change_num_ == grid_cache().z_grid_change_num() )
			return _z_grid_;
	}
	_z_grid_ = brgastro::get_z_grid( z() );
	_z_grid_cached_ = true;
	_local_z_grid_change_num_ = grid_cache().z_grid_change_num();
	return _z_grid_;
}

const BRG_UNITS brgastro::redshift_obj::H() const
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
			// Friedmann equation, assuming omega = -1
			_H_cache_ = brgastro::H(_z_);
		}
		_H_cached_ = true;
	}
	return _H_cache_;
}

#endif// brgastro::dfa_cache class methods
#if (1)
const int brgastro::dfa_cache::_calculate( const double in_params, double & out_params ) const
{
	try
	{
		out_params = brgastro::integrate_add( 0, in_params );
	}
	catch( const std::exception &e)
	{
		return LOWER_LEVEL_ERROR;
	}

	return 0;
}

#endif // end brgastro::dfa_cache functions

// brgastro::add_cache class methods
const int brgastro::add_cache::_calculate( const double in_param_1, const double in_param_2,
		double  & out_param ) const
{
	out_param = brgastro::integrate_add(in_param_1,in_param_2);
	return 0;
}

const int brgastro::tfa_cache::_calculate( const double in_params, double & out_params ) const
{
	try
	{
		out_params = -brgastro::integrate_ltd( 0, brgastro::zfa( in_params ) ) / c;
	}
	catch(const std::exception &e)
	{
		std::cerr << "ERROR: Could not calculate cache for " << _name_base() << "\n"
				<< "Exception: " << e.what() << "\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

#endif // end class methods

/** Global Function Definitions **/
#if (1)

const BRG_UNITS brgastro::H( const double test_z )
{
	// Friedmann equation, assuming omega = -1
	if(test_z==0) return H_0;
	double zp1 = 1.+test_z;
	return H_0
			* std::sqrt( Omega_r * quart( zp1 )
							+ Omega_m * cube( zp1 )
							+ Omega_k * square( zp1 ) + Omega_l );
}

// grid functions
#if (1)
const int brgastro::get_ra_grid( const BRG_ANGLE &ra )
{
	grid_cache cache;
	return (int)floor(
			( ra - cache.ra_grid_min() ) / safe_d( cache.ra_grid_step() ) );
}
const int brgastro::get_dec_grid( const BRG_ANGLE &dec )
{
	grid_cache cache;
	return (int)floor( ( dec - cache.dec_grid_min() ) / cache.dec_grid_step() );
}
const int brgastro::get_z_grid( const double z )
{
	grid_cache cache;
	return (int)floor( ( z - cache.z_grid_min() ) / cache.z_grid_step() );
}

const BRG_ANGLE brgastro::get_ra_grid_lower( const int ra_grid )
{
	grid_cache cache;
	return cache.ra_grid_min() + cache.ra_grid_step() * ra_grid;
}
const BRG_ANGLE brgastro::get_dec_grid_lower( const int dec_grid )
{
	grid_cache cache;
	return cache.dec_grid_min() + cache.dec_grid_step() * dec_grid;
}
const double brgastro::get_z_grid_lower( const int z_grid )
{
	grid_cache cache;
	return cache.z_grid_min() + cache.z_grid_step() * z_grid;
}

const BRG_ANGLE brgastro::get_ra_grid_upper( const int ra_grid )
{
	grid_cache cache;
	return cache.ra_grid_min() + cache.ra_grid_step() * ( ra_grid + 1 );
}
const BRG_ANGLE brgastro::get_dec_grid_upper( const int dec_grid )
{
	grid_cache cache;
	return cache.dec_grid_min() + cache.dec_grid_step() * ( dec_grid + 1 );
}
const double brgastro::get_z_grid_upper( const int z_grid )
{
	grid_cache cache;
	return cache.z_grid_min() + cache.z_grid_step() * ( z_grid + 1 );
}

const BRG_ANGLE brgastro::get_ra_grid_mid( const int ra_grid )
{
	grid_cache cache;
	return cache.ra_grid_min() + cache.ra_grid_step() * ( ra_grid + 0.5 );
}
const BRG_ANGLE brgastro::get_dec_grid_mid( const int dec_grid )
{
	grid_cache cache;
	return cache.dec_grid_min() + cache.dec_grid_step() * ( dec_grid + 0.5 );
}
const double brgastro::get_z_grid_mid( const int z_grid )
{
	grid_cache cache;
	return cache.z_grid_min() + cache.z_grid_step() * ( z_grid + 0.5 );
}
#endif // end grid functions

// dfa and afd functions
#if (1)

const BRG_DISTANCE brgastro::dfa( const BRG_ANGLE &da, const double z )
{
	return da * dfa_cache().get( z );
}
const BRG_DISTANCE brgastro::dfa( const BRG_ANGLE &a1, const BRG_ANGLE &a2,
		const double z )
{
	return brgastro::dfa( a2 - a1, z );
}
const BRG_DISTANCE brgastro::dfa( const BRG_ANGLE &a1x, const BRG_ANGLE &a1y,
		const BRG_ANGLE &a2x, const BRG_ANGLE &a2y, const double z )
{
	return brgastro::dfa( skydist2d( a1x, a1y, a2x, a2y ), z );
}
const BRG_ANGLE brgastro::afd( const BRG_DISTANCE &dd, const double z )
{
	return dd / safe_d( dfa_cache().get( z ) );
}
const BRG_ANGLE brgastro::afd( const BRG_DISTANCE &d1, const BRG_DISTANCE &d2,
		const double z )
{
	return brgastro::afd( fabs( d2 - d1 ), z );
}
const BRG_ANGLE brgastro::afd( const BRG_DISTANCE &d1x,
		const BRG_DISTANCE &d1y, const BRG_DISTANCE &d2x,
		const BRG_DISTANCE &d2y, const double z )
{
	return brgastro::afd( dist2d( d1x, d1y, d2x, d2y ), z );
}

const double brgastro::zfa( const double a )
{
	return 1. / safe_d( a ) - 1.;
}
const double brgastro::afz( const double z )
{
	return 1. / safe_d( 1 + z );
}

const BRG_TIME brgastro::tfz( const double z )
{
	return brgastro::tfa( afz( z ) );
}
const BRG_TIME brgastro::tfa( const double a )
{
	return tfa_cache().get( a );
}
const double brgastro::zft( const BRG_TIME &t )
{
	return brgastro::zfa( brgastro::aft( t ) );
}
const double brgastro::aft( const BRG_TIME &t )
{
	brgastro::tfa_cache cache;
	return cache.inverse_get( t );
}

#endif

// Functions to integrate out distances
#if (1)
const double brgastro::integrate_add( const double z1, const double z2 )
{
	return brgastro::integrate_distance( z1, z2, 0, 100000 );
}
const double brgastro::integrate_cmd( const double z1, const double z2 )
{
	return brgastro::integrate_distance( z1, z2, 1, 10000000 );
}
const double brgastro::integrate_Ld( const double z1, const double z2 )
{
	return brgastro::integrate_distance( z1, z2, 2, 10000000 );
}
const double brgastro::integrate_ltd( const double z1, const double z2 )
{
	return brgastro::integrate_distance( z1, z2, 3, 10000000 );
}
const double brgastro::integrate_add( const double z )
{
	return brgastro::integrate_distance( 0, z, 0, 100000 );
}
const double brgastro::integrate_cmd( const double z )
{
	return brgastro::integrate_distance( 0, z, 1, 10000000 );
}
const double brgastro::integrate_Ld( const double z )
{
	return brgastro::integrate_distance( 0, z, 2, 10000000 );
}
const double brgastro::integrate_ltd( const double z )
{
	return brgastro::integrate_distance( 0, z, 3, 10000000 );
}
const double brgastro::integrate_distance( const double z1_init,
		const double z2_init, const int mode, const int n )
{
	// Function that will help calculate cosmic distances, thanks to Richard Powell - http://www.atlasoftheuniverse.com/
	// NOTE: Will return a negative value if z2 < z1. This is by design.

	double OK0;
	double OM0 = Omega_m, OR0 = Omega_r, OL0 = Omega_l;
	double OM, OR, OL, OK;
	double z1 = z1_init, z2 = z2_init;
	double HD; //Hubble distance in billions of lightyears
	double z, a, a1, a2, adot, h1;
	double DC = DBL_MAX, DCC = 0, DT = DBL_MAX, DTT = 0, DM;
	//double age, size;
	int i;
	short int sign = 1;

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

const BRG_DISTANCE brgastro::ad_distance( double z1, double z2 )
{
	if ( z2 < z1 )
		std::swap( z1, z2 );
	return brgastro::add_cache().get( z1, z2 );
}

const BRG_UNITS brgastro::sigma_crit( const double z_lens,
		const double z_source )
{
	return square( c ) / ( 4. * pi * Gc )
			* brgastro::ad_distance( 0, z_source )
			/ ( brgastro::ad_distance( 0, z_lens )
					* brgastro::ad_distance( z_lens, z_source ) );

}
#endif // end other functions

#endif // end Global function definitions
