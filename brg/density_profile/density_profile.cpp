/**       @file density_profile.cpp
 *
 *     Project: brg
 *        Path: /brg/density_profile/density_profile.cpp
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#include <iostream>

#include "../brg_global.h"

#include "../brg_calculus.hpp"
#include "../brg_solvers.hpp"
#include "../brg_units.h"
#include "density_profile_functors.h"
#include "density_profile.h"

// brgastro::density_profile class methods
#if (1)
const BRG_UNITS brgastro::density_profile::Daccel( CONST_BRG_DISTANCE_REF r,
		const bool silent ) const
{
	BRG_DISTANCE dr;
	BRG_UNITS a1, a2;
	// It's simplest, and a ton faster to just manually differentiate here.
	dr = max( r * SMALL_FACTOR, SMALL_FACTOR );
	a1 = accel( r, silent );
	a2 = accel( r + dr, silent );
	return ( a2 - a1 ) / safe_d( dr );
}

const BRG_DISTANCE brgastro::density_profile::rhmtot( const bool silent ) const
{
	// If cached, return the cached value
	if ( hmtot_cached )
		return _rhmtot_cache_;

	// Not cached, so we'll have to calculate it
	double target_mass = hmtot();
	solve_rhm_functor func( this, target_mass );
	solve_rhm_functor *func_ptr = &func;

	BRG_UNITS max_r( default_tau_factor * rvir() );
	BRG_UNITS rhm_test( 0 );

	// First check for zero mass/radius/density
	if ( ( mvir() <= 0 ) || ( rvir() <= 0 ) || ( dens( rvir() / 2 ) < 0 ) )
	{
		hmtot_cached = true;
		return _rhmtot_cache_ = 0;
	}

	if ( solve_grid( func_ptr, 1, 0., max_r, 10, 0., rhm_test ) ) // If we can't solve it
	{
		if ( !silent )
			std::cerr << "WARNING: Could not solve half-mass radius. Assuming it's zero.\n";

		_rhmtot_cache_ = 0;
		hmtot_cached = true;
		return _rhmtot_cache_;
	}
	else
	{
		_rhmtot_cache_ = std::fabs( rhm_test );
		hmtot_cached = true;
	}
	return _rhmtot_cache_;
}

const BRG_DISTANCE brgastro::density_profile::rhmvir( const bool silent ) const
{
	// If cached, return the cached value
	if ( hmvir_cached )
		return _rhmvir_cache_;

	// Not cached, so we'll have to calculate it
	double target_mass = hmvir();
	solve_rhm_functor func( this, target_mass );
	solve_rhm_functor *func_ptr = &func;

	BRG_UNITS max_r( default_tau_factor * rvir() );
	BRG_UNITS rhm_test( 0 );

	// First check for zero mass/radius/density
	if ( ( mvir() <= 0 ) || ( rvir() <= 0 ) || ( dens( rvir() / 2 ) < 0 ) )
	{
		hmvir_cached = true;
		return _rhmvir_cache_ = 0;
	}

	if ( solve_grid( func_ptr, 1, 0., max_r, 10, 0., rhm_test ) ) // If we can't solve it
	{
		if ( !silent )
			std::cerr << "WARNING: Could not solve half-mass radius.\n";
		return -1;
	}
	else
	{
		_rhmvir_cache_ = max(0,std::fabs( rhm_test ));
		hmvir_cached = true;
	}
	return _rhmvir_cache_;
}

const BRG_MASS brgastro::density_profile::enc_mass( CONST_BRG_DISTANCE_REF r,
		const bool silent ) const
{
	if ( r == 0 )
		return 0;
	BRG_DISTANCE r_to_use = std::fabs( r );
	brgastro::spherical_density_functor func( this );
	unsigned int num_in_params = 1, num_out_params = 0;
	BRG_UNITS min_in_params( 0 ), max_in_params( r_to_use ), out_params( 0 );
	if ( brgastro::integrate_Romberg( &func, num_in_params, min_in_params,
			max_in_params, num_out_params, out_params ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Could not integrate enclosed mass of density profile.\n";
	}
	return out_params;
}

#endif

const BRG_TIME brgastro::period( const brgastro::density_profile *host,
		CONST_BRG_DISTANCE_REF r, CONST_BRG_VELOCITY_REF vr, CONST_BRG_VELOCITY_REF vt )
{
	BRG_UNITS mu = host->enc_mass( r ) * Gc;
	BRG_VELOCITY v = quad_add( vr, vt );
	BRG_DISTANCE a = -mu / 2 / safe_d( v * v / 2 - mu / safe_d( r ) );
	BRG_TIME result = (
			a > 0 ? 2 * pi * std::sqrt( cube(a) / mu ) : BRG_TIME( 0 ) );
	return result;
}
