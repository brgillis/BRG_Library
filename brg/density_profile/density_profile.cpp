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
const BRG_UNITS brgastro::density_profile::Daccel( const BRG_DISTANCE &r,
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

const BRG_MASS brgastro::density_profile::enc_mass( const BRG_DISTANCE &r,
		const bool silent ) const
{
	if ( r == 0 )
		return 0;
	BRG_DISTANCE r_to_use = std::fabs( r );
	brgastro::spherical_density_functor func( this );
	unsigned int num_in_params = 1, num_out_params = 0;
	BRG_UNITS min_in_params( 0 ), max_in_params( r_to_use ), out_params( 0 );
	if ( brgastro::integrate_Rhomberg( &func, num_in_params, min_in_params,
			max_in_params, num_out_params, out_params ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Could not integrate enclosed mass of density profile.\n";
	}
	return out_params;
}

const BRG_UNITS brgastro::density_profile::proj_dens( const BRG_DISTANCE &R,
		const bool silent ) const
{
	BRG_DISTANCE R_to_use = std::fabs( R );
	double inf_factor = 20;
	brgastro::cylindrical_density_functor func( this );
	unsigned int num_in_params = 1, num_out_params = 0;
	BRG_UNITS min_in_params( 0 ), max_in_params( inf_factor * rvir() ),
			out_params( 0 );
	if ( R_to_use == 0 )
	{
		// In this case, we might be integrating over a singularity, so the trapezoid method is safer
		const int num_steps = 10000;

		if ( brgastro::integrate( &func, num_in_params, min_in_params,
				max_in_params, num_steps, num_out_params, out_params ) )
		{
			if ( !silent )
				std::cerr
						<< "WARNING: Could not integrate projected density of density profile.\n";
		}
	}
	else
	{
		if ( brgastro::integrate_Rhomberg( &func, num_in_params, min_in_params,
				max_in_params, num_out_params, out_params ) )
		{
			if ( !silent )
				std::cerr
						<< "WARNING: Could not integrate projected density of density profile.\n";
		}
	}
	return 2 * out_params;
}

const BRG_MASS brgastro::density_profile::proj_enc_mass( const BRG_DISTANCE &R,
		const bool silent ) const
{
	if ( R == 0 )
		return 0;
	BRG_DISTANCE R_to_use = std::fabs( R );
	brgastro::cylindrical_density_functor func( this );
	unsigned int num_in_params = 1, num_out_params = 0;
	BRG_UNITS min_in_params( 0 ), max_in_params( R_to_use ), out_params( 0 );
	if ( brgastro::integrate_Rhomberg( &func, num_in_params, min_in_params,
			max_in_params, num_out_params, out_params ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Could not integrate projected enclosed mass of density profile.\n";
	}
	return out_params;
}

const BRG_UNITS brgastro::density_profile::offset_WLsig( const BRG_DISTANCE &R,
		const BRG_DISTANCE &offset_R, const bool silent ) const
{
	unsigned int num_out_params = 1;
	if ( offset_R == 0 )
		return deltasigma( R );
	BRG_DISTANCE R_to_use = std::fabs( R );
	BRG_DISTANCE offset_R_to_use = std::fabs( R );
	offset_ring_dens_functor ringfunc( this, offset_R_to_use, R_to_use );
	offset_circ_dens_functor circfunc( this, offset_R_to_use );

	BRG_UNITS min_in_params_ring( 0 );
	BRG_UNITS max_in_params_ring( 0 );
	BRG_UNITS out_params_ring( 0 );
	std::vector< BRG_UNITS > min_in_params_circ( 2, 0 );
	std::vector< BRG_UNITS > max_in_params_circ( 2, 0 );
	std::vector< BRG_UNITS > out_params_circ( 1, 0 );
	BRG_UNITS circmean;
	BRG_UNITS ringmean;
	BRG_UNITS result;

	double precision = 0.001;

	min_in_params_ring = 0;
	max_in_params_ring = pi; // We'll double it later, due to symmetry

	min_in_params_circ[0] = R_to_use * precision;
	max_in_params_circ[0] = R_to_use;
	min_in_params_circ[1] = 0;
	max_in_params_circ[1] = pi; // We'll double it later, due to symmetry

	if ( brgastro::integrate_Rhomberg( &ringfunc, 1, min_in_params_ring,
			max_in_params_ring, num_out_params, out_params_ring, precision ) )
	{
		if ( !silent )
			std::cerr << "ERROR: Could not integrate in offset_WLsig" << std::endl;
		return 0;
	}

	ringmean = 2 * out_params_ring / pi;

	if ( brgastro::integrate_Rhomberg( &circfunc, 2, min_in_params_circ,
			max_in_params_circ, num_out_params, out_params_circ, precision ) )

	{
		if ( !silent )
			std::cerr
					<< "ERROR: Could not integrate in unitless_NFW_offset_WLsig"
					<< std::endl;
		return 0;
	}

	circmean = out_params_circ[0] / ( pi * R_to_use );

	result = circmean - ringmean;

	return result;
}
const BRG_UNITS brgastro::density_profile::group_WLsig( const BRG_DISTANCE &r,
		const double group_c, const bool silent ) const
{
	BRG_DISTANCE r_to_use = std::fabs( r );
	brgastro::offset_WLsig_functor func( this, r_to_use );
	brgastro::offset_WLsig_weight_functor weight_func( this, group_c );
	unsigned int num_in_params = 1, num_out_params = 0;
	BRG_UNITS min_in_params( SMALL_FACTOR ), max_in_params( rvir() ),
			out_params( 0 );
	if ( brgastro::integrate_weighted_Rhomberg( &func, &weight_func,
			num_in_params, min_in_params, max_in_params, num_out_params,
			out_params ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Could not integrate group deltasigma of density profile.\n";
	}
	return out_params;
}

#endif

const BRG_TIME brgastro::period( const brgastro::density_profile *host,
		const BRG_DISTANCE &r, const BRG_VELOCITY &vr, const BRG_VELOCITY &vt )
{
	BRG_UNITS mu = host->enc_mass( r ) * Gc;
	BRG_VELOCITY v = quad_add( vr, vt );
	BRG_DISTANCE a = -mu / 2 / safe_d( v * v / 2 - mu / safe_d( r ) );
	BRG_TIME result = (
			a > 0 ? 2 * pi * std::sqrt( cube(a) / mu ) : BRG_TIME( 0 ) );
	return result;
}
