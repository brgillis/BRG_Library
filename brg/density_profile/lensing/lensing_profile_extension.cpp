/**       @file lensing_profile_extension.cpp
 *
 *     Project: brg
 *        Path: /brg/density_profile/lensing_profile_extension.cpp
 *
 *  Created on: 12 Aug 2014
 *      Author: brg
 */

#include <iostream>
#include <vector>

#include "../../brg_calculus.hpp"
#include "../../brg_solvers.hpp"
#include "../../brg_units.h"
#include "lensing_profile_extension.h"
#include "lensing_profile_extension_functors.h"

// Lensing profile extension private methods

#if(1)

// Gives the expected shift in physical coordinates for a given angular separation
const BRG_DISTANCE brgastro::lensing_profile_extension::_shift_sigma( CONST_BRG_DISTANCE_REF R,
		const bool silent) const
{
	const BRG_ANGLE theta_separation = afd(R,z());
	const BRG_ANGLE theta_sigma = 0.01*theta_separation; // TODO: Replace placeholder here with actual function
	const BRG_DISTANCE result = dfa(theta_sigma,z());
	return result;
}

#endif // Lensing profile extension private methods

const BRG_UNITS brgastro::lensing_profile_extension::proj_dens( const BRG_DISTANCE &R,
		const bool silent ) const
{
	BRG_DISTANCE R_to_use = std::fabs( R );
	double inf_factor = 20;
	brgastro::projected_density_functor func( this, R_to_use );
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

const BRG_MASS brgastro::lensing_profile_extension::proj_enc_mass( const BRG_DISTANCE &R,
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

const BRG_UNITS brgastro::lensing_profile_extension::offset_WLsig( const BRG_DISTANCE &R,
		const BRG_DISTANCE &offset_R, const bool silent ) const
{
	unsigned int num_out_params = 1;
	if ( offset_R == 0 )
		return WLsig( R );
	BRG_DISTANCE R_to_use = std::fabs( R );
	BRG_DISTANCE offset_R_to_use = std::fabs( offset_R );
	offset_ring_dens_functor ringfunc( this, offset_R_to_use, R_to_use );
	offset_circ_dens_functor circfunc( this, offset_R_to_use, R_to_use );

	BRG_UNITS out_param_ring = 0;
	BRG_UNITS out_param_circ = 0;
	BRG_UNITS circmean;
	BRG_UNITS ringmean;
	BRG_UNITS result;

	double precision = 0.001;

	BRG_UNITS min_in_params_ring = 0;
	BRG_UNITS max_in_params_ring = pi;

	if ( brgastro::integrate_Rhomberg( &ringfunc, 1, min_in_params_ring,
			max_in_params_ring, num_out_params, out_param_ring, precision ) )
	{
		if ( !silent )
			std::cerr << "ERROR: Could not integrate in offset_WLsig" << std::endl;
		return 0;
	}

	ringmean = out_param_ring / pi;

	BRG_UNITS min_in_params_circ = max(offset_R_to_use-R_to_use,0.);
	BRG_UNITS max_in_params_circ = offset_R_to_use+R_to_use;

	if ( brgastro::integrate_Rhomberg( &circfunc, 1, min_in_params_circ,
			max_in_params_circ, num_out_params, out_param_circ, precision ) )

	{
		if ( !silent )
			std::cerr
					<< "ERROR: Could not integrate in offset_WLsig"
					<< std::endl;
		return 0;
	}

	circmean = out_param_circ / ( pi * square(R_to_use) );

	result = circmean - ringmean;

	return result;
}

const BRG_UNITS brgastro::lensing_profile_extension::group_WLsig( const BRG_DISTANCE &R,
		const double group_c, const bool silent ) const
{
	BRG_DISTANCE R_to_use = std::fabs( R );
	brgastro::offset_WLsig_functor func( this, R_to_use );
	brgastro::group_WLsig_weight_functor weight_func( this, group_c );
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

const BRG_UNITS brgastro::lensing_profile_extension::semiquick_group_WLsig( const BRG_DISTANCE &R,
		const double group_c, const bool silent ) const
{
	BRG_DISTANCE R_to_use = std::fabs( R );
	brgastro::quick_offset_WLsig_functor func( this, R_to_use );
	brgastro::group_WLsig_weight_functor weight_func( this, group_c );
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

const BRG_UNITS brgastro::lensing_profile_extension::shifted_WLsig( CONST_BRG_DISTANCE_REF R,
		const bool silent ) const
{
	BRG_DISTANCE R_to_use = std::fabs( R );
	BRG_DISTANCE sigma = _shift_sigma(R_to_use);

	brgastro::offset_WLsig_functor func( this, R_to_use );
	brgastro::shifted_WLsig_weight_functor weight_func( sigma );

	unsigned int num_in_params = 1, num_out_params = 0;

	BRG_UNITS min_in_params( 0 ), max_in_params( 4*sigma ),
			out_params( 0 );

	if ( brgastro::integrate_weighted_Rhomberg( &func, &weight_func,
			num_in_params, min_in_params, max_in_params, num_out_params,
			out_params ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Could not integrate shifted deltasigma of density profile.\n";
	}
	return out_params;
}
const BRG_UNITS brgastro::lensing_profile_extension::semiquick_shifted_WLsig( CONST_BRG_DISTANCE_REF R,
		const bool silent ) const
{
	BRG_DISTANCE R_to_use = std::fabs( R );
	BRG_DISTANCE sigma = _shift_sigma(R_to_use);

	brgastro::quick_offset_WLsig_functor func( this, R_to_use );
	brgastro::shifted_WLsig_weight_functor weight_func( sigma );

	unsigned int num_in_params = 1, num_out_params = 0;

	BRG_UNITS min_in_params( 0 ), max_in_params( 4*sigma ),
			out_params( 0 );

	if ( brgastro::integrate_weighted_Rhomberg( &func, &weight_func,
			num_in_params, min_in_params, max_in_params, num_out_params,
			out_params ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Could not integrate shifted deltasigma of density profile.\n";
	}
	return out_params;
}


