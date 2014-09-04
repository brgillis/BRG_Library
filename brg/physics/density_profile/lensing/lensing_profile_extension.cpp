/**********************************************************************\
  @file lensing_profile_extension.cpp

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

#include <iostream>
#include <vector>

#include "brg/math/calculus/integrate.hpp"
#include "brg/math/solvers/solvers.hpp"
#include "brg/physics/units/unit_obj.h"
#include "brg/physics/density_profile/lensing/lensing_profile_extension_functors.h"
#include "brg/physics/density_profile/lensing/shifting/shifting_cache.h"

#include "lensing_profile_extension.h"

// Lensing profile extension private methods

#if(1)

const double brgastro::lensing_profile_extension::shift_factor( CONST_BRG_DISTANCE_REF R,
		const bool silent) const
{
	if(R==0) return 0;
	const BRG_ANGLE theta_separation = afd(R,z());
	return shifting_cache().get(theta_separation,z())/theta_separation;
}

// Gives the expected shift in physical coordinates for a given physical separation
const BRG_DISTANCE brgastro::lensing_profile_extension::_shift_sigma( CONST_BRG_DISTANCE_REF R,
		const bool silent) const
{
	return R*shift_factor(R,silent);
}

#endif // Lensing profile extension private methods

const BRG_UNITS brgastro::lensing_profile_extension::proj_dens( CONST_BRG_DISTANCE_REF R,
		const bool silent ) const
{
	BRG_DISTANCE R_to_use = std::fabs( R );
	double inf_factor = 20;
	brgastro::projected_density_functor func( this, R_to_use );
	BRG_UNITS min_in_params( 0 ), max_in_params( inf_factor * rvir() ),
			out_params( 0 );
	if ( R_to_use == 0 )
	{
		// In this case, we might be integrating over a singularity, so the trapezoid method is safer
		const int num_steps = 10000;

		if ( brgastro::integrate_trapezoid( &func, min_in_params,
				max_in_params, num_steps, out_params ) )
		{
			if ( !silent )
				std::cerr
						<< "WARNING: Could not integrate projected density of density profile.\n";
		}
	}
	else
	{
		out_params = brgastro::integrate_Romberg( &func, min_in_params,
				max_in_params, 0.00001, false, silent );
	}
	return 2 * out_params;
}

const BRG_MASS brgastro::lensing_profile_extension::proj_enc_mass( CONST_BRG_DISTANCE_REF R,
		const bool silent ) const
{
	if ( R == 0 )
		return 0;
	BRG_DISTANCE R_to_use = std::fabs( R );
	brgastro::cylindrical_density_functor func( this );
	BRG_UNITS min_in_params( 0 ), max_in_params( R_to_use ), out_params( 0 );
	out_params = brgastro::integrate_Romberg( &func, min_in_params,
			max_in_params, 0.00001, false, silent );
	return out_params;
}

const BRG_UNITS brgastro::lensing_profile_extension::offset_WLsig( CONST_BRG_DISTANCE_REF R,
		CONST_BRG_DISTANCE_REF offset_R, const bool silent ) const
{
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

	out_param_ring = brgastro::integrate_Romberg( &ringfunc, min_in_params_ring,
			max_in_params_ring, precision, false, silent );

	ringmean = out_param_ring / pi;

	BRG_UNITS min_in_params_circ = max(offset_R_to_use-R_to_use,0.);
	BRG_UNITS max_in_params_circ = offset_R_to_use+R_to_use;

	out_param_circ = brgastro::integrate_Romberg( &circfunc, min_in_params_circ,
			max_in_params_circ, precision, false, silent );

	circmean = out_param_circ / ( pi * square(R_to_use) );

	result = circmean - ringmean;

	return result;
}

const BRG_UNITS brgastro::lensing_profile_extension::group_WLsig( CONST_BRG_DISTANCE_REF R,
		const double group_c, const bool silent ) const
{
	BRG_DISTANCE R_to_use = std::fabs( R );
	brgastro::offset_WLsig_functor func( this, R_to_use );
	brgastro::group_WLsig_weight_functor weight_func( this, group_c );
	BRG_UNITS min_in_params( SMALL_FACTOR ), max_in_params( rvir() ),
			out_params( 0 );
	out_params = brgastro::integrate_weighted_Romberg( &func, &weight_func,
			min_in_params, max_in_params, 0.00001, false, silent);
	return out_params;
}

const BRG_UNITS brgastro::lensing_profile_extension::semiquick_group_WLsig( CONST_BRG_DISTANCE_REF R,
		const double group_c, const bool silent ) const
{
	BRG_DISTANCE R_to_use = std::fabs( R );
	brgastro::quick_offset_WLsig_functor func( this, R_to_use );
	brgastro::group_WLsig_weight_functor weight_func( this, group_c );
	BRG_UNITS min_in_params( SMALL_FACTOR ), max_in_params( 2.5*rvir() ),
			out_params( 0 );
	out_params = brgastro::integrate_weighted_Romberg( &func, &weight_func,
			min_in_params, max_in_params, 0.00001, false, silent);
	return out_params;
}

const BRG_UNITS brgastro::lensing_profile_extension::shifted_WLsig( CONST_BRG_DISTANCE_REF R,
		const bool silent ) const
{
	BRG_DISTANCE R_to_use = std::fabs( R );
	BRG_DISTANCE sigma = _shift_sigma(R_to_use);

	brgastro::shifted_WLsig_functor func( this, R_to_use );
	brgastro::shifted_WLsig_weight_functor weight_func( sigma );

	double precision = 0.00001;

	BRG_UNITS min_in_params( 0 ), max_in_params( 4*sigma ),
			out_params( 0 );

	out_params = brgastro::integrate_weighted_Romberg( &func, &weight_func,
			min_in_params, max_in_params, precision, false, silent );
	return out_params;
}
const BRG_UNITS brgastro::lensing_profile_extension::semiquick_shifted_WLsig( CONST_BRG_DISTANCE_REF R,
		const bool silent ) const
{
	BRG_DISTANCE R_to_use = std::fabs( R );
	BRG_DISTANCE sigma = _shift_sigma(R_to_use);

	brgastro::shifted_WLsig_functor func( this, R_to_use );
	brgastro::shifted_WLsig_weight_functor weight_func( sigma );

	BRG_UNITS min_in_params( 0 ), max_in_params( 4*sigma ),
			out_params( 0 );

	out_params = brgastro::integrate_weighted_Romberg( &func, &weight_func,
			min_in_params, max_in_params, 0.00001, false, silent);
	return out_params;
}


