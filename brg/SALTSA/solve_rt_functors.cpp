/**       @file solve_rt_functors.cpp
 *
 *     Project: brg
 *        Path: /brg/SALTSA/solve_rt_functors.cpp
 *
 *  Created on: 29 Aug 2014
 *      Author: brg
 */

#include <cstdlib>

#include "../brg_global.h"

#include "../density_profile/density_profile.h"
#include "../brg_misc_functions.hpp"
#include "../brg_units.h"

#include "solve_rt_functors.h"

// brgastro::solve_rt_grid_functor class method implementations
#if (1)

BRG_UNITS brgastro::solve_rt_grid_functor::operator()(
		CONST_BRG_UNITS_REF  in_param, const bool silent ) const
{
	BRG_DISTANCE r;
	BRG_MASS delta_M;
	r = std::fabs( in_param );

	delta_M = sum_delta_rho * 4. / 3. * pi * cube( r );

	if ( r == 0 )
	{
		return DBL_MAX;
	}
	else
	{
		return std::fabs( Gc * ( satellite_ptr->enc_mass( r ) - delta_M )
						/ safe_d(( omega * omega + Daccel ) * cube( r ) ) - 1 );
	}
	return 0;
}

brgastro::solve_rt_grid_functor::solve_rt_grid_functor()
{
	omega = 0;
	Daccel = 0;
	sum_delta_rho = 0;
	satellite_ptr = NULL;
	return;
}
brgastro::solve_rt_grid_functor::solve_rt_grid_functor(
		const BRG_UNITS init_omega, const density_profile *init_satellite,
		const BRG_UNITS init_Daccel, const long double init_sum_delta_rho )
{
	omega = init_omega;
	satellite_ptr = init_satellite;
	Daccel = init_Daccel;
	sum_delta_rho = init_sum_delta_rho;
	return;
}

#endif // end brgastro::solve_rt_grid_functor class function definitions

// brgastro::solve_rt_grid_functor class method implementations
#if(1)

BRG_UNITS brgastro::solve_rt_it_functor::operator()(
		CONST_BRG_UNITS_REF in_param, const bool silent ) const
{
	BRG_DISTANCE r = 0;
	BRG_MASS delta_M = 0;
	BRG_UNITS r3 = 0;

	r = std::fabs( in_param );

	delta_M = sum_delta_rho * 4. / 3. * pi * cube( r );

	r3 = Gc * ( satellite_ptr->enc_mass( r ) - delta_M )
			/ safe_d( omega * omega + Daccel );
	if ( r3 <= 0 )
	{
		return r * 0.9;
	}
	else
	{
		return std::pow( r3, 1. / 3. );
	}
	return 0;
}

brgastro::solve_rt_it_functor::solve_rt_it_functor()
{
	omega = 0;
	satellite_ptr = NULL;
	Daccel = 0;
	sum_delta_rho = 0;
	return;
}
brgastro::solve_rt_it_functor::solve_rt_it_functor(
		const BRG_UNITS init_omega, const density_profile *init_satellite,
		const BRG_UNITS init_Daccel, const long double init_sum_delta_rho )
{
	omega = init_omega;
	satellite_ptr = init_satellite;
	Daccel = init_Daccel;
	sum_delta_rho = init_sum_delta_rho;
	return;
}

#endif // end brgastro::solve_rt_it_functor class function definitions
