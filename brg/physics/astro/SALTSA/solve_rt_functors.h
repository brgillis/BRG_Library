/**       @file solve_rt_functors.h
 *
 *     Project: brg
 *        Path: /brg/SALTSA/solve_rt_functors.h
 *
 *  Created on: 29 Aug 2014
 *      Author: brg
 */

#ifndef _SOLVE_RT_FUNCTORS_H_INCLUDED_
#define _SOLVE_RT_FUNCTORS_H_INCLUDED_

#include "brg/global.h"

#include "brg/physics/astro/density_profile/density_profile.h"
#include "brg/physics/units/unit_obj.h"

namespace brgastro {

class solve_rt_it_functor // Always uses one param, returns new value for that parameter for iteration.
{
	/************************************************************
	 solve_rt_it_functor
	 --------------------

	 Provides a functor to be used to iteratively solve
	 for tidal radius.

	 \************************************************************/
public:
	const density_profile *satellite_ptr;

	long double sum_delta_rho, Daccel, omega;
	BRG_UNITS operator()( CONST_BRG_UNITS_REF in_param, const bool silent = false ) const;
	solve_rt_it_functor( const BRG_UNITS init_omega,
			const density_profile *init_satellite, const BRG_UNITS init_Daccel,
			const long double init_sum_delta_rho = 0 );

	solve_rt_it_functor();
}; // class solve_rt_it_functor

class solve_rt_grid_functor
{
	/************************************************************
	 solve_rt_grid_functor
	 ----------------------

	 Provides a functor to be used with the grid solver.

	 \************************************************************/
public:
	const density_profile *satellite_ptr;

	long double sum_delta_rho, Daccel, omega;
	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param, const bool silent = false ) const;
	solve_rt_grid_functor( const BRG_UNITS init_omega,
			const density_profile *init_satellite, const BRG_UNITS init_Daccel,
			const long double init_sum_delta_rho = 0 );

	solve_rt_grid_functor();
}; // class solve_rt_grid_functor

} // end namespace brgastro

#endif /* SOLVE_RT_FUNCTORS_H_INCLUDED_ */
