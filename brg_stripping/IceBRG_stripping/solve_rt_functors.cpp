/**********************************************************************\
  @file solve_rt_functors.cpp

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

#include "IceBRG_main/common.hpp"

#include "IceBRG_physics/density_profile/detail/density_profile.hpp"
#include "IceBRG_main/units/units.hpp"
#include "IceBRG_main/utility.hpp"

#include "solve_rt_functors.h"

using namespace IceBRG;

// IceBRG::solve_rt_grid_functor class method implementations
#if (1)

flt_t IceBRG::solve_rt_grid_functor::operator()(
		flt_t const &  in_param, const bool silent ) const
{
	distance_type r;
	mass_type delta_M;
	r = std::fabs( in_param );

	delta_M = sum_delta_rho * 4. / 3. * pi * cube( r );

	if ( r == 0 )
	{
		return std::numeric_limits<double>::max();
	}
	else
	{
		return std::fabs( Gc * ( satellite_ptr->enc_mass( r ) - delta_M )
						/ safe_d(( omega * omega + Daccel ) * cube( r ) ) - 1 );
	}
	return 0;
}

IceBRG::solve_rt_grid_functor::solve_rt_grid_functor()
{
	omega = 0;
	Daccel = 0;
	sum_delta_rho = 0;
	satellite_ptr = NULL;
	return;
}
IceBRG::solve_rt_grid_functor::solve_rt_grid_functor(
		const flt_t init_omega, const density_profile *init_satellite,
		const flt_t init_Daccel, const long double init_sum_delta_rho )
{
	omega = init_omega;
	satellite_ptr = init_satellite;
	Daccel = init_Daccel;
	sum_delta_rho = init_sum_delta_rho;
	return;
}

#endif // end IceBRG::solve_rt_grid_functor class function definitions

// IceBRG::solve_rt_grid_functor class method implementations
#if(1)

flt_t IceBRG::solve_rt_it_functor::operator()(
		flt_t const & in_param, const bool silent ) const
{
	distance_type r = 0;
	mass_type delta_M = 0;
	flt_t r3 = 0;

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

IceBRG::solve_rt_it_functor::solve_rt_it_functor()
{
	omega = 0;
	satellite_ptr = NULL;
	Daccel = 0;
	sum_delta_rho = 0;
	return;
}
IceBRG::solve_rt_it_functor::solve_rt_it_functor(
		const flt_t init_omega, const density_profile *init_satellite,
		const flt_t init_Daccel, const long double init_sum_delta_rho )
{
	omega = init_omega;
	satellite_ptr = init_satellite;
	Daccel = init_Daccel;
	sum_delta_rho = init_sum_delta_rho;
	return;
}

#endif // end IceBRG::solve_rt_it_functor class function definitions
