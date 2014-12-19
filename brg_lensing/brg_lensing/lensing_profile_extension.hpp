/**********************************************************************\
  @file lensing_profile_extension.hpp

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

#ifndef _BRG_LENSING_PROFILE_EXTENSION_HPP_
#define _BRG_LENSING_PROFILE_EXTENSION_HPP_

#include "brg/global.h"

#include "brg_lensing/lensing_profile_extension_functors.hpp"
#include "brg_lensing/shifting/shifting_cache.h"

#include "brg_physics/density_profile/density_profile.h"
#include "brg_physics/units/unit_obj.h"

// Macro definitions

// SPCP: "Static Polymorphic Const Pointer"
#define SPCP(name) static_cast<const name*>(this)

// SPP: "Static Polymorphic Pointer"
#define SPP(name) static_cast<name*>(this)

namespace brgastro {

template<typename name>
class lensing_profile_extension {
private:

// Calculation functions
#if(1)

// Basic calculation functions which should be overridden if possible
#if(1) // Basic calculation functions which should be overridden if possible

	const BRG_UNITS _proj_dens( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const // Projected surface density at radius R
	{
		BRG_DISTANCE R_to_use = std::fabs( R );
		double inf_factor = 20;
		brgastro::projected_density_functor<name> func( SPCP(name), R_to_use );
		BRG_UNITS min_in_params( 0 ), max_in_params( inf_factor * SPCP(name)->rvir() ),
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

	const BRG_MASS _proj_enc_mass( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const // Mass enclosed within a cylinder of radius R
	{
		if ( R == 0 )
			return 0;
		BRG_DISTANCE R_to_use = std::fabs( R );
		brgastro::cylindrical_density_functor<name> func( SPCP(name) );
		BRG_UNITS min_in_params( 0 ), max_in_params( R_to_use ), out_params( 0 );
		out_params = brgastro::integrate_Romberg( &func, min_in_params,
				max_in_params, 0.00001, false, silent );
		return out_params;
	}

#endif

// Advanced calculations which usually can't be overridden, but should be if possible
#if(1)

	const BRG_UNITS _offset_WLsig( CONST_BRG_DISTANCE_REF R,
			CONST_BRG_DISTANCE_REF offset_R, const bool silent = true ) const // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from position offset by offset_R
	{
		if ( offset_R == 0 )
			return SPCP(name)->WLsig( R );
		BRG_DISTANCE R_to_use = std::fabs( R );
		BRG_DISTANCE offset_R_to_use = std::fabs( offset_R );
		offset_ring_dens_functor<name> ringfunc( SPCP(name), offset_R_to_use, R_to_use );
		offset_circ_dens_functor<name> circfunc( SPCP(name), offset_R_to_use, R_to_use );

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
	const BRG_UNITS _group_WLsig( CONST_BRG_DISTANCE_REF R,
			const double group_c, const bool silent = true ) const // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from ensemble of satellites in group with satellite concentration group_c
	{
		BRG_DISTANCE R_to_use = std::fabs( R );
		brgastro::offset_WLsig_functor<name> func( SPCP(name), R_to_use );
		brgastro::group_WLsig_weight_functor<name> weight_func( SPCP(name), group_c );
		BRG_UNITS min_in_params( SMALL_FACTOR ), max_in_params( 2.5*SPCP(name)->rvir() ),
				out_params( 0 );
		out_params = brgastro::integrate_weighted_Romberg( &func, &weight_func,
				min_in_params, max_in_params, 0.00001, false, silent);
		return out_params;
	}
	const BRG_UNITS _semiquick_group_WLsig( CONST_BRG_DISTANCE_REF R,
			const double group_c, const bool silent = true ) const // As group_WLsig, but uses offset_WLsig cache to speed it up if overwritten
	{
		BRG_DISTANCE R_to_use = std::fabs( R );
		brgastro::quick_offset_WLsig_functor<name> func( SPCP(name), R_to_use );
		brgastro::group_WLsig_weight_functor<name> weight_func( SPCP(name), group_c );
		BRG_UNITS min_in_params( SMALL_FACTOR ), max_in_params( 2.5*SPCP(name)->rvir() ),
				out_params( 0 );
		out_params = brgastro::integrate_weighted_Romberg( &func, &weight_func,
				min_in_params, max_in_params, 0.00001, false, silent);
		return out_params;
	}

	const BRG_UNITS _shifted_WLsig( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const
	{
		BRG_DISTANCE R_to_use = std::fabs( R );
		BRG_DISTANCE sigma = SPCP(name)->_shift_sigma(R_to_use);

		brgastro::shifted_WLsig_functor<name> func( SPCP(name), R_to_use );
		brgastro::shifted_WLsig_weight_functor<name> weight_func( sigma );

		double precision = 0.00001;

		BRG_UNITS min_in_params( 0 ), max_in_params( 4*sigma ),
				out_params( 0 );

		out_params = brgastro::integrate_weighted_Romberg( &func, &weight_func,
				min_in_params, max_in_params, precision, false, silent );
		return out_params;
	}
	const BRG_UNITS _semiquick_shifted_WLsig( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const
	{
		BRG_DISTANCE R_to_use = std::fabs( R );
		BRG_DISTANCE sigma = SPCP(name)->_shift_sigma(R_to_use);

		brgastro::shifted_WLsig_functor<name> func( SPCP(name), R_to_use );
		brgastro::shifted_WLsig_weight_functor<name> weight_func( sigma );

		BRG_UNITS min_in_params( 0 ), max_in_params( 4*sigma ),
				out_params( 0 );

		out_params = brgastro::integrate_weighted_Romberg( &func, &weight_func,
				min_in_params, max_in_params, 0.00001, false, silent);
		return out_params;
	}

#endif // Advanced calculations which usually can't be overridden, but should be if possible

// Quick functions - should be overridden if a cache is implemented for the halo
#if(1)

	const BRG_UNITS _quick_offset_WLsig( CONST_BRG_DISTANCE_REF R,
			CONST_BRG_DISTANCE_REF offset_R, const bool silent = true ) const // As offset_WLsig, but uses cache to speed it up if overwritten
	{
		return SPCP(name)->offset_WLsig( R, offset_R, silent );
	}
	const BRG_UNITS _quick_group_WLsig( CONST_BRG_DISTANCE_REF R,
			const double group_c, const bool silent = true ) const // As deltasigma, but uses group_WLsig cache to speed it up if overwritten
	{
		return SPCP(name)->semiquick_group_WLsig( R, group_c, silent );
	}
	const BRG_UNITS _quick_shifted_WLsig( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const
	{
		return SPCP(name)->semiquick_shifted_WLsig( R, silent );
	}

#endif // Quick functions - should be overridden if a cache is implemented for the halo

// Simple calculation functions that should rarely need to be overridden
#if (1)
	const BRG_UNITS _proj_enc_dens( CONST_BRG_DISTANCE_REF R,
			const bool silent = false ) const // Mean surface density enclosed within a cylinder of radius R
	{
		BRG_DISTANCE R_to_use = max( std::fabs( R ), SMALL_FACTOR );
		return SPCP(name)->proj_enc_mass( R_to_use, silent )
				/ ( pi * square( R_to_use ) );
	}

	const BRG_DISTANCE _shift_sigma( CONST_BRG_DISTANCE_REF R, const bool silent = true ) const
	{
		return R*shift_factor(R,silent);
	}
	const double _shift_factor( CONST_BRG_DISTANCE_REF R, const bool silent = true ) const
	{
		if(R==0) return 0;
		const BRG_ANGLE theta_separation = afd(R,SPCP(name)->z());
		return shifting_cache().get(theta_separation,SPCP(name)->z())/theta_separation;
	}
	virtual const BRG_UNITS _WLsig( CONST_BRG_DISTANCE_REF R, const bool silent =
			false ) const // Weak lensing signal in tangential shear Delta-Sigma at radius R
	{
		return SPCP(name)->proj_enc_dens( R, silent ) - SPCP(name)->proj_dens( R, silent );
	}
#endif // Simple calculation functions that should rarely need to be overridden

#endif // Calculation functions

public:

	// Constructors and destructors
#if (1)
	lensing_profile_extension()
	{
	}
	virtual ~lensing_profile_extension()
	{
	}
#endif

	// Virtual functions that will need to be implemented by base classes
#if (1)

	// Clone function for this
	name* lensing_profile_extension_clone() const
	{
		return new name(*SPCP(name));
	}

#endif // Virtual functions that will need to be implemented by base classes

	// Projected mass and density functions
#if (1)

	// These ones should be overridden if at all possible, as otherwise they'll have to be
	// integrated
#if (1)
	const BRG_UNITS proj_dens( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const // Projected surface density at radius R
	{
		return SPCP(name)->_proj_dens(R,silent);
	}
	const BRG_MASS proj_enc_mass( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const // Mass enclosed within a cylinder of radius R
	{
		return SPCP(name)->_proj_enc_mass(R,silent);
	}
#endif

	// These ones typically don't need to be overridden, but they can be if it's convenient
#if (1)
	const BRG_UNITS proj_enc_dens( CONST_BRG_DISTANCE_REF R,
			const bool silent = false ) const // Mean surface density enclosed within a cylinder of radius R
	{
		return SPCP(name)->_proj_enc_dens(R,silent);
	}
#endif

#endif // Projected mass and density functions

	// Weak lensing functions
#if (1)

	// WLsig (delta sigma)
#if (1)
	const BRG_UNITS WLsig( CONST_BRG_DISTANCE_REF R, const bool silent =
			false ) const // Weak lensing signal in tangential shear Delta-Sigma at radius R
	{
		return SPCP(name)->_WLsig(R,silent);
	}
	const BRG_UNITS quick_WLsig( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const // As deltasigma, but uses cache to speed it up if overwritten
	{
		return SPCP(name)->_quick_WLsig( R, silent );
	}
#endif

	// Offset WLsig
#if (1)
	const BRG_UNITS offset_WLsig( CONST_BRG_DISTANCE_REF R,
			CONST_BRG_DISTANCE_REF offset_R, const bool silent = true ) const // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from position offset by offset_R
	{
		return SPCP(name)->_offset_WLsig( R, offset_R, silent );
	}
	const BRG_UNITS quick_offset_WLsig( CONST_BRG_DISTANCE_REF R,
			CONST_BRG_DISTANCE_REF offset_R, const bool silent = true ) const // As offset_WLsig, but uses cache to speed it up if overwritten
	{
		return SPCP(name)->_quick_offset_WLsig( R, offset_R, silent );
	}
#endif

	// Group WLsig
#if (1)
	const BRG_UNITS group_WLsig( CONST_BRG_DISTANCE_REF R,
			const double group_c=2.5, const bool silent = true ) const // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from ensemble of satellites in group with satellite concentration group_c
	{
		return SPCP(name)->_group_WLsig( R, group_c, silent );
	}
	const BRG_UNITS semiquick_group_WLsig( CONST_BRG_DISTANCE_REF R,
			const double group_c=2.5, const bool silent = true ) const // As group_WLsig, but uses offset_WLsig cache to speed it up if overwritten
	{
		return SPCP(name)->_semiquick_group_WLsig( R, group_c, silent );
	}
	const BRG_UNITS quick_group_WLsig( CONST_BRG_DISTANCE_REF R,
			const double group_c=2.5, const bool silent = true ) const // As deltasigma, but uses group_WLsig cache to speed it up if overwritten
	{
		return SPCP(name)->_quick_group_WLsig( R, group_c, silent );
	}
#endif

	// Shifted WLsig
#if (1)
	// This represents the weak-lensing signal after being corrected for errors due to relative
	// shifting of the lens and source due to an intervening mass distribution

	const double shift_factor( CONST_BRG_DISTANCE_REF R, const bool silent = true ) const
	{
		return SPCP(name)->_shift_factor( R, silent );
	}

	const BRG_UNITS shifted_WLsig( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const
	{
		return SPCP(name)->_shifted_WLsig( R, silent );
	}
	const BRG_UNITS semiquick_shifted_WLsig( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const
	{
		return SPCP(name)->_semiquick_shifted_WLsig( R, silent );
	}
	const BRG_UNITS quick_shifted_WLsig( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const
	{
		return SPCP(name)->_quick_shifted_WLsig( R, silent );
	}
#endif

#endif

};

} // end namespace brgastro

// Undef macros
#undef SPP
#undef SPCP

#endif /* _BRG_LENSING_PROFILE_EXTENSION_H_ */
