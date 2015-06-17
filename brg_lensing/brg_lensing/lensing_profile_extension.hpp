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

#include "brg/common.h"

#include "brg/math/interpolator/interpolator.h"

#include "brg_lensing/lensing_profile_extension_functors.hpp"
#include "brg_lensing/shifting/shifting_cache.h"
#include "brg_lensing/two_halo_term_functions.hpp"

#include "brg_physics/density_profile/density_profile.h"
#include "brg_physics/units/unit_obj.h"

// TODO: Update lensing tNFW_profile along with updates here

// Macro definitions

// SPCP: "Static Polymorphic Const Pointer"
#define SPCP(name) static_cast<const name*>(this)

namespace brgastro {

constexpr flt_type delta_c = 1.686;

template<typename name>
class lensing_profile_extension {
private:

// Calculation functions
#if(1)

// Basic calculation functions which should be overridden if possible
#if(1)

	const flt_type _proj_dens( const flt_type & R,
			const bool silent = false ) const // Projected surface density at radius R
	{
		flt_type R_to_use = std::fabs( R );
		flt_type inf_factor = 20;
		brgastro::projected_density_functor<name> func( SPCP(name), R_to_use );
		flt_type min_in_params( 0 ), max_in_params( inf_factor * SPCP(name)->rvir() ),
				out_params( 0 );
		if ( R_to_use == 0 )
		{
			// In this case, we might be integrating over a singularity, so the trapezoid method is safer
			const int_type num_steps = 10000;

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

	const flt_type _proj_enc_mass( const flt_type & R,
			const bool silent = true ) const // Mass enclosed within a cylinder of radius R
	{
		if ( R == 0 )
			return 0;
		flt_type R_to_use = std::fabs( R );
		brgastro::cylindrical_density_functor<name> func( SPCP(name) );
		flt_type min_in_params( 0 ), max_in_params( R_to_use ), out_params( 0 );
		out_params = brgastro::integrate_Romberg( &func, min_in_params,
				max_in_params, 0.00001, false, silent );
		return out_params;
	}

#endif

// Advanced calculations which usually can't be overridden, but should be if possible
#if(1)

	flt_type _offset_Delta_Sigma( const flt_type & R,
			const flt_type & offset_R, const bool silent = false ) const // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from position offset by offset_R
	{
		if ( offset_R == 0 )
			return SPCP(name)->Delta_Sigma( R );
		flt_type R_to_use = std::fabs( R );
		flt_type offset_R_to_use = std::fabs( offset_R );
		offset_ring_dens_functor<name> ringfunc( SPCP(name), offset_R_to_use, R_to_use );
		offset_circ_dens_functor<name> circfunc( SPCP(name), offset_R_to_use, R_to_use );

		flt_type out_param_ring = 0;
		flt_type out_param_circ = 0;
		flt_type circmean;
		flt_type ringmean;
		flt_type result;

		flt_type precision = 0.001;

		flt_type min_in_params_ring = 0;
		flt_type max_in_params_ring = pi;

		out_param_ring = brgastro::integrate_Romberg( &ringfunc, min_in_params_ring,
				max_in_params_ring, precision, false );

		ringmean = out_param_ring / pi;

		flt_type min_in_params_circ = max(offset_R_to_use-R_to_use,0.);
		flt_type max_in_params_circ = offset_R_to_use+R_to_use;

		out_param_circ = brgastro::integrate_Romberg( &circfunc, min_in_params_circ,
				max_in_params_circ, precision, false );

		circmean = out_param_circ / ( pi * square(R_to_use) );

		result = circmean - ringmean;

		return result;
	}
	flt_type _group_Delta_Sigma( const flt_type & R,
			const flt_type & group_c, const bool silent = false ) const // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from ensemble of satellites in group with satellite concentration group_c
	{
		flt_type R_to_use = std::fabs( R );
		brgastro::offset_Delta_Sigma_functor<name> func( SPCP(name), R_to_use );
		brgastro::group_Delta_Sigma_weight_functor<name> weight_func( SPCP(name), group_c );
		flt_type min_in_params( SMALL_FACTOR ), max_in_params( 2.5*SPCP(name)->rvir() ),
				out_params( 0 );
		out_params = brgastro::integrate_weighted_Romberg( &func, &weight_func,
				min_in_params, max_in_params, 0.00001, false, silent);
		return out_params;
	}
	flt_type _semiquick_group_Delta_Sigma( const flt_type & R,
			const flt_type & group_c, const bool silent = false ) const // As _group_Delta_Sigma, but uses offset_Delta_sigma cache to speed it up if overwritten
	{
		flt_type R_to_use = std::fabs( R );
		brgastro::quick_offset_Delta_Sigma_functor<name> func( SPCP(name), R_to_use );
		brgastro::group_Delta_Sigma_weight_functor<name> weight_func( SPCP(name), group_c );
		flt_type min_in_params( SMALL_FACTOR ), max_in_params( 2.5*SPCP(name)->rvir() ),
				out_params( 0 );
		out_params = brgastro::integrate_weighted_Romberg( &func, &weight_func,
				min_in_params, max_in_params, 0.00001, false);
		return out_params;
	}
	flt_type _two_halo_Delta_Sigma( const flt_type & R,
			const flt_type & group_c, const bool silent = false ) const
	{
		assert(false); // TODO: Fill this out
		return flt_type();
	}

	flt_type _offset_Sigma( const flt_type & R,
			const flt_type & offset_R, const bool silent = false ) const // Expected Sigma at radius R from position offset by offset_R
	{
		if ( offset_R == 0 )
			return SPCP(name)->proj_dens( R );
		flt_type R_to_use = std::fabs( R );
		flt_type offset_R_to_use = std::fabs( offset_R );
		offset_ring_dens_functor<name> ringfunc( SPCP(name), offset_R_to_use, R_to_use );

		flt_type out_param_ring = 0;
		flt_type ringmean;

		flt_type precision = 0.001;

		flt_type min_in_params_ring = 0;
		flt_type max_in_params_ring = pi;

		out_param_ring = brgastro::integrate_Romberg( &ringfunc, min_in_params_ring,
				max_in_params_ring, precision, false );

		ringmean = out_param_ring / pi;

		return ringmean;
	}
	flt_type _group_Sigma( const flt_type & R,
			const flt_type & group_c, const bool silent = false ) const // Expected Sigma at radius R from ensemble of satellites in group with satellite concentration group_c
	{
		flt_type R_to_use = std::fabs( R );
		auto func = [=, &R_to_use] (const flt_type & offset_R, const bool silent = true)
		{
			return SPCP(name)->offset_Sigma( R_to_use, offset_R, silent );
		};
		auto weight_func = [=, &group_c] (const flt_type & R, const bool silent = true)
		{
			if ( group_c == 0 )
			{
				return 2 * pi * R * SPCP(name)->proj_dens( R );
			}
			else
			{
				std::unique_ptr<name> group_profile(SPCP(name)->lensing_profile_extension_clone());
				group_profile->set_c(group_c);
				group_profile->set_tau(group_profile->tau()*group_c/SPCP(name)->c());
				return 2 * pi * R * group_profile->proj_dens(R);
			}
		};
		flt_type min_in_param( SMALL_FACTOR ), max_in_param( 2.5*SPCP(name)->rvir() ),
				out_param( 0 );
		out_param = brgastro::integrate_weighted_Romberg( &func, &weight_func,
				min_in_param, max_in_param, 0.00001, false, silent);
		return out_param;
	}
	flt_type _semiquick_group_Sigma( const flt_type & R,
			const flt_type & group_c, const bool silent = false ) const // As group_Delta_Sigma, but uses offset_Delta_Sigma cache to speed it up if overwritten
	{
		flt_type R_to_use = std::fabs( R );
		auto func = [=, &R] (const flt_type & offset_R, const bool silent = true)
		{
			return SPCP(name)->quick_offset_Sigma( R_to_use, offset_R, silent );
		};
		auto weight_func = [=, &group_c] (const flt_type & R, const bool silent = true)
		{
			if ( group_c == 0 )
			{
				return 2 * pi * R * SPCP(name)->proj_dens( R );
			}
			else
			{
				std::unique_ptr<name> group_profile(SPCP(name)->lensing_profile_extension_clone());
				group_profile->set_c(group_c);
				group_profile->set_tau(group_profile->tau()*group_c/SPCP(name)->c());
				return 2 * pi * R * group_profile->proj_dens(R);
			}
		};
		flt_type min_in_param( SMALL_FACTOR ), max_in_param( 2.5*SPCP(name)->rvir() ),
				out_param( 0 );
		out_param = brgastro::integrate_weighted_Romberg( &func, &weight_func,
				min_in_param, max_in_param, 0.00001, false);
		return out_param;
	}
	flt_type _two_halo_Sigma( const flt_type & R,
			const flt_type & group_c, const bool silent = false ) const
	{
		assert(false); // TODO: Fill this out

		return flt_type();
	}

	flt_type _shifted_Delta_Sigma( const flt_type & R,
			const bool silent = false ) const
	{
		flt_type R_to_use = std::fabs( R );
		flt_type sigma = SPCP(name)->_shift_sigma(R_to_use);

		brgastro::shifted_Delta_Sigma_functor<name> func( SPCP(name), R_to_use );
		brgastro::shifted_Delta_Sigma_weight_functor<name> weight_func( sigma );

		flt_type precision = 0.00001;

		flt_type min_in_params( 0 ), max_in_params( 4*sigma ),
				out_params( 0 );

		out_params = brgastro::integrate_weighted_Romberg( &func, &weight_func,
				min_in_params, max_in_params, precision, false );
		return out_params;
	}
	flt_type _semiquick_shifted_Delta_Sigma( const flt_type & R,
			const bool silent = false ) const
	{
		flt_type R_to_use = std::fabs( R );
		flt_type sigma = SPCP(name)->_shift_sigma(R_to_use);

		brgastro::shifted_Delta_Sigma_functor<name> func( SPCP(name), R_to_use );
		brgastro::shifted_Delta_Sigma_weight_functor<name> weight_func( sigma );

		flt_type min_in_params( 0 ), max_in_params( 4*sigma ),
				out_params( 0 );

		out_params = brgastro::integrate_weighted_Romberg( &func, &weight_func,
				min_in_params, max_in_params, 0.00001, false);
		return out_params;
	}

#endif // Advanced calculations which usually can't be overridden, but should be if possible

// Quick functions - should be overridden if a cache is implemented for the halo
#if(1)

	flt_type _quick_Delta_Sigma( const flt_type & R, const bool silent =
			false ) const
	{
		return SPCP(name)->_Delta_Sigma(R,silent);
	}
	flt_type _quick_offset_Delta_Sigma( const flt_type & R,
			const flt_type & offset_R, const bool silent = false ) const
	{
		return SPCP(name)->_offset_Delta_Sigma( R, offset_R, silent );
	}
	flt_type _quick_group_Delta_Sigma( const flt_type & R,
			const flt_type & group_c, const bool silent = false ) const
	{
		return SPCP(name)->_group_Delta_Sigma( R, group_c, silent );
	}
	flt_type _quick_two_halo_Delta_Sigma( const flt_type & R, const bool silent =
			false ) const
	{
		return SPCP(name)->_two_halo_Delta_Sigma(R,silent);
	}
	flt_type _quick_Sigma( const flt_type & R, const bool silent =
			false ) const
	{
		return SPCP(name)->_Sigma(R,silent);
	}
	flt_type _quick_offset_Sigma( const flt_type & R,
			const flt_type & offset_R, const bool silent = false ) const // As offset_Sigma, but uses cache to speed it up if overwritten
	{
		return SPCP(name)->_offset_Sigma( R, offset_R, silent );
	}
	flt_type _quick_group_Sigma( const flt_type & R,
			const flt_type & group_c, const bool silent = false ) const // As Sigma, but uses group_Delta_Sigma cache to speed it up if overwritten
	{
		return SPCP(name)->_group_Sigma( R, group_c, silent );
	}
	flt_type _quick_two_halo_Sigma( const flt_type & R, const bool silent =
			false ) const
	{
		return SPCP(name)->_two_halo_Sigma(R,silent);
	}
	flt_type _quick_shifted_Delta_Sigma( const flt_type & R,
			const bool silent = false ) const
	{
		return SPCP(name)->_shifted_Delta_Sigma( R, silent );
	}

#endif // Quick functions - should be overridden if a cache is implemented for the halo

// Simple calculation functions that should rarely need to be overridden
#if (1)
	flt_type _proj_enc_dens( const flt_type & R,
			const bool silent = false ) const // Mean surface density enclosed within a cylinder of radius R
	{
		flt_type R_to_use = max( std::fabs( R ), SMALL_FACTOR );
		return SPCP(name)->proj_enc_mass( R_to_use, silent )
				/ ( pi * square( R_to_use ) );
	}

	flt_type _shift_sigma( const flt_type & R, const bool silent = false ) const
	{
		return R*shift_factor(R,silent);
	}
	flt_type _shift_factor( const flt_type & R, const bool silent = false ) const
	{
		if(R==0) return 0;
		const flt_type theta_separation = afd(R,SPCP(name)->z());
		return shifting_cache().get(theta_separation,SPCP(name)->z())/theta_separation;
	}
	flt_type _Delta_Sigma( const flt_type & R, const bool silent =
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
	~lensing_profile_extension()
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
	flt_type proj_dens( const flt_type & R,
			const bool silent = false ) const // Projected surface density at radius R
	{
		return SPCP(name)->_proj_dens(R,silent);
	}
	flt_type proj_enc_mass( const flt_type & R,
			const bool silent = false ) const // Mass enclosed within a cylinder of radius R
	{
		return SPCP(name)->_proj_enc_mass(R,silent);
	}
#endif

	// These ones typically don't need to be overridden, but they can be if it's convenient
#if (1)
	flt_type proj_enc_dens( const flt_type & R,
			const bool silent = false ) const // Mean surface density enclosed within a cylinder of radius R
	{
		return SPCP(name)->_proj_enc_dens(R,silent);
	}
#endif

#endif // Projected mass and density functions

	// Weak lensing functions
#if (1)

	// Delta Sigma
#if (1)
	flt_type Delta_Sigma( const flt_type & R, const bool silent =
			false ) const // Weak lensing signal in tangential shear Delta-Sigma at radius R
	{
		return SPCP(name)->_Delta_Sigma(R,silent);
	}
	flt_type quick_Delta_Sigma( const flt_type & R,
			const bool silent = false ) const // As deltasigma, but uses cache to speed it up if overwritten
	{
		return SPCP(name)->_quick_Delta_Sigma( R, silent );
	}
#endif

	// Offset Delta Sigma
#if (1)
	flt_type offset_Delta_Sigma( const flt_type & R,
			const flt_type & offset_R, const bool silent = false ) const // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from position offset by offset_R
	{
		return SPCP(name)->_offset_Delta_Sigma( R, offset_R, silent );
	}
	flt_type quick_offset_Delta_Sigma( const flt_type & R,
			const flt_type & offset_R, const bool silent = false ) const // As offset_Delta_Sigma, but uses cache to speed it up if overwritten
	{
		return SPCP(name)->_quick_offset_Delta_Sigma( R, offset_R, silent );
	}
#endif

	// Group Delta Sigma
#if (1)
	flt_type group_Delta_Sigma( const flt_type & R,
			const flt_type & group_c=2.5, const bool silent = false ) const // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from ensemble of satellites in group with satellite concentration group_c
	{
		return SPCP(name)->_group_Delta_Sigma( R, group_c, silent );
	}
	flt_type semiquick_group_Delta_Sigma( const flt_type & R,
			const flt_type & group_c=2.5, const bool silent = false ) const // As group_Delta_Sigma, but uses offset_Delta_Sigma cache to speed it up if overwritten
	{
		return SPCP(name)->_semiquick_group_Delta_Sigma( R, group_c, silent );
	}
	flt_type quick_group_Delta_Sigma( const flt_type & R,
			const flt_type & group_c=2.5, const bool silent = false ) const // As deltasigma, but uses group_Delta_Sigma cache to speed it up if overwritten
	{
		return SPCP(name)->_quick_group_Delta_Sigma( R, group_c, silent );
	}
#endif

	// Two Halo Delta Sigma
#if (1)
	flt_type two_halo_Delta_Sigma( const flt_type & R, const bool silent = false ) const // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from ensemble of satellites in group with satellite concentration group_c
	{
		return SPCP(name)->_two_halo_Delta_Sigma( R, silent );
	}
	flt_type quick_two_halo_Delta_Sigma( const flt_type & R, const bool silent = false ) const // As deltasigma, but uses group_Delta_Sigma cache to speed it up if overwritten
	{
		return SPCP(name)->_quick_two_halo_Delta_Sigma( R, silent );
	}
#endif

	// Shifted Delta Sigma
#if (1)
	// This represents the weak-lensing signal after being corrected for errors due to relative
	// shifting of the lens and source due to an intervening mass distribution

	flt_type shift_factor( const flt_type & R, const bool silent = false ) const
	{
		return SPCP(name)->_shift_factor( R, silent );
	}

	flt_type shifted_Delta_Sigma( const flt_type & R,
			const bool silent = false ) const
	{
		return SPCP(name)->_shifted_Delta_Sigma( R, silent );
	}
	flt_type semiquick_shifted_Delta_Sigma( const flt_type & R,
			const bool silent = false ) const
	{
		return SPCP(name)->_semiquick_shifted_Delta_Sigma( R, silent );
	}
	flt_type quick_shifted_Delta_Sigma( const flt_type & R,
			const bool silent = false ) const
	{
		return SPCP(name)->_quick_shifted_Delta_Sigma( R, silent );
	}
#endif

	// Sigma
#if (1)
	flt_type Sigma( const flt_type & R, const bool silent =
			false ) const // Magnification signal Sigma at radius R
	{
		return SPCP(name)->_proj_dens(R,silent);
	}
	flt_type quick_Sigma( const flt_type & R,
			const bool silent = false ) const // As Sigma, but uses cache to speed it up if overridden
	{
		return SPCP(name)->_quick_Sigma( R, silent );
	}
#endif

	// Offset Sigma
#if (1)
	flt_type offset_Sigma( const flt_type & R,
			const flt_type & offset_R, const bool silent = false ) const // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from position offset by offset_R
	{
		return SPCP(name)->_offset_Sigma( R, offset_R, silent );
	}
	flt_type quick_offset_Sigma( const flt_type & R,
			const flt_type & offset_R, const bool silent = false ) const // As offset_Delta_Sigma, but uses cache to speed it up if overwritten
	{
		return SPCP(name)->_quick_offset_Sigma( R, offset_R, silent );
	}
#endif

	// Group Sigma
#if (1)
	flt_type group_Sigma( const flt_type & R,
			const flt_type & group_c=2.5, const bool silent = false ) const // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from ensemble of satellites in group with satellite concentration group_c
	{
		return SPCP(name)->_group_Sigma( R, group_c, silent );
	}
	flt_type semiquick_group_Sigma( const flt_type & R,
			const flt_type & group_c=2.5, const bool silent = false ) const // As group_Delta_Sigma, but uses offset_Delta_Sigma cache to speed it up if overwritten
	{
		return SPCP(name)->_semiquick_group_Sigma( R, group_c, silent );
	}
	flt_type quick_group_Sigma( const flt_type & R,
			const flt_type & group_c=2.5, const bool silent = false ) const // As deltasigma, but uses group_Delta_Sigma cache to speed it up if overwritten
	{
		return SPCP(name)->_quick_group_Sigma( R, group_c, silent );
	}
#endif

	// Two Halo Sigma
#if (1)
	flt_type two_halo_Sigma( const flt_type & R, const bool silent = false ) const // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from ensemble of satellites in group with satellite concentration group_c
	{
		return SPCP(name)->_two_halo_Delta_Sigma( R, silent );
	}
	flt_type quick_two_halo_Sigma( const flt_type & R, const bool silent = false ) const // As deltasigma, but uses group_Delta_Sigma cache to speed it up if overwritten
	{
		return SPCP(name)->_quick_two_halo_Delta_Sigma( R, silent );
	}
#endif

#include "brg_lensing/two_halo_term_methods.hpp"

#endif

};

} // end namespace brgastro

// Undef macros
#undef SPP
#undef SPCP

#endif /* _BRG_LENSING_PROFILE_EXTENSION_H_ */
