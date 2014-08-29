/**       @file density_profile_functors.cpp
 *
 *     Project: brg
 *        Path: /brg/density_profile/density_profile_functors.cpp
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#include <iostream>

#include "brg/brg_global.h"

#include "brg/brg_misc_functions.hpp"
#include "brg/brg_units.h"

#include "density_profile_functors.h"

// brgastro::accel_functor class methods
#if (1)

void brgastro::accel_functor::set_host_ptr(
		const density_profile *new_host )
{
	_host_ptr_ = new_host;
}

BRG_UNITS brgastro::accel_functor::operator()( CONST_BRG_UNITS_REF  in_param,
		const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		throw std::runtime_error("ERROR: Host must be assigned to accel_functor before function can be called.\n");
	}
	return _host_ptr_->accel( in_param, silent );
}

brgastro::accel_functor::accel_functor()
{
	_host_ptr_ = NULL;
}
brgastro::accel_functor::accel_functor( const density_profile *init_host )
{
	set_host_ptr( init_host );
}

#endif // end brgastro::accel_functor function implementations

// brgastro::spherical_density_functor class methods
#if (1)

void brgastro::spherical_density_functor::set_host_ptr(
		const density_profile *new_host )
{
	_host_ptr_ = new_host;
}

BRG_UNITS brgastro::spherical_density_functor::operator()(
		CONST_BRG_UNITS_REF  in_param, const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		throw std::runtime_error("ERROR: Host must be assigned to spherical_density_functor before function can be called.\n");
	}
	return 4 * pi * square( in_param )
			* _host_ptr_->dens( in_param );
}

brgastro::spherical_density_functor::spherical_density_functor()
{
	_host_ptr_ = NULL;
}
brgastro::spherical_density_functor::spherical_density_functor(
		const density_profile *init_host )
{
	set_host_ptr( init_host );
}

#endif // end brgastro::spherical_density_functor class methods

// brgastro::solve_rhm_functor class methods
#if (1)
void brgastro::solve_rhm_functor::set_host_ptr(
		const density_profile *new_host )
{
	_host_ptr_ = new_host;
}

void brgastro::solve_rhm_functor::set_target_mass(
		const BRG_MASS &new_target_mass )
{
	_target_mass_ = new_target_mass;
}

BRG_UNITS brgastro::solve_rhm_functor::operator()( CONST_BRG_UNITS_REF  in_param,
		const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		throw std::runtime_error("ERROR: Host must be assigned to solve_rhm_functor before function can be called.\n");
	}
	return std::fabs(_target_mass_ - _host_ptr_->enc_mass( fabs( in_param ), silent ) );
}

brgastro::solve_rhm_functor::solve_rhm_functor()
{
	_host_ptr_ = NULL;
	_target_mass_ = 0;
}
;
brgastro::solve_rhm_functor::solve_rhm_functor(
		const density_profile *init_host, const BRG_MASS &new_target_mass )
{
	set_host_ptr( init_host );
	set_target_mass( new_target_mass );
}
#endif // end brgastro::solve_rhm_functor function implementations
