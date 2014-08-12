/**       @file lensing_profile_extension_functors.cpp
 *
 *     Project: brg
 *        Path: /brg/lensing_profile_extension/lensing_profile_extension_functors.cpp
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#include <iostream>

#include "../../brg_global.h"

#include "../../brg_misc_functions.hpp"
#include "../../brg_units.h"
#include "lensing_profile_extension_functors.h"

// brgastro::projected_density_functor class methods
#if (1)

const int brgastro::projected_density_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

const int brgastro::projected_density_functor::set_offset_R(
		const BRG_DISTANCE &new_offset_R )
{
	_offset_R_ = new_offset_R;
	return 0;
}

const int brgastro::projected_density_functor::operator()(
		const BRG_UNITS & in_param,
		BRG_UNITS & out_param, const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to projected_density_functor before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}
	BRG_DISTANCE r = quad_add( in_param, _offset_R_ );
	out_param = _host_ptr_->dens( r );
	return 0;
}

brgastro::projected_density_functor::projected_density_functor()
{
	_host_ptr_ = NULL;
	_offset_R_ = 0;
}

brgastro::projected_density_functor::projected_density_functor(
		const lensing_profile_extension *init_host, const BRG_DISTANCE &init_offset_R )
{
	set_host_ptr( init_host );
	set_offset_R( init_offset_R );
}

#endif // end brgastro::projected_density_functor class methods

// brgastro::cylindrical_density_functor class methods
#if (1)

const int brgastro::cylindrical_density_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

const int brgastro::cylindrical_density_functor::operator()(
		const BRG_UNITS & in_param,
		BRG_UNITS & out_param, const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to cylindrical_density_functor before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}
	out_param = 2 * pi * in_param * _host_ptr_->proj_dens( in_param );
	return 0;
}

brgastro::cylindrical_density_functor::cylindrical_density_functor()
{
	_host_ptr_ = NULL;
}
brgastro::cylindrical_density_functor::cylindrical_density_functor(
		const lensing_profile_extension *init_host )
{
	set_host_ptr( init_host );
}

#endif // end brgastro::cylindrical_density_functor class methods

// brgastro::offset_ring_dens_functor class methods
#if (1)

const int brgastro::offset_ring_dens_functor::set_R0(
		const BRG_DISTANCE &new_R0 )
{
	_R0_ = new_R0;
	return 0;
}

const int brgastro::offset_ring_dens_functor::set_R(
		const BRG_DISTANCE &new_R )
{
	_R_ = new_R;
	return 0;
}

const int brgastro::offset_ring_dens_functor::operator()(
		const BRG_UNITS &in_param,
		BRG_UNITS & out_param, const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to offset_ring_dens_functor before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}

	BRG_DISTANCE d = lc_add( _R0_, _R_, in_param );

	out_param = _host_ptr_->proj_dens( d, silent );

	return 0;
}

const int brgastro::offset_ring_dens_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

brgastro::offset_ring_dens_functor::offset_ring_dens_functor()
{
	_host_ptr_ = NULL;
	_R_ = 0;
	_R0_ = 0;
}

brgastro::offset_ring_dens_functor::offset_ring_dens_functor(
		const lensing_profile_extension *new_host, const BRG_DISTANCE &new_R_0,
		const BRG_DISTANCE &new_R )
{
	_host_ptr_ = new_host;
	_R_ = new_R;
	_R0_ = new_R_0;
}

#endif

// brgastro::offset_circ_dens_functor class methods
#if (1)

const int brgastro::offset_circ_dens_functor::set_R0(
		const BRG_DISTANCE &new_R_0 )
{
	_R0_ = new_R_0;
	return 0;
}

const int brgastro::offset_circ_dens_functor::operator()(
		const std::vector< BRG_UNITS > &in_params,
		std::vector< BRG_UNITS > & out_params, const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to offset_circ_dens_functor before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}
	if ( in_params.size() != 2 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: in_params.size() must == 2 in offset_circ_dens_functor.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	out_params.clear();

	BRG_DISTANCE R = in_params[0];
	BRG_DISTANCE d = lc_add( _R0_, R, in_params[1] );

	BRG_UNITS result = _host_ptr_->proj_dens( d, silent );
	out_params.push_back( result );

	return 0;
}

const int brgastro::offset_circ_dens_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

brgastro::offset_circ_dens_functor::offset_circ_dens_functor()
{
	_host_ptr_ = NULL;
	_R0_ = 0;
}

brgastro::offset_circ_dens_functor::offset_circ_dens_functor(
		const lensing_profile_extension *new_host, const BRG_DISTANCE &new_R_0 )
{
	_host_ptr_ = new_host;
	_R0_ = new_R_0;
}
#endif // end brgastro::offset_circ_dens_functor function implementations

// brgastro::offset_WLsig_functor class methods
#if (1)

const int brgastro::offset_WLsig_functor::set_R( const BRG_DISTANCE &new_R )
{
	_R_ = new_R;
	return 0;
}

const int brgastro::offset_WLsig_functor::operator()(
		const BRG_UNITS &in_param,
		BRG_UNITS & out_param, const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to offset_WLsig_functor before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}

	out_param = _host_ptr_->offset_WLsig( _R_, in_param, silent );

	return 0;
}

const int brgastro::offset_WLsig_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

brgastro::offset_WLsig_functor::offset_WLsig_functor()
{
	_host_ptr_ = NULL;
	_R_ = 0;
}

brgastro::offset_WLsig_functor::offset_WLsig_functor(
		const lensing_profile_extension *new_host, const BRG_DISTANCE &new_R )
{
	_host_ptr_ = new_host;
	_R_ = new_R;
}

#endif

// brgastro::quick_offset_WLsig_functor class methods
#if (1)

const int brgastro::quick_offset_WLsig_functor::set_R( const BRG_DISTANCE &new_R )
{
	_R_ = new_R;
	return 0;
}

const int brgastro::quick_offset_WLsig_functor::operator()(
		const BRG_UNITS &in_param,
		BRG_UNITS & out_param, const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to offset_WLsig_functor before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}

	out_param = _host_ptr_->quick_offset_WLsig( _R_, in_param, silent );

	return 0;
}

const int brgastro::quick_offset_WLsig_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

brgastro::quick_offset_WLsig_functor::quick_offset_WLsig_functor()
{
	_host_ptr_ = NULL;
	_R_ = 0;
}

brgastro::quick_offset_WLsig_functor::quick_offset_WLsig_functor(
		const lensing_profile_extension *new_host, const BRG_DISTANCE &new_R )
{
	_host_ptr_ = new_host;
	_R_ = new_R;
}

#endif

// brgastro::offset_WLsig_weight_functor class methods
#if (1)

const int brgastro::offset_WLsig_weight_functor::set_c( const double new_c )
{
	_c_ = new_c;
	return 0;
}

const int brgastro::offset_WLsig_weight_functor::operator()(
		const BRG_UNITS &in_param,
		BRG_UNITS & out_param, const bool silent ) const
{
	BRG_UNITS result;
	if ( _host_ptr_ == NULL )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to offset_WLsig_functor before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}

	if ( _c_ == 0 )
	{
		result = 2 * pi * in_param * _host_ptr_->proj_dens( in_param );
	}
	else
	{
		BRG_UNIQUE_PTR<brgastro::lensing_profile_extension> group_profile(_host_ptr_->lensing_profile_extension_clone());
		group_profile->set_c(_c_);
		result = 2 * pi * in_param * group_profile->proj_dens(in_param);
	}

	out_param = result;

	return 0;
}

const int brgastro::offset_WLsig_weight_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

brgastro::offset_WLsig_weight_functor::offset_WLsig_weight_functor()
{
	_host_ptr_ = NULL;
	_c_ = 0;
}

brgastro::offset_WLsig_weight_functor::offset_WLsig_weight_functor(
		const lensing_profile_extension *new_host, const double init_c )
{
	_host_ptr_ = new_host;
	_c_ = init_c;
}

#endif
