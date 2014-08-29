/**       @file lensing_profile_extension_functors.cpp
 *
 *     Project: brg
 *        Path: /brg/lensing_profile_extension/lensing_profile_extension_functors.cpp
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "brg/brg_global.h"

#include "brg/brg_calculus.hpp"
#include "brg/brg_misc_functions.hpp"
#include "brg/brg_units.h"

#include "lensing_profile_extension_functors.h"

// brgastro::projected_density_functor class methods
#if (1)

void brgastro::projected_density_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
}

void brgastro::projected_density_functor::set_offset_R(
		CONST_BRG_DISTANCE_REF new_offset_R )
{
	_offset_R_ = new_offset_R;
}

BRG_UNITS brgastro::projected_density_functor::operator()(
		CONST_BRG_UNITS_REF in_param,
		const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		throw std::runtime_error("ERROR: Host must be assigned to projected_density_functor before function can be called.\n");
	}
	BRG_DISTANCE r = quad_add( in_param, _offset_R_ );
	return _host_ptr_->dens( r );
}

brgastro::projected_density_functor::projected_density_functor()
{
	_host_ptr_ = NULL;
	_offset_R_ = 0;
}

brgastro::projected_density_functor::projected_density_functor(
		const lensing_profile_extension *init_host, CONST_BRG_DISTANCE_REF init_offset_R )
{
	set_host_ptr( init_host );
	set_offset_R( init_offset_R );
}

#endif // end brgastro::projected_density_functor class methods

// brgastro::cylindrical_density_functor class methods
#if (1)

void brgastro::cylindrical_density_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
}

BRG_UNITS brgastro::cylindrical_density_functor::operator()(
		CONST_BRG_UNITS_REF in_param,const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		throw std::runtime_error("ERROR: Host must be assigned to cylindrical_density_functor before function can be called.\n");
	}
	return 2 * pi * in_param * _host_ptr_->proj_dens( in_param );
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

void brgastro::offset_ring_dens_functor::set_R0(
		CONST_BRG_DISTANCE_REF new_R0 )
{
	_R0_ = new_R0;
}

void brgastro::offset_ring_dens_functor::set_R(
		CONST_BRG_DISTANCE_REF new_R )
{
	_R_ = new_R;
}

BRG_UNITS brgastro::offset_ring_dens_functor::operator()(
		CONST_BRG_UNITS_REF in_param,
		const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		throw std::runtime_error("ERROR: Host must be assigned to offset_ring_dens_functor before function can be called.\n");
	}

	BRG_DISTANCE d = lc_add( _R0_, _R_, in_param );

	return _host_ptr_->proj_dens( d, silent );
}

void brgastro::offset_ring_dens_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
}

brgastro::offset_ring_dens_functor::offset_ring_dens_functor()
{
	_host_ptr_ = NULL;
	_R_ = 0;
	_R0_ = 0;
}

brgastro::offset_ring_dens_functor::offset_ring_dens_functor(
		const lensing_profile_extension *new_host, CONST_BRG_DISTANCE_REF new_R_0,
		CONST_BRG_DISTANCE_REF new_R )
{
	_host_ptr_ = new_host;
	_R_ = new_R;
	_R0_ = new_R_0;
}

#endif

// brgastro::offset_circ_dens_functor class methods
#if (1)

CONST_BRG_DISTANCE_REF brgastro::offset_circ_dens_functor::_arc_length_in_circle(
		CONST_BRG_DISTANCE_REF R2 ) const
{
	// Check for complete enclosure
	if( _R0_ + R2 <= _R_)
	{
		return 2.*pi*R2;
	}
	else
	{
		BRG_DISTANCE result = 2.*R2 * std::acos( (square(_R0_)-square(_R_)+square(R2)) / (2.*_R0_*R2) );
		if(result>0) return result;
		return 0.; // We'll get here only in cases due to round-off error, where it should actually be zero
	}
}

void brgastro::offset_circ_dens_functor::set_R0(
		CONST_BRG_DISTANCE_REF new_R_0 )
{
	_R0_ = new_R_0;
}

void brgastro::offset_circ_dens_functor::set_R(
		CONST_BRG_DISTANCE_REF new_R )
{
	_R_ = new_R;
}

BRG_UNITS brgastro::offset_circ_dens_functor::operator()(
		CONST_BRG_UNITS_REF in_param, const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		throw std::runtime_error("ERROR: Host must be assigned to offset_circ_dens_functor before function can be called.\n");
	}

	BRG_DISTANCE L = _arc_length_in_circle(in_param);

	return L * _host_ptr_->proj_dens( in_param, silent );
}

void brgastro::offset_circ_dens_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
}

brgastro::offset_circ_dens_functor::offset_circ_dens_functor()
{
	_host_ptr_ = NULL;
	_R0_ = 0;
	_R_ = 0;
}

brgastro::offset_circ_dens_functor::offset_circ_dens_functor(
		const lensing_profile_extension *new_host, CONST_BRG_DISTANCE_REF new_R_0,
		CONST_BRG_DISTANCE_REF new_R )
{
	_host_ptr_ = new_host;
	_R0_ = new_R_0;
	_R_ = new_R;
}
#endif // end brgastro::offset_circ_dens_functor function implementations

// brgastro::offset_WLsig_functor class methods
#if (1)

void brgastro::offset_WLsig_functor::set_R( CONST_BRG_DISTANCE_REF new_R )
{
	_R_ = new_R;
}

BRG_UNITS brgastro::offset_WLsig_functor::operator()(
		CONST_BRG_UNITS_REF in_param, const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		throw std::runtime_error("ERROR: Host must be assigned to offset_WLsig_functor before function can be called.\n");
	}

	return _host_ptr_->offset_WLsig( _R_, in_param, silent );
}

void brgastro::offset_WLsig_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
}

brgastro::offset_WLsig_functor::offset_WLsig_functor()
{
	_host_ptr_ = NULL;
	_R_ = 0;
}

brgastro::offset_WLsig_functor::offset_WLsig_functor(
		const lensing_profile_extension *new_host, CONST_BRG_DISTANCE_REF new_R )
{
	_host_ptr_ = new_host;
	_R_ = new_R;
}

#endif

// brgastro::quick_offset_WLsig_functor class methods
#if (1)

void brgastro::quick_offset_WLsig_functor::set_R( CONST_BRG_DISTANCE_REF new_R )
{
	_R_ = new_R;
}

BRG_UNITS brgastro::quick_offset_WLsig_functor::operator()(
		CONST_BRG_DISTANCE_REF in_param, const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		throw std::runtime_error("ERROR: Host must be assigned to offset_WLsig_functor before function can be called.\n");
	}

	return _host_ptr_->quick_offset_WLsig( _R_, in_param, silent );
}

void brgastro::quick_offset_WLsig_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
}

brgastro::quick_offset_WLsig_functor::quick_offset_WLsig_functor()
{
	_host_ptr_ = NULL;
	_R_ = 0;
}

brgastro::quick_offset_WLsig_functor::quick_offset_WLsig_functor(
		const lensing_profile_extension *new_host, CONST_BRG_DISTANCE_REF new_R )
{
	_host_ptr_ = new_host;
	_R_ = new_R;
}

#endif

// brgastro::group_WLsig_weight_functor class methods
#if (1)

void brgastro::group_WLsig_weight_functor::set_c( const double new_c )
{
	_c_ = new_c;
}

BRG_UNITS brgastro::group_WLsig_weight_functor::operator()(
		CONST_BRG_DISTANCE_REF in_param, const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		throw std::runtime_error("ERROR: Host must be assigned to offset_WLsig_functor before function can be called.\n");
	}

	if ( _c_ == 0 )
	{
		return 2 * pi * in_param * _host_ptr_->proj_dens( in_param );
	}
	else
	{
		BRG_UNIQUE_PTR<brgastro::lensing_profile_extension> group_profile(_host_ptr_->lensing_profile_extension_clone());
		group_profile->set_c(_c_);
		group_profile->set_tau(group_profile->tau()*_c_/_host_ptr_->c());
		return 2 * pi * in_param * group_profile->proj_dens(in_param);
	}
}

void brgastro::group_WLsig_weight_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
}

brgastro::group_WLsig_weight_functor::group_WLsig_weight_functor()
{
	_host_ptr_ = NULL;
	_c_ = 0;
}

brgastro::group_WLsig_weight_functor::group_WLsig_weight_functor(
		const lensing_profile_extension *new_host, const double init_c )
{
	_host_ptr_ = new_host;
	_c_ = init_c;
}

#endif

// brgastro::shifted_WLsig_weight_functor class methods
#if (1)

BRG_UNITS brgastro::shifted_WLsig_weight_functor::operator()(
		CONST_BRG_UNITS_REF in_param, const bool silent ) const
{
	// The output here is the height of a Rayleigh distribution at in_param
	return in_param/square(_sigma_) * std::exp(-square(in_param)/(2*square(_sigma_)));
}

void brgastro::shifted_WLsig_weight_functor::set_sigma(
		CONST_BRG_DISTANCE_REF new_sigma )
{
	_sigma_ = new_sigma;
}

brgastro::shifted_WLsig_weight_functor::shifted_WLsig_weight_functor()
{
	_sigma_ = 0;
}

brgastro::shifted_WLsig_weight_functor::shifted_WLsig_weight_functor(
		CONST_BRG_DISTANCE_REF new_sigma )
{
	_sigma_ = new_sigma;
}

#endif

// brgastro::shifted_WLsig_circ_functor class methods
#if (1)

BRG_UNITS brgastro::shifted_WLsig_circ_functor::operator()(
		CONST_BRG_UNITS_REF in_param, const bool silent ) const
{
	// in_param here will be angle theta in radians

	const BRG_DISTANCE R_actual(lc_add(_R_, _R_shift_, in_param));

	const BRG_ANGLE theta(asin(_R_shift_/R_actual * sin(in_param)));
	const double angle_factor = cos(theta);

	double extra_shear_factor;
	if(_R_shift_==0)
		extra_shear_factor = 0;
	else
		extra_shear_factor = (R_actual-_R_)/_R_shift_*
			_host_ptr_->shift_factor(0);

	if(isbad(theta))
	{
		return _host_ptr_->WLsig(R_actual,silent) +
			extra_shear_factor*sigma_crit(_host_ptr_->z(),2*_host_ptr_->z());
	}
	else
	{
		return _host_ptr_->WLsig(R_actual,silent)*angle_factor +
				extra_shear_factor*sigma_crit(_host_ptr_->z(),2*_host_ptr_->z()); //TODO Should there be a factor of 1/2 here?
	}
}

void brgastro::shifted_WLsig_circ_functor::set_R_shift(
		CONST_BRG_DISTANCE_REF new_R_shift )
{
	_R_shift_ = new_R_shift;
}

void brgastro::shifted_WLsig_circ_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
}

void brgastro::shifted_WLsig_circ_functor::set_R(
		CONST_BRG_DISTANCE_REF new_R )
{
	_R_ = new_R;
}

brgastro::shifted_WLsig_circ_functor::shifted_WLsig_circ_functor( const lensing_profile_extension *new_host,
		CONST_BRG_DISTANCE_REF new_R_shift, CONST_BRG_DISTANCE_REF new_R )
{
	_host_ptr_ = new_host;
	_R_shift_ = new_R_shift;
	_R_ = new_R;
}

#endif

// brgastro::shifted_WLsig_functor class methods
#if (1)

BRG_UNITS brgastro::shifted_WLsig_functor::operator()(
		CONST_BRG_UNITS_REF in_param,const bool silent ) const
{
	// in_param here will be R_shift
	shifted_WLsig_circ_functor func(_host_ptr_,in_param,_R_);

	const double min_in_param = 0;
	const double max_in_param = pi;
	const double precision = 0.000001;
	BRG_UNITS out_param(0);

	out_param = brgastro::integrate_Romberg( &func, min_in_param, max_in_param,
			precision, false, silent );

	return out_param /= pi;
}

void brgastro::shifted_WLsig_functor::set_host_ptr(
		const lensing_profile_extension *new_host )
{
	_host_ptr_ = new_host;
}

void brgastro::shifted_WLsig_functor::set_R(
		CONST_BRG_DISTANCE_REF new_R )
{
	_R_ = new_R;
}

brgastro::shifted_WLsig_functor::shifted_WLsig_functor( const lensing_profile_extension *new_host,
		CONST_BRG_DISTANCE_REF new_R)
{
	_host_ptr_ = new_host;
	_R_ = new_R;
}

#endif
