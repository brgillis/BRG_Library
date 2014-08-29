/**       @file interpolator_derivative_functors.cpp
 *
 *     Project: brg
 *        Path: /brg/interpolator/interpolator_derivative_functors.cpp
 *
 *  Created on: 29 Aug 2014
 *      Author: brg
 */

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <utility>

#include "brg/math/calculus/calculus.hpp"
#include "brg/math/random/distributions.hpp"

#include "interpolator_derivative_functors.h"

// brgastro::interpolator_functor class method implementations
#if (1)

// Swap functions
void brgastro::interpolator_functor::swap(interpolator_functor & other)
{
	std::swap(_interpolator_ptr_,other._interpolator_ptr_);
	std::swap(_interpolator_ptr_set_up_,other._interpolator_ptr_set_up_);
}

// Constructors
brgastro::interpolator_functor::interpolator_functor()
{
	_interpolator_ptr_ = NULL;
	_interpolator_ptr_set_up_ = false;
}
brgastro::interpolator_functor::interpolator_functor(
		const interpolator_functor& other)
{
	_interpolator_ptr_ = other._interpolator_ptr_;
	_interpolator_ptr_set_up_ = other._interpolator_ptr_set_up_;
}
brgastro::interpolator_functor::interpolator_functor(
		const brgastro::interpolator *init_interpolator_ptr )
{
	set_interpolator_ptr( init_interpolator_ptr );
}

// Operator=
brgastro::interpolator_functor & brgastro::interpolator_functor::operator=(
		brgastro::interpolator_functor other)
{
    swap(other);

    return *this;
}

// Set functions
const int brgastro::interpolator_functor::set_interpolator_ptr(
		const brgastro::interpolator *new_spline_ptr )
{
	_interpolator_ptr_ = new_spline_ptr;
	_interpolator_ptr_set_up_ = true;
	return 0;
}

// Function method
BRG_UNITS brgastro::interpolator_functor::operator()( CONST_BRG_UNITS_REF  in_param,
 const bool silent ) const
{
	if ( !_interpolator_ptr_set_up_ )
	{
		throw std::logic_error("ERROR: Interpolator pointer must be defined in spline_functor.\n");
	}

	return ( *_interpolator_ptr_ )( in_param );
	return 0;
}
#endif

// brgastro::interpolator_derivative_functor class method implementations
#if (1)

// Swap functions
void brgastro::interpolator_derivative_functor::swap(
		interpolator_derivative_functor& other)
{
    using std::swap;
	swap(_interpolator_functor_set_up_,other._interpolator_functor_set_up_);
	_interpolator_functor_.swap(other._interpolator_functor_);
}
namespace std
{
	template <>
	void swap(brgastro::interpolator_derivative_functor &same,
			brgastro::interpolator_derivative_functor &other)
	{
		same.swap(other);
	}
}

// Constructors
brgastro::interpolator_derivative_functor::interpolator_derivative_functor()
{
	_interpolator_functor_set_up_ = false;
}
brgastro::interpolator_derivative_functor::interpolator_derivative_functor(
		const interpolator_derivative_functor& other)
{
	_interpolator_functor_set_up_ = other._interpolator_functor_set_up_;
	_interpolator_functor_ = other._interpolator_functor_;
}
brgastro::interpolator_derivative_functor::interpolator_derivative_functor(
		brgastro::interpolator *init_spline_ptr )
{
	set_interpolator_ptr( init_spline_ptr );
}

// Operator=
brgastro::interpolator_derivative_functor& brgastro::interpolator_derivative_functor::operator=(
		brgastro::interpolator_derivative_functor other)
{
	swap(other);
	return *this;
}

// Set functions
const int brgastro::interpolator_derivative_functor::set_interpolator_ptr(
		const brgastro::interpolator *new_interpolator_ptr )
{
	_interpolator_functor_.set_interpolator_ptr( new_interpolator_ptr );
	_interpolator_functor_set_up_ = true;
	return 0;
}

// Function method
BRG_UNITS brgastro::interpolator_derivative_functor::operator()(
		CONST_BRG_UNITS_REF  in_param, const bool silent ) const
{
	if ( !_interpolator_functor_set_up_ )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Spline function must be set up in spline_derivative_functor.\n";
		return NOT_SET_UP_ERROR;
	}
	BRG_UNITS Jacobian=0;

	Jacobian = differentiate( &_interpolator_functor_, in_param );

	return Jacobian;
}

#endif

// brgastro::interpolator_derivative_weight_functor class method implementations
#if (1)
// Swap functions
void brgastro::interpolator_derivative_weight_functor::swap(
		interpolator_derivative_weight_functor &other)
{
	using std::swap;
	swap(_sample_scale_,other._sample_scale_);
	swap(_sample_max_width_,other._sample_max_width_);
	swap(_t_max_,other._t_max_);
	swap(_t_min_,other._t_min_);
	swap(_centre_point_,other._centre_point_);

}
namespace std
{
	template <>
	void swap(brgastro::interpolator_derivative_weight_functor &same,
			brgastro::interpolator_derivative_weight_functor &other)
	{
		same.swap(other);
	}
}


// Constructors
brgastro::interpolator_derivative_weight_functor::interpolator_derivative_weight_functor()
{
	_sample_scale_ = 0;
	_sample_max_width_ = 0;
	_t_max_ = -( 0.99 * DBL_MIN );
	_t_min_ = DBL_MAX;
	_centre_point_ = 0;
}
brgastro::interpolator_derivative_weight_functor::interpolator_derivative_weight_functor(
		const interpolator_derivative_weight_functor &other)
{
	_sample_scale_ = other._sample_scale_;
	_sample_max_width_ = other._sample_max_width_;
	_t_max_ = other._t_max_;
	_t_min_ = other._t_min_;
	_centre_point_ = other._centre_point_;
}

// Operator=
brgastro::interpolator_derivative_weight_functor &
	brgastro::interpolator_derivative_weight_functor::operator=(
			brgastro::interpolator_derivative_weight_functor other)
{
	swap(other);
	return *this;
}

// Set functions
const int brgastro::interpolator_derivative_weight_functor::set_sample_scale(
		double new_sample_scale )
{
	_sample_scale_ = new_sample_scale;
	return 0;
}

const int brgastro::interpolator_derivative_weight_functor::set_sample_max_width(
		double new_sample_max_width )
{
	_sample_max_width_ = new_sample_max_width;
	return 0;
}

const int brgastro::interpolator_derivative_weight_functor::set_center_point(
		double new_center_point )
{
	_centre_point_ = new_center_point;
	return 0;
}

const int brgastro::interpolator_derivative_weight_functor::set_t_min(
		double new_t_min )
{
	_t_min_ = new_t_min;
	return 0;
}

const int brgastro::interpolator_derivative_weight_functor::set_t_max(
		double new_t_max )
{
	_t_max_ = new_t_max;
	return 0;
}

// Function method
BRG_UNITS brgastro::interpolator_derivative_weight_functor::operator()(
		CONST_BRG_UNITS_REF in_param, const bool silent ) const
{

	// Check bounds
	if ( std::fabs( in_param - _centre_point_ )
			> _sample_max_width_ * std::fabs( _t_max_ - _t_min_ ) )
	{
		return 0;
	}
	else
	{
		return Gaus_pdf( in_param, _centre_point_,
				_sample_scale_ * std::fabs( _t_max_ - _t_min_ ) );
	}
}
#endif


