/**********************************************************************\
  @file interpolator_derivative.cpp

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

#include <utility>

#include "brg/math/calculus/differentiate.hpp"
#include "brg/math/calculus/integrate.hpp"
#include "brg/math/random/distributions.hpp"
#include "brg/utility.hpp"

#include "interpolator_derivative.h"

// interpolator_derivative static member initialisations
#if (1)
double brgastro::interpolator_derivative::_default_sample_scale_ = 0.02; // As a fraction of t_max-t_min
double brgastro::interpolator_derivative::_default_sample_max_width_ = 0.05; // As a fraction of t_max-t_min
double brgastro::interpolator_derivative::_default_sample_precision_ = 0.01;
#endif

// brgastro::interpolator_derivative class method implementations
#if (1)

// Swap functions
void brgastro::interpolator_derivative::swap(interpolator_derivative &other)
{
	using std::swap;
	swap(_interpolation_type_,other._interpolation_type_);
	swap(_interpolator_ptr_,other._interpolator_ptr_);
	swap(_interpolator_ptr_set_up_,other._interpolator_ptr_set_up_);

	_known_interpolator_.swap(other._known_interpolator_);
	swap(_unknown_t_list_,other._unknown_t_list_);

	swap(_t_max_,other._t_max_);
	swap(_t_min_,other._t_min_);

	swap(_calculated_,other._calculated_);

	swap(_sample_scale_,other._sample_scale_);
	swap(_sample_max_width_,other._sample_max_width_);
	swap(_sample_precision_,other._sample_precision_);
}
namespace std
{
	template <>
	void swap(brgastro::interpolator_derivative &same,
			brgastro::interpolator_derivative &other)
	{
		same.swap(other);
	}
}

// Constructors
brgastro::interpolator_derivative::interpolator_derivative()
{
	clear();
}
brgastro::interpolator_derivative::interpolator_derivative( const interpolator_derivative &other )
{
	_interpolation_type_ = other._interpolation_type_;
	_interpolator_ptr_ = other._interpolator_ptr_;
	_interpolator_ptr_set_up_ = other._interpolator_ptr_set_up_;

	_known_interpolator_ = other._known_interpolator_;
	_unknown_t_list_ = other._unknown_t_list_;

	_t_max_ = other._t_max_;
	_t_min_ = other._t_min_;

	_calculated_ = other._calculated_;

	_sample_scale_ = other._sample_scale_;
	_sample_max_width_ = other._sample_max_width_;
	_sample_precision_ = other._sample_precision_;
}
brgastro::interpolator_derivative::interpolator_derivative(
		brgastro::interpolator *init_spline_ptr )
{
	clear();
	set_spline_ptr( init_spline_ptr );
}

// Operator=
brgastro::interpolator_derivative & brgastro::interpolator_derivative::operator=(
		interpolator_derivative other)
{
	swap(other);
	return *this;
}

// Set functions
void brgastro::interpolator_derivative::set_spline_ptr(
		brgastro::interpolator *new_spline_ptr )
{
	_interpolator_ptr_ = new_spline_ptr;
	_interpolator_ptr_set_up_ = true;
}

void brgastro::interpolator_derivative::clear_spline_ptr()
{
	_interpolator_ptr_ = NULL;
	_interpolator_ptr_set_up_ = false;
}

void brgastro::interpolator_derivative::set_default_sample_scale(
		double new_default_sample_scale )
{
	_default_sample_scale_ = new_default_sample_scale;
}

void brgastro::interpolator_derivative::set_default_sample_max_width(
		double new_default_sample_max_width )
{
	_default_sample_max_width_ = new_default_sample_max_width;
}

void brgastro::interpolator_derivative::set_sample_scale(
		double new_sample_scale )
{
	_sample_scale_ = new_sample_scale;
	_calculated_ = false;
}

void brgastro::interpolator_derivative::set_sample_max_width(
		double new_sample_max_width )
{
	_sample_max_width_ = new_sample_max_width;
	_calculated_ = false;
}

void brgastro::interpolator_derivative::reset_sample_scale() // Sets it to default
{
	_sample_scale_ = _default_sample_scale_;
	_calculated_ = false;
}

void brgastro::interpolator_derivative::reset_sample_max_width() // Sets it to default
{
	_sample_max_width_ = _default_sample_max_width_;
	_calculated_ = false;
}

void brgastro::interpolator_derivative::set_interpolation_type(
		brgastro::interpolator::allowed_interpolation_type new_type)
{
	_known_interpolator_.set_interpolation_type(new_type);
	_interpolation_type_ = new_type;
	_calculated_ = false;
}

void brgastro::interpolator_derivative::reset_interpolation_type()
{
	_known_interpolator_.set_interpolation_type(_known_interpolator_.default_interpolation_type());
	_interpolation_type_ = _known_interpolator_.default_interpolation_type();
	_calculated_ = false;
}

void brgastro::interpolator_derivative::clear_known_points()
{
	_known_interpolator_.clear();
	if(_unknown_t_list_.empty())
	{
		_t_max_ = ( -std::numeric_limits<double>::max() );
		_t_min_ = std::numeric_limits<double>::max();
	}
	_calculated_ = false;
}

void brgastro::interpolator_derivative::clear_unknown_points()
{
	_unknown_t_list_.clear();
	if(_known_interpolator_.empty())
	{
		_t_max_ = ( -std::numeric_limits<double>::max() );
		_t_min_ = std::numeric_limits<double>::max();
	}
	_calculated_ = false;
}

void brgastro::interpolator_derivative::clear_points()
{
	_known_interpolator_.clear();
	_unknown_t_list_.clear();
	_calculated_ = false;
}

void brgastro::interpolator_derivative::clear()
{
	_interpolation_type_ = _known_interpolator_.default_interpolation_type();
	_interpolator_ptr_ = 0;
	_interpolator_ptr_set_up_ = false;

	clear_points();

	_sample_scale_ = _default_sample_scale_;
	_sample_max_width_ = _default_sample_max_width_;
	_sample_precision_ = _default_sample_precision_;
}

void brgastro::interpolator_derivative::add_point( const double t,
		const double x )
{
	_known_interpolator_.add_point( t, x );
	_calculated_ = false;
}

void brgastro::interpolator_derivative::add_unknown_point( const double t )
{
	_unknown_t_list_.push_back( t );
	_calculated_ = false;
}

// Get functions
double brgastro::interpolator_derivative::operator()( double xval, bool silent ) const
{
	if ( !_interpolator_ptr_set_up_ )
	{
		if ( ssize(_known_interpolator_) >= 2 ) // We can use the known spline for everything
		{
			return _known_interpolator_( xval );
		} // if(ssize(known_spline) >= 2)
		else // We don't know enough to get any points
		{
			throw std::logic_error("ERROR: Spline_derivative called without spline assigned to it.\n");
		} //  if(ssize(known_spline) >= 2) ... else
	} // if(!spline_set_up)

	if ( _calculated_ )
	{
		return _estimated_interpolator_( xval );
	} // if(calculated)
	else // We'll have to calculate
	{
		// Get t_min and t_max
		_t_min_ = std::numeric_limits<double>::max();
		_t_max_ = ( -std::numeric_limits<double>::max() );

		for ( ssize_t i = 0; i < ssize(_known_interpolator_.sorted_data()); i++ )
		{
			if ( _known_interpolator_.sorted_data().at( i ).first < _t_min_ )
				_t_min_ = _known_interpolator_.sorted_data().at( i ).first;
			if ( _known_interpolator_.sorted_data().at( i ).first > _t_max_ )
				_t_max_ = _known_interpolator_.sorted_data().at( i ).first;
		}

		for ( ssize_t i = 0, size=ssize(_unknown_t_list_); i < size; i++ )
		{
			if ( _unknown_t_list_[i] < _t_min_ )
				_t_min_ = _unknown_t_list_[i];
			if ( _unknown_t_list_[i] > _t_max_ )
				_t_max_ = _unknown_t_list_[i];
		}

		// Set up the estimated spline, starting by making a copy of the known spline
		_estimated_interpolator_ = _known_interpolator_;
		_estimated_interpolator_.set_interpolation_type(_interpolation_type_);
		ssize_t num_points_to_calculate = ssize(_unknown_t_list_);

		auto interpolator_functor = [&] (const double & in_param, bool silent=false)
		{
			return (*_interpolator_ptr_)(in_param);
		};

		auto spline_derivative_functor_val = [&] (const double & in_param, bool silent=false)
		{
			return differentiate( &interpolator_functor, in_param );
		};

		double t;

		auto spline_derivative_weight_functor_val = [&] (const double & in_param, bool silent=false)
		{

			if ( std::fabs( in_param - t )
					> _sample_max_width_ * std::fabs( _t_max_ - _t_min_ ) )
			{
				return 0.;
			}
			else
			{
				return Gaus_pdf( in_param, t,
						_sample_scale_ * std::fabs( _t_max_ - _t_min_ ) );
			}
		};

		double delta_t = fabs( _t_max_ - _t_min_ ) * _sample_max_width_;

		for ( ssize_t i = 0; i < num_points_to_calculate; i++ ) // For each point we need to calculate
		{
			t = _unknown_t_list_[i];
			BRG_UNITS min_in_params( t - delta_t );
			BRG_UNITS max_in_params( t + delta_t );
			BRG_UNITS out_params( 0 );
			BRG_UNITS Jacobian( 0 );

			if ( delta_t <= 0 )
			{
				Jacobian = differentiate( &spline_derivative_functor_val, min_in_params);

				_estimated_interpolator_.add_point( t, Jacobian );
			}
			else
			{

				out_params = integrate_weighted_Romberg(
						&spline_derivative_functor_val,
						&spline_derivative_weight_functor_val,
						min_in_params, max_in_params, _sample_precision_,
						false );

				_estimated_interpolator_.add_point( t, out_params );

			}

		} // for(ssize_t i = 0; i < num_points_to_calculate; i++ )
	} // if(calculated) ... else

	_calculated_ = true;

	return _estimated_interpolator_( xval );
} // const double operator()(double xval)

#endif
