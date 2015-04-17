/**********************************************************************\
  @file name_functors.h

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

// body file: name_functors.cpp

#ifndef _BRG_LENSING_PROFILE_EXTENSION_FUNCTORS_H_
#define _BRG_LENSING_PROFILE_EXTENSION_FUNCTORS_H_

#include "brg/global.h"

#include "brg/math/calculus/integrate.hpp"

#include "brg_physics/astro.h"
#include "brg_physics/units/unit_conversions.hpp"

namespace brgastro {

template<typename name>
class projected_density_functor
{
	/**********************************
	 projected_density_functor class
	 -----------------------------

	 Function class integrating density along a projected line

	 Parent class: function_class (from brg_functors)

	 **********************************/
private:
	const name *_host_ptr_;
	BRG_UNITS _offset_R_;

public:

	void set_host_ptr( const name *new_host_ptr )
	{
		_host_ptr_ = new_host_ptr;
	}
	const name * host_ptr()
	{
		return _host_ptr_;
	}

	void set_offset_R( CONST_BRG_DISTANCE_REF new_offset_R )
	{
		_offset_R_ = new_offset_R;
	}
	const BRG_DISTANCE offset_R()
	{
		return _offset_R_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF in_param,
			const bool silent = false ) const
	{
		if ( _host_ptr_ == NULL )
		{
			throw std::runtime_error("ERROR: Host must be assigned to projected_density_functor before function can be called.\n");
		}
		BRG_DISTANCE r = quad_add( in_param, _offset_R_ );
		return _host_ptr_->dens( r );
	}

	projected_density_functor( const name *init_host=NULL,
			CONST_BRG_DISTANCE_REF init_offset_R=0 )
	{
		set_host_ptr( init_host );
		set_offset_R( init_offset_R );
	}
	virtual ~projected_density_functor()
	{
	}
};

template<typename name>
class cylindrical_density_functor
{
	/**********************************
	 cylindrical_density_functor class
	 -----------------------------

	 Function class for integrating density in a cylinder

	 Parent class: function_class (from brg_functors)

	 **********************************/
private:
	const name *_host_ptr_;

public:

	void set_host_ptr( const name *new_host_ptr )
	{
		_host_ptr_ = new_host_ptr;
	}
	const name * host_ptr()
	{
		return _host_ptr_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF in_param, const bool silent = false ) const
	{
		if ( _host_ptr_ == NULL )
		{
			throw std::runtime_error("ERROR: Host must be assigned to cylindrical_density_functor before function can be called.\n");
		}
		return 2 * pi * in_param * _host_ptr_->proj_dens( in_param );
	}

	cylindrical_density_functor( const name *init_host = NULL )
	{
		set_host_ptr( init_host );
	}
	virtual ~cylindrical_density_functor()
	{
	}
};

template<typename name>
class offset_ring_dens_functor
{

	const name *_host_ptr_;
	BRG_DISTANCE _R0_, _R_;

public:

	void set_host_ptr( const name *new_host_ptr )
	{
		_host_ptr_ = new_host_ptr;
	}
	const name * host_ptr()
	{
		return _host_ptr_;
	}

	void set_R0( CONST_BRG_DISTANCE_REF new_R0 )
	{
		_R0_ = new_R0;
	}
	CONST_BRG_DISTANCE_REF  R0()
	{
		return _R0_;
	}

	void set_R( CONST_BRG_DISTANCE_REF new_R )
	{
		_R_ = new_R;
	}
	CONST_BRG_DISTANCE_REF  R()
	{
		return _R_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param,
			const bool silent = false ) const
	{
		if ( _host_ptr_ == NULL )
		{
			throw std::runtime_error("ERROR: Host must be assigned to offset_ring_dens_functor before function can be called.\n");
		}

		BRG_DISTANCE d = lc_add( _R0_, _R_, in_param );

		return _host_ptr_->proj_dens( d, silent );
	}

	offset_ring_dens_functor( const name *new_host=NULL,
			CONST_BRG_DISTANCE_REF new_R_0 = 0, CONST_BRG_DISTANCE_REF new_R = 0 )
	{
		_host_ptr_ = new_host;
		_R_ = new_R;
		_R0_ = new_R_0;
	}

};

template<typename name>
class offset_circ_dens_functor
{
private:
	const name *_host_ptr_;
	BRG_DISTANCE _R0_, _R_;

	BRG_DISTANCE _arc_length_in_circle( CONST_BRG_DISTANCE_REF R2 ) const
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

public:

	void set_host_ptr( const name *new_host_ptr )
	{
		_host_ptr_ = new_host_ptr;
	}
	const name * host_ptr()
	{
		return _host_ptr_;
	}

	void set_R0( CONST_BRG_DISTANCE_REF new_R0 )
	{
		_R0_ = new_R0;
	}
	CONST_BRG_DISTANCE_REF R0()
	{
		return _R0_;
	}

	void set_R( CONST_BRG_DISTANCE_REF new_R )
	{
		_R_ = new_R;
	}
	CONST_BRG_DISTANCE_REF R()
	{
		return _R_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF in_param,
			const bool silent = false ) const
	{
		if ( _host_ptr_ == NULL )
		{
			throw std::runtime_error("ERROR: Host must be assigned to offset_circ_dens_functor before function can be called.\n");
		}

		BRG_DISTANCE L = _arc_length_in_circle(in_param);

		return L * _host_ptr_->proj_dens( in_param, silent );
	}

	offset_circ_dens_functor( const name *new_host=NULL,
			CONST_BRG_DISTANCE_REF new_R0 = 0, CONST_BRG_DISTANCE_REF new_R = 0  )
	{
		_host_ptr_ = new_host;
		_R0_ = new_R0;
		_R_ = new_R;
	}
};

template<typename name>
class offset_Delta_Sigma_functor
{

private:

	const name *_host_ptr_;
	BRG_DISTANCE _R_;

public:

	void set_host_ptr( const name *new_host_ptr )
	{
		_host_ptr_ = new_host_ptr;
	}
	const name * host_ptr()
	{
		return _host_ptr_;
	}

	void set_R( CONST_BRG_DISTANCE_REF new_R )
	{
		_R_ = new_R;
	}
	CONST_BRG_DISTANCE_REF  R()
	{
		return _R_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param,
			const bool silent = false ) const
	{
		if ( _host_ptr_ == NULL )
		{
			throw std::runtime_error("ERROR: Host must be assigned to offset_Delta_Sigma_functor before function can be called.\n");
		}

		return _host_ptr_->offset_Delta_Sigma( _R_, in_param, silent );
	}

	offset_Delta_Sigma_functor( const name *init_host=NULL,
			CONST_BRG_DISTANCE_REF init_R = 0 )
	{
		_host_ptr_ = init_host;
		_R_ = init_R;
	}

};

template<typename name>
class quick_offset_Delta_Sigma_functor
{

private:

	const name *_host_ptr_;
	BRG_DISTANCE _R_;

public:

	void set_host_ptr( const name *new_host_ptr )
	{
		_host_ptr_ = new_host_ptr;
	}
	const name * host_ptr()
	{
		return _host_ptr_;
	}

	void set_R( CONST_BRG_DISTANCE_REF new_R )
	{
		_R_ = new_R;
	}
	CONST_BRG_DISTANCE_REF  R()
	{
		return _R_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param,
			const bool silent = false ) const
	{
		if ( _host_ptr_ == NULL )
		{
			throw std::runtime_error("ERROR: Host must be assigned to offset_Delta_Sigma_functor before function can be called.\n");
		}

		return _host_ptr_->quick_offset_Delta_Sigma( _R_, in_param, silent );
	}

	quick_offset_Delta_Sigma_functor( const name *init_host=NULL,
			CONST_BRG_DISTANCE_REF init_R = 0 )
	{
		_host_ptr_ = init_host;
		_R_ = init_R;
	}

};

template<typename name>
class group_Delta_Sigma_weight_functor
{

private:

	const name *_host_ptr_;
	double _c_;

public:

	void set_host_ptr( const name *new_host_ptr )
	{
		_host_ptr_ = new_host_ptr;
	}
	const name * host_ptr()
	{
		return _host_ptr_;
	}

	void set_c( const double new_c )
	{
		_c_ = new_c;
	}
	const double c()
	{
		return _c_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param,
			const bool silent = false ) const
	{
		if ( _host_ptr_ == NULL )
		{
			throw std::runtime_error("ERROR: Host must be assigned to offset_Delta_Sigma_functor before function can be called.\n");
		}

		if ( _c_ == 0 )
		{
			return 2 * pi * in_param * _host_ptr_->proj_dens( in_param );
		}
		else
		{
			BRG_UNIQUE_PTR<name> group_profile(_host_ptr_->lensing_profile_extension_clone());
			group_profile->set_c(_c_);
			group_profile->set_tau(group_profile->tau()*_c_/_host_ptr_->c());
			return 2 * pi * in_param * group_profile->proj_dens(in_param);
		}
	}

	group_Delta_Sigma_weight_functor( const name *init_host=NULL,
			const double init_c = -1 )
	{
		_host_ptr_ = init_host;
		_c_ = init_c;
	}

};

template<typename name>
class shifted_Delta_Sigma_weight_functor
{

private:

	BRG_DISTANCE _sigma_;

public:

	void set_sigma( CONST_BRG_DISTANCE_REF new_sigma )
	{
		_sigma_ = new_sigma;
	}
	CONST_BRG_DISTANCE_REF sigma()
	{
		return _sigma_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param,
			const bool silent = false ) const
	{
		// The output here is the height of a Rayleigh distribution at in_param
		return in_param/square(_sigma_) * std::exp(-square(in_param)/(2*square(_sigma_)));
	}

	shifted_Delta_Sigma_weight_functor( CONST_BRG_DISTANCE_REF new_sigma=0 )
	{
		_sigma_ = new_sigma;
	}

};

template<typename name>
class shifted_Delta_Sigma_circ_functor
{

private:

	const name *_host_ptr_;
	BRG_DISTANCE _R_, _R_shift_;

public:

	void set_host_ptr( const name *new_host_ptr )
	{
		_host_ptr_ = new_host_ptr;
	}
	const name * host_ptr()
	{
		return _host_ptr_;
	}

	void set_R_shift( CONST_BRG_DISTANCE_REF new_R_shift )
	{
		_R_shift_ = new_R_shift;
	}
	CONST_BRG_DISTANCE_REF R_shift()
	{
		return _R_shift_;
	}

	void set_R( CONST_BRG_DISTANCE_REF new_R )
	{
		_R_ = new_R;
	}
	CONST_BRG_DISTANCE_REF R()
	{
		return _R_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param, const bool silent = false ) const
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
				_host_ptr_->shift_factor(1*brgastro::unitconv::kpctom);

		if(isbad(theta))
		{
			return _host_ptr_->Delta_Sigma(R_actual,silent) +
				extra_shear_factor*sigma_crit(_host_ptr_->z(),2*_host_ptr_->z());
		}
		else
		{
			return _host_ptr_->Delta_Sigma(R_actual,silent)*angle_factor +
					extra_shear_factor*sigma_crit(_host_ptr_->z(),2*_host_ptr_->z()); //TODO Should there be a factor of 1/2 here?
		}
	}

	shifted_Delta_Sigma_circ_functor( const name *new_host=NULL,
			CONST_BRG_DISTANCE_REF new_R_shift=0, CONST_BRG_DISTANCE_REF new_R=0 )
	{
		_host_ptr_ = new_host;
		_R_shift_ = new_R_shift;
		_R_ = new_R;
	}

};

template<typename name>
class shifted_Delta_Sigma_functor
{

private:

	const name *_host_ptr_;
	BRG_DISTANCE _R_;

public:

	void set_host_ptr( const name *new_host_ptr )
	{
		_host_ptr_ = new_host_ptr;
	}
	const name * host_ptr()
	{
		return _host_ptr_;
	}

	void set_R( CONST_BRG_DISTANCE_REF new_R )
	{
		_R_ = new_R;
	}
	CONST_BRG_DISTANCE_REF R()
	{
		return _R_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF in_param, const bool silent = false ) const
	{
		// in_param here will be R_shift
		shifted_Delta_Sigma_circ_functor<name> func(_host_ptr_,in_param,_R_);

		const double min_in_param = 0;
		const double max_in_param = pi;
		const double precision = 0.000001;
		BRG_UNITS out_param(0);

		out_param = brgastro::integrate_Romberg( &func, min_in_param, max_in_param,
				precision, false );

		return out_param /= pi;
	}

	shifted_Delta_Sigma_functor( const name *new_host=NULL,
			CONST_BRG_DISTANCE_REF new_R=0 )
	{
		_host_ptr_ = new_host;
		_R_ = new_R;
	}

};

} // namespace brgastro

#endif /* _BRG_name_FUNCTORS_H_ */
