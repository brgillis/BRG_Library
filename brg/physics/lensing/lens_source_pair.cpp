/**********************************************************************\
 @file lens_source_pair.cpp
 ------------------

 TODO <Insert file description here>

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


#include <cassert>

#ifndef _BRG_USE_CPP_11_STD_
#include <boost/shared_ptr.hpp>
#else
#include <memory>
#endif

#include "brg/global.h"

#include "brg/physics/lensing/source_obj.hpp"
#include "brg/physics/sky_obj/sky_obj.h"
#include "brg/math/misc_math.hpp"
#include "brg/physics/units/unit_obj.h"

#include "lens_source_pair.h"

namespace brgastro {

#if(1) // Private methods

void lens_source_pair::_calc_gamma(const double gamma_1, const double gamma_2) const
{
	double rot_cos = cos(-2*_theta_);
	double rot_sin = sin(-2*_theta_);
	_gamma_t_ = -(rot_cos*gamma_1-rot_sin*gamma_2);
	_gamma_x_ = rot_sin*gamma_1+rot_cos*gamma_2;
}

#endif // Private methods

#if(1) // Public methods

// Calculation/recalculation function
#if(1)
void lens_source_pair::store_data() const
{
	assert(lens()!=NULL);
	assert(source()!=NULL);
	const sky_obj *lens_ptr = lens();
	const source_obj *source_ptr = source();

	_z_lens_ = lens_ptr->z();
	_z_source_ = source_ptr->z();
	_m_lens_ = lens_ptr->m();
	_mag_lens_ = lens_ptr->mag();
	_mag_source_ = source_ptr->mag();

	// Note minus sign in dra to correct for ra's reverse orientation in the sky
	BRG_ANGLE dra = -(source_ptr->ra()-lens_ptr->ra())*cos(lens_ptr->dec());
	BRG_ANGLE ddec = source_ptr->dec()-lens_ptr->dec();

	_R_proj_ = dfa(dist2d(dra,ddec),_z_lens_);
	_theta_ = atan2(ddec,dra);

	_calc_gamma(source_ptr->gamma_1(),source_ptr->gamma_2());
}
#endif

// Constructors and destructor
#if(1)
lens_source_pair::lens_source_pair()
:	_using_lens_clone_(false),
 	_using_source_clone_(false),
	_init_lens_ptr_(NULL),
	_init_source_ptr_(NULL),
	_data_stored_(false),
	_z_lens_(0),
	_z_source_(0),
	_m_lens_(0),
	_mag_lens_(0),
	_mag_source_(0),
	_R_proj_(0),
	_theta_(0),
	_gamma_t_(0),
	_gamma_x_(0)
{
}
lens_source_pair::lens_source_pair( const sky_obj* lens_ptr, const source_obj* source_ptr,
		bool make_clones)
:	_using_lens_clone_(make_clones),
 	_using_source_clone_(make_clones),
 	_init_lens_ptr_(lens_ptr),
	_init_source_ptr_(source_ptr),
	_data_stored_(false),
	_z_lens_(0),
	_z_source_(0),
	_m_lens_(0),
	_mag_lens_(0),
	_mag_source_(0),
	_R_proj_(0),
	_theta_(0),
	_gamma_t_(0),
	_gamma_x_(0)
{
	_init_lens_ptr_ = lens_ptr;
	_init_source_ptr_ = source_ptr;
	if(make_clones)
	{
		_lens_clone_ = BRG_SHARED_PTR<const sky_obj>(lens_ptr->sky_obj_clone());
		_source_clone_ = BRG_SHARED_PTR<const source_obj>(source_ptr->source_obj_clone());
	}
	store_data();
}
lens_source_pair::~lens_source_pair()
{
}
#endif

// Set lens and source
#if(1)

void lens_source_pair::set_lens( const sky_obj *lens_ptr, const bool make_clone)
{
	assert(lens_ptr!=NULL);

	if(make_clone)
	{
		_lens_clone_ = BRG_SHARED_PTR<const sky_obj>(lens_ptr->sky_obj_clone());
	}
	else
	{
		_init_lens_ptr_ = lens_ptr;
	}
	_using_lens_clone_ = make_clone;
}

void lens_source_pair::set_source( const source_obj *source_ptr, const bool make_clone)
{
	assert(source_ptr!=NULL);

	if(make_clone)
	{
		_source_clone_ = BRG_SHARED_PTR<const source_obj>(source_ptr->source_obj_clone());
	}
	else
	{
		_init_source_ptr_ = source_ptr;
	}
	_using_source_clone_ = make_clone;
}

#endif


// Lens and source access
#if(1)
const sky_obj *lens_source_pair::lens() const
{
	if(_using_lens_clone_)
	{
		return _lens_clone_.get();
	}
	else
	{
		return _init_lens_ptr_;
	}
}
const source_obj *lens_source_pair::source() const
{
	if(_using_source_clone_)
	{
		return _source_clone_.get();
	}
	else
	{
		return _init_source_ptr_;
	}
}
#endif // Lens and source access

// Access to stored values
#if(1)

double lens_source_pair::z_lens() const
{
	_conditional_store_data();
	return _z_lens_;
}
double lens_source_pair::z_source() const
{
	_conditional_store_data();
	return _z_source_;
}
BRG_DISTANCE lens_source_pair::R_proj()
{
	_conditional_store_data();
	return _R_proj_;
}
BRG_ANGLE lens_source_pair::theta()
{
	_conditional_store_data();
	return _theta_;
}
double lens_source_pair::gamma_t()
{
	_conditional_store_data();
	return _gamma_t_;
}
double lens_source_pair::gamma_x()
{
	_conditional_store_data();
	return _gamma_x_;
}

#endif // Access to stored values

// Calculated values
#if(1)

BRG_UNITS lens_source_pair::sigma_crit()
{
	_conditional_store_data();
	return square( c ) / ( 4. * pi * Gc )
			* ad_distance( 0, _z_source_ )
			/ ( ad_distance( 0, _z_lens_ ) * ad_distance( _z_lens_, _z_source_ ) );
}

BRG_UNITS lens_source_pair::delta_Sigma_t()
{
	return sigma_crit()*gamma_t();
}

BRG_UNITS lens_source_pair::delta_Sigma_x()
{
	return sigma_crit()*gamma_x();
}
#endif // Calculated values

#endif // Public methods

} // end namespace brgastro
