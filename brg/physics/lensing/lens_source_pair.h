/**********************************************************************\
 @file lens_source_pair.h
 ------------------

 A class representing a pair of a lens and source object.

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

// body file: lens_source_pair.cpp

#ifndef _BRG_LENS_SOURCE_PAIR_H_INCLUDED_
#define _BRG_LENS_SOURCE_PAIR_H_INCLUDED_

#include <cassert>

#ifndef _BRG_USE_CPP_11_STD_
#include <boost/shared_ptr.hpp>
#else
#include <memory>
#endif

#include "brg/global.h"

#include "brg/math/misc_math.hpp"
#include "brg/physics/lensing/source_obj.hpp"
#include "brg/physics/sky_obj/sky_obj.h"

namespace brgastro {

/**
 *
 */
class lens_source_pair {
private:

	bool _using_clones_;
	const sky_obj *_init_lens_ptr_;
	const source_obj *_init_source_ptr_;
	BRG_SHARED_PTR<const sky_obj> _lens_clone_;
	BRG_SHARED_PTR<const source_obj> _source_clone_;

	mutable bool _data_stored_;
	mutable double _z_lens_, _z_source_;
	mutable BRG_DISTANCE _R_proj_;
	mutable BRG_ANGLE _theta_;
	mutable double _gamma_t_, _gamma_x_;

	void _conditional_store_data()
	{
		if(!_data_stored_) store_data();
	}

	void _calc_gamma(const double gamma_1, const double gamma_2) const
	{
		double rot_cos = cos(-2*_theta_);
		double rot_sin = sin(-2*_theta_);
		_gamma_t_ = -(rot_cos*gamma_1-rot_sin*gamma_2);
		_gamma_x_ = rot_sin*gamma_1+rot_cos*gamma_2;
	}

public:

	// Calculation/recalculation function
#if(1)
	void store_data() const
	{
		assert(lens()!=NULL);
		assert(source()!=NULL);
		const sky_obj *lens_ptr = lens();
		const source_obj *source_ptr = source();

		_z_lens_ = lens_ptr->z();
		_z_source_ = source_ptr->z();

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
	lens_source_pair()
	:	_using_clones_(false),
		_init_lens_ptr_(NULL),
		_init_source_ptr_(NULL),
		_data_stored_(false),
		_z_lens_(0),
		_z_source_(0),
		_R_proj_(0),
		_theta_(0),
		_gamma_t_(0),
		_gamma_x_(0)
	{
	}
	lens_source_pair( const sky_obj* lens_ptr, const source_obj* source_ptr, bool make_clones=false)
	:	_using_clones_(make_clones),
	 	_init_lens_ptr_(lens_ptr),
		_init_source_ptr_(source_ptr)
	{
		_init_lens_ptr_ = lens_ptr;
		_init_source_ptr_ = source_ptr;
		if(make_clones)
		{
			_lens_clone_ = BRG_SHARED_PTR<const sky_obj>(lens_ptr->sky_obj_clone());
			_source_clone_ = BRG_SHARED_PTR<const source_obj>(source_ptr->source_obj_clone());
			_using_clones_ = true;
		}
		else
		{
			_using_clones_ = false;
		}
		store_data();
	}
	virtual ~lens_source_pair()
	{
	}
#endif

	// Lens and source access
#if(1)
	const sky_obj *lens() const
	{
		if(_using_clones_)
		{
			return _lens_clone_.get();
		}
		else
		{
			return _init_lens_ptr_;
		}
	}
	const source_obj *source() const
	{
		if(_using_clones_)
		{
			return _source_clone_.get();
		}
		else
		{
			return _init_source_ptr_;
		}
	}
#endif
};

} // end namespace brgastro

#endif // _BRG_LENS_SOURCE_PAIR_H_INCLUDED_
