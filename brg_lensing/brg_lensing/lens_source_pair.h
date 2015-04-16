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

#ifndef _BRG_USE_CPP_11_STD_
#include <boost/shared_ptr.hpp>
#else
#include <memory>
#endif

#include "brg/global.h"

#include "brg_lensing/source_obj.hpp"

#include "brg_physics/sky_obj/sky_obj.h"
#include "brg_physics/units/unit_obj.h"

namespace brgastro {

/**
 *
 */
class lens_source_pair {
private:

	bool _using_lens_clone_;
	bool _using_source_clone_;
	const sky_obj *_init_lens_ptr_;
	const source_obj *_init_source_ptr_;
	BRG_SHARED_PTR<const sky_obj> _lens_clone_;
	BRG_SHARED_PTR<const source_obj> _source_clone_;

	mutable bool _data_stored_;
	mutable double _z_lens_, _z_source_;
	mutable BRG_MASS _m_lens_;
	mutable size_t _id_lens_;
	mutable double _mag_lens_, _mag_source_;
	mutable double _weight_lens_, _weight_source_;
	double _weight_pair_;
	mutable BRG_DISTANCE _R_proj_;
	mutable BRG_ANGLE _theta_;
	mutable double _gamma_t_, _gamma_x_;

	void _conditional_store_data() const
	{
		if(!_data_stored_) store_data();
	}

	void _calc_gamma(const double gamma_1, const double gamma_2) const;

public:

	// Calculation/recalculation function
#if(1)
	void store_data() const;
#endif

	// Constructors and destructor
#if(1)
	lens_source_pair();
	lens_source_pair( const sky_obj* lens_ptr, const source_obj* source_ptr,
			double init_weight_pair=1, bool make_clones=false);
	virtual ~lens_source_pair();
#endif

	// Set lens and source
#if(1)
	void set_lens( const sky_obj *lens_ptr, const bool make_clone=false);
	void set_source( const source_obj *source_ptr, const bool make_clone=false);
#endif

	// Set pair weight
#if(1)
	void set_weight_pair( double new_weight_pair );
#endif

	// Lens and source access
#if(1)
	const sky_obj *lens() const;
	const source_obj *source() const;
#endif // Lens and source access

	// Access to stored values
#if(1)

	const double & z_lens() const;
	const double & z_source() const;
	double z_diff() const;
	const BRG_MASS & m_lens() const;
	const size_t & id_lens() const;
	const double & mag_lens() const;
	const double & mag_source() const;
	const BRG_DISTANCE & R_proj() const;
	const BRG_ANGLE & theta() const;
	const double & gamma_t() const;
	const double & gamma_x() const;
	const double & weight_lens() const;
	const double & weight_source() const;
	const double & weight_pair() const;
	double shear_weight() const;
	double mag_weight() const;

#endif // Access to stored values

	// Calculated values
#if(1)

	BRG_UNITS sigma_crit() const;
	BRG_UNITS delta_Sigma_t() const;
	BRG_UNITS delta_Sigma_x() const;
#endif
};

} // end namespace brgastro

#endif // _BRG_LENS_SOURCE_PAIR_H_INCLUDED_