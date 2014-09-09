/**********************************************************************\
 @file source_obj.hpp
 ------------------

 A class representing a source object (for lensing) in the sky. It's an
 extension of sky_obj which adds gamma and kappa.

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

// body file: source_obj.cpp
#ifndef _BRG_SOURCE_OBJ_H_INCLUDED_
#define _BRG_SOURCE_OBJ_H_INCLUDED_

#include <brg/physics/sky_obj/sky_obj.h>

// SPCP: "Static Polymorphic Const Pointer"
#define SPCP(name) static_cast<const name*>(this)

// SPP: "Static Polymorphic Pointer"
#define SPP(name) static_cast<name*>(this)

namespace brgastro {

/**
 *
 */
class source_obj: public sky_obj {
private:

	double _gamma_1_, _gamma_2_; // Gamma components in ra/dec coordinates
	double _kappa_; // Convergence

public:

	// Constructors and destructor
#if(1)
	source_obj( CONST_BRG_ANGLE_REF init_ra = 0, CONST_BRG_ANGLE_REF init_dec = 0, double init_z = 0,
			double init_gamma_1 = 0, double init_gamma_2 = 0, double init_kappa = 0 )
	: sky_obj(init_ra,init_dec,init_z),
	  _gamma_1_(init_gamma_1),
	  _gamma_2_(init_gamma_2),
	  _kappa_(init_kappa)
	{
	}
	virtual ~source_obj()
	{
	}
#endif

	// Clone function
	virtual source_obj *source_obj_clone() const =0;

	// Accessors
#if (1)

	double gamma_1() const
	{
		return _gamma_1_;
	}
	double gamma_2() const
	{
		return _gamma_2_;
	}
	double kappa() const
	{
		return _kappa_;
	}

#endif

	// Setters
#if (1)

	void set_gamma_1(double new_gamma_1)
	{
		_gamma_1_ = new_gamma_1;
	}
	void set_gamma_2(double new_gamma_2)
	{
		_gamma_2_ = new_gamma_2;
	}
	void set_kappa(double new_kappa)
	{
		_kappa_ = new_kappa;
	}

#endif

};

} // end namespace brgastro

// Undef macros
#undef SPP
#undef SPCP

#endif // _BRG_SOURCE_OBJ_H_INCLUDED_
