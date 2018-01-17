/**********************************************************************\
  @file gabdt.h

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

// body file: brg/physics/astro/SALTSA/gabdt.cpp

#ifndef _GABDT_H_INCLUDED_
#define _GABDT_H_INCLUDED_

#include <vector>

#include "IceBRG_main/common.hpp"

#include "IceBRG_physics/density_profile/detail/density_profile.hpp"
#include "IceBRG_main/units/units.hpp"

namespace IceBRG {

class gabdt
{
	/************************************************************
	 gabdt
	 -----

	 This class represents a useful physical construct, of how
	 much particles in halos have been disrupted by tidal shocking.

	 See the g_ab object from equation 10 of Taylor and Babul
	 (2001). This represents that multiplied by the timestep.

	 \************************************************************/

private:
	mutable bool _is_cached_;
	const density_profile *_host_ptr_;

	distance_type _x_, _y_, _z_, _r_;
	time_type _dt_;
	mutable std::vector< std::vector< long double > > _dv_;

public:

	// Constructors
	gabdt();
	gabdt( const density_profile *init_host, distance_type const & init_x,
			distance_type const & init_y, distance_type const & init_z,
			time_type const & init_dt );

	// Destructor
	virtual ~gabdt()
	{
	}

	// Full clear function
	void clear();

	// Set functions
	void set( const IceBRG::density_profile *new_host_ptr,
			distance_type const & new_x, distance_type const & new_y,
			distance_type const & new_z, time_type const & new_dt );
	void set_pos( distance_type const & new_x, distance_type const & new_y,
			distance_type const & new_z );
	void set_dt( time_type const & dt );
	void set_host_ptr( const density_profile *new_host_ptr );
	void override_zero();

	// Calculation function
	void calc_dv( const bool silent = false ) const;

	// Get functions
	const density_profile * host() const;
	distance_type x() const;
	distance_type y() const;
	distance_type z() const;
	distance_type r() const;
	std::vector< std::vector< long double > > dv() const; // involves calculation if necessary
	long double dv( const int x_i, const int y_i ) const; // involves calculation if necessary

	// Operator overloading
	flt_t operator*( const gabdt & other_gabdt ) const; // Dot-product(ish) operator

	gabdt & operator+=( const gabdt & other_gabdt ); // Addition
	gabdt operator+( const gabdt & other_gabdt ) const;

	gabdt & operator*=( const double scale_fraction ); // Multiplication by a double
	gabdt operator*( const double scale_fraction ) const;

};

class gabdt_functor
{
	/************************************************************
	 gabdt_functor
	 --------------

	 Child of functor

	 This class provides a functor * for getting the 3-D
	 acceleration within a halo.

	 \************************************************************/
public:

	// Constructors
	gabdt_functor();

	// Destructor
	virtual ~gabdt_functor()
	{
	}

	// Host accessor
	const density_profile *host_ptr;

	// Function method
	std::vector< flt_t > operator()( const std::vector< flt_t > & in_params,
			const bool silent = false ) const;

};

} // end namespace IceBRG

#endif /* _GABDT_H_INCLUDED_ */
