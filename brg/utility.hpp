/**********************************************************************\
  @file utility.hpp

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

#ifndef _BRG_UTILITY_HPP_INCLUDED_
#define _BRG_UTILITY_HPP_INCLUDED_

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "global.h"

namespace brgastro
{

// Generic functions
#if (1)

// Set_zero function - a way for other template functions to "clear" or initialize a value in various ways
inline void set_zero( int & obj )
{
	obj = 0;
}
inline void set_zero( long int & obj )
{
	obj = 0;
}
inline void set_zero( long long int & obj )
{
	obj = 0;
}
inline void set_zero( short int & obj )
{
	obj = 0;
}
inline void set_zero( unsigned int & obj )
{
	obj = 0;
}
inline void set_zero( unsigned long int & obj )
{
	obj = 0;
}
inline void set_zero( unsigned long long int & obj )
{
	obj = 0;
}
inline void set_zero( unsigned short int & obj )
{
	obj = 0;
}
inline void set_zero( double & obj )
{
	obj = 0;
}
inline void set_zero( long double & obj )
{
	obj = 0;
}
inline void set_zero( float & obj )
{
	obj = 0;
}
template< typename T >
inline void set_zero( std::vector< T > & vec )
{
	vec.clear();
}
template< typename T >
inline void set_zero( T * &obj )
{
	obj = NULL;
}
template< typename obj_type >
inline void set_zero( obj_type & obj )
{
	obj = obj_type();
}

// Various "make" functions, to allocate dynamic memory.
// After allocating memory, these functions initialize the new variables using the
// set_zero function (see above).

template< typename obj_type >
inline void make_obj( BRG_UNIQUE_PTR<obj_type> & obj_pointer )
{
	obj_pointer = BRG_UNIQUE_PTR<obj_type>(new obj_type);
	set_zero(*obj_pointer);
}

#endif // Ending functions


} // end namespace brgastro

#endif // _BRG_UTILITY_HPP_INCLUDED_
