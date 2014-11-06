/**********************************************************************\
 @file coerce.hpp
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

#ifndef _BRG_CONTAINER_COERCE_HPP_INCLUDED_
#define _BRG_CONTAINER_COERCE_HPP_INCLUDED_

#include <iterator>

namespace brgastro {

// Helper containers
#if(1)

// Assignment coercer - uses assignment to coerce
#if(1)

/**
 * assignment_coercer - a structure which performs coercion of one type to
 * another as a side-effect of its constructor.
 */
template<unsigned short d, typename container, typename other_container>
struct assignment_coercer
{
	assignment_coercer(container & obj, const other_container & other_obj)
	{
		obj.resize(other_obj.size());
		auto o_it = begin(other_obj);

		for(auto it=begin(obj); it!=end(obj); ++it, ++o_it)
		{
			assignment_coercer<d-1,decltype(*it),decltype(*o_it)>(*it,*o_it);
		}
	}
};

/**
 * d=0 partial specialization of assignment_coercer. Uses simple assignment to
 * coerce.
 */
template<typename container, typename other_container>
struct assignment_coercer<0,container,other_container>
{
	assignment_coercer(container & obj, const other_container & other_obj)
	{
		obj = other_obj;
	}
};

#endif

// Range coercer - uses range constructor to coerce
#if(1)

/**
 * range_coercer - a structure which performs coercion of one type to
 * another as a side-effect of its constructor.
 */
template<unsigned short d, typename container, typename other_container>
struct range_coercer
{
	range_coercer(container & obj, const other_container & other_obj)
	{
		obj.resize(other_obj.size());
		auto o_it = begin(other_obj);

		for(auto it=begin(obj); it!=end(obj); ++it, ++o_it)
		{
			assignment_coercer<d-1,decltype(*it),decltype(*o_it)>(*it,*o_it);
		}
	}
};

/**
 * d=1 partial specialization of range_coercer. Uses range constructor to
 * coerce.
 */
template<typename container, typename other_container>
struct range_coercer<1,container,other_container>
{
	range_coercer(container & obj, const other_container & other_obj)
	{
		obj = container(begin(other_obj),end(other_obj));
	}
};

/**
 * d=0 partial specialization of range_coercer. Uses simple assignment to
 * coerce.
 */
template<typename container, typename other_container>
struct range_coercer<0,container,other_container>
{
	range_coercer(container & obj, const other_container & other_obj)
	{
		obj = other_obj;
	}
};

#endif

#endif

/**
 * coerce - helper function which returns a coercion of old into NewType.
 *
 * @param old
 * @return
 */
template<typename NewType, short unsigned d=1, typename OldType=NewType>
NewType coerce(const OldType & old)
{
	NewType result;
	assignment_coercer<d,NewType,OldType>(result,old);
	return result;
}

/**
 * range_coerce - helper function which returns a coercion of old into NewType.
 * Use this version if the lowest level container has a range constructor for
 * a moderate performance gain.
 *
 * @param old
 * @return
 */
template<typename NewType, short unsigned d=1, typename OldType=NewType>
NewType range_coerce(const OldType & old)
{
	NewType result;
	range_coercer<d,NewType,OldType>(result,old);
	return result;
}

} // namespace brgastro


#endif // _BRG_CONTAINER_COERCE_HPP_INCLUDED_
