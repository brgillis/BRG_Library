/**********************************************************************\
 @file limit_vector.hpp
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

#ifndef _BRG_LIMIT_VECTOR_HPP_INCLUDED_
#define _BRG_LIMIT_VECTOR_HPP_INCLUDED_

#include <cassert>
#include <stdexcept>
#include <vector>

#include "brg/global.h"

#include "brg/math/misc_math.hpp"

namespace brgastro {

template<typename T>
std::vector<T> make_limit_vector(T min, T max, T step)
{
	assert(max>min);
	assert(step>0);

	std::vector<T> result(1,min);

	if(!isinf(step))
	{
		if( isinf(min) || isinf(max) )
			throw std::logic_error("Cannot generate limit vector with finite step and infinite limits.");

		T cur_val = min+step;
		while(cur_val<max)
		{
			result.push_back(cur_val);
			cur_val += step;
		}
	}
	result.push_back(max);

	return result;
}

} // namespace brgastro



#endif // _BRG_LIMIT_VECTOR_HPP_INCLUDED_
