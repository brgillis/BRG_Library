/**********************************************************************\
 @file expected_count_cache.h
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

// body file: expected_count_cache.cpp

#ifndef _BRG_EXPECTED_COUNT_CACHE_H_INCLUDED_
#define _BRG_EXPECTED_COUNT_CACHE_H_INCLUDED_

#include <string>

#include "brg/global.h"

#include "brg/math/cache/cache_2d.hpp"
#include "brg/physics/units/unit_obj.h"

namespace brgastro {

/**
 *
 */
class expected_count_cache: public brg_cache_2d<expected_count_cache> {
	DECLARE_BRG_CACHE_2D_STATIC_VARS();

	friend class brg_cache_2d<expected_count_cache>;

protected:

	const std::string _name_base() const throw()
	{
		char name_base[BRG_CACHE_ND_NAME_SIZE] = "magexcnt";
		return name_base;
	}

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0);
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function.
	const double _calculate( const double in_param_1, const double in_param_2 ) const;
	void _load_cache_dependencies() const;

public:
	expected_count_cache()
	{
	}
	virtual ~expected_count_cache()
	{
	}
};

} // end namespace brgastro

#endif // _BRG_EXPECTED_COUNT_CACHE_H_INCLUDED_
