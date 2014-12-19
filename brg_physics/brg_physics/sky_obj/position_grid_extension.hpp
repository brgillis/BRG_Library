/**********************************************************************\
 @file position_grid_extension.hpp
 ------------------

 An extension for the sky_obj class, adding grid functions for ra, dec,
 and z.

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

// body file: position_grid_extension.cpp
#ifndef _BRG_POSITION_GRID_EXTENSION_HPP_INCLUDED_
#define _BRG_POSITION_GRID_EXTENSION_HPP_INCLUDED_

#include <cassert>

#include "brg/global.h"

#include "brg_physics/astro.h"
#include "brg_physics/sky_obj/position_grid_cache.h"

// Macro definitions
#if(1)

// SPCP: "Static Polymorphic Const Pointer"
#define SPCP(name) static_cast<const name*>(this)

// SPP: "Static Polymorphic Pointer"
#define SPP(name) static_cast<name*>(this)

#endif // Macro definitions

namespace brgastro {

/**
 *
 */
template<typename name>
class position_grid_extension {
private:

	mutable int _ra_grid_, _dec_grid_, _z_grid_;

	mutable bool _ra_grid_cached_, _dec_grid_cached_, _z_grid_cached_;
	mutable unsigned int _local_ra_grid_change_num_, _local_dec_grid_change_num_,
		_local_z_grid_change_num_;

	// Functions to cache values
#if(1)
	void cache_ra_grid() const
	{
		_ra_grid_ = grid_cache::get_ra_grid( SPCP(name)->ra() );
		_ra_grid_cached_ = true;
		_local_ra_grid_change_num_ = grid_cache::ra_grid_change_num();
	}
	void cache_dec_grid() const
	{
		_dec_grid_ = grid_cache::get_dec_grid( SPCP(name)->dec() );
		_dec_grid_cached_ = true;
		_local_dec_grid_change_num_ = grid_cache::dec_grid_change_num();
	}
	void cache_z_grid() const
	{
		_z_grid_ = grid_cache::get_z_grid( SPCP(name)->z() );
		_z_grid_cached_ = true;
		_local_z_grid_change_num_ = grid_cache::z_grid_change_num();
	}
#endif

public:

	// Access methods
#if (1)

	// Safe versions - will check, then recache if necessary
	const int ra_grid() const
	{
		if ( _ra_grid_cached_ )
		{
			if ( _local_ra_grid_change_num_ == grid_cache::ra_grid_change_num() )
				return _ra_grid_;
		}
		cache_ra_grid();
		return _ra_grid_;
	}
	const int dec_grid() const
	{
		if ( _dec_grid_cached_ )
		{
			if ( _local_dec_grid_change_num_ == grid_cache::dec_grid_change_num() )
				return _dec_grid_;
		}
		cache_dec_grid();
		return _dec_grid_;
	}
	const int z_grid() const
	{
		if ( _z_grid_cached_ )
		{
			if ( _local_z_grid_change_num_ == grid_cache::z_grid_change_num() )
				return _z_grid_;
		}
		cache_z_grid();
		return _z_grid_;
	}

	// Unsafe versions, won't check, won't cache (but will do assert check on debug builds)
	const int quick_ra_grid() const
	{
		assert(_ra_grid_cached_);
		assert(_local_ra_grid_change_num_ == grid_cache::ra_grid_change_num());
		return _ra_grid_;
	}
	const int quick_dec_grid() const
	{
		assert(_dec_grid_cached_);
		assert(_local_dec_grid_change_num_ == grid_cache::dec_grid_change_num());
		return _dec_grid_;
	}
	const int quick_z_grid() const
	{
		assert(_z_grid_cached_);
		assert(_local_z_grid_change_num_ == grid_cache::z_grid_change_num());
		return _z_grid_;
	}
#endif

	position_grid_extension()
	{
		_ra_grid_ = _dec_grid_ = _z_grid_ = 0;
		_ra_grid_cached_ = _dec_grid_cached_ = _z_grid_cached_ = false;
		_local_ra_grid_change_num_ = 0;
		_local_dec_grid_change_num_ = 0;
		_local_z_grid_change_num_ = 0;

	}
	virtual ~position_grid_extension()
	{
	}
};

} // end namespace brgastro

// Undef macros
#undef SPP
#undef SPCP

#endif // _BRG_POSITION_GRID_EXTENSION_HPP_INCLUDED_
