/**********************************************************************\
 @file pair_binner.cpp
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

#include <stdexcept>
#include <vector>

#include "brg/global.h"

#include "brg/physics/units/unit_obj.h"
#include "brg/utility.hpp"
#include "brg/vector/summary_functions.hpp"

#include "pair_binner.h"

namespace brgastro {

	// Private methods
#if(1)

void pair_binner::_check_limits()
{
	// Now check they're all monotonically increasing
	// Note that this function returns false if they're too small as well
	if(!is_monotonically_increasing(_R_bin_limits_))
	{
		_valid_limits_ = false;
		return;
	}
	if(!is_monotonically_increasing(_m_bin_limits_))
	{
		_valid_limits_ = false;
		return;
	}
	if(!is_monotonically_increasing(_z_bin_limits_))
	{
		_valid_limits_ = false;
		return;
	}
	if(!is_monotonically_increasing(_mag_bin_limits_))
	{
		_valid_limits_ = false;
		return;
	}

	_valid_limits_ = true;
}

void pair_binner::_sort() const
{
	if(_sorted_) return;

	if(!_valid_limits_) throw std::logic_error("Pairs can't be sorted without valid bin limits.");

	// Check if the pair_bins vector has been generated and is the proper size
	if(_pair_bins_.empty())
	{
		make_array4d(_pair_bins_,_R_bin_limits_.size()-1,_m_bin_limits_.size()-1,
				_z_bin_limits_.size()-1,_mag_bin_limits_.size()-1);
	}
	else
	{
		if(_pair_bins_.size() != _R_bin_limits_.size()-1)
		{
			_resort();
			return;
		}
		if(_pair_bins_[0].size() != _m_bin_limits_.size()-1)
		{
			_resort();
			return;
		}
		if(_pair_bins_[0][0].size() != _z_bin_limits_.size()-1)
		{
			_resort();
			return;
		}
		if(_pair_bins_[0][0][0].size() != _mag_bin_limits_.size()-1)
		{
			_resort();
			return;
		}
	}

	// Now loop through and sort each pair into the proper bin
	// Note that we don't zero the sorting index here. This is intentional, so
	// we don't have to resort previously sorted pairs when new ones are added.
	for(; _sorting_index_<_pairs_.size(); ++_sorting_index_)
	{
		size_t i;

		BRG_DISTANCE R_proj = _pairs_[_sorting_index_].R_proj();
		if((R_proj < _R_bin_limits_.front()) || (R_proj > _R_bin_limits_.back()))
			continue;
		for(i=0; i<_R_bin_limits_.size(); ++i)
		{
			if(_R_bin_limits_[i]>R_proj)
			{
				--i;
				break;
			}
		}
		size_t R_i = i;

		BRG_MASS m = _pairs_[_sorting_index_].m_lens();
		if((m < _m_bin_limits_.front()) || (m > _m_bin_limits_.back()))
			continue;
		for(i=0; i<_m_bin_limits_.size(); ++i)
		{
			if(_m_bin_limits_[i]>m)
			{
				--i;
				break;
			}
		}
		size_t m_i = i;

		double z = _pairs_[_sorting_index_].z_lens();
		if((z < _z_bin_limits_.front()) || (z > _z_bin_limits_.back()))
			continue;
		for(i=0; i<_z_bin_limits_.size(); ++i)
		{
			if(_z_bin_limits_[i]>z)
			{
				--i;
				break;
			}
		}
		size_t z_i = i;

		double mag = _pairs_[_sorting_index_].mag_lens();
		if((mag < _mag_bin_limits_.front()) || (mag > _mag_bin_limits_.back()))
			continue;
		for(i=0; i<_mag_bin_limits_.size(); ++i)
		{
			if(_mag_bin_limits_[i]>mag)
			{
				--i;
				break;
			}
		}
		size_t mag_i = i;

		_pair_bins_[R_i][m_i][z_i][mag_i].add_pair(_pairs_[_sorting_index_]);

	}
}

void pair_binner::_resort() const
{
	_sorting_index_ = 0;
	_sorted_ = false;
	_pair_bins_.clear();
	_sort();
}

#endif

	// Constructors
#if(1)

	// Set limits by vectors
pair_binner::pair_binner(std::vector< BRG_DISTANCE > R_bin_limits,
				std::vector< BRG_MASS > m_bin_limits,
				std::vector< double > z_bin_limits,
				std::vector< double > mag_bin_limits)
:	_R_bin_limits_(R_bin_limits),
 	_m_bin_limits_(m_bin_limits),
 	_z_bin_limits_(z_bin_limits),
 	_mag_bin_limits_(mag_bin_limits),
 	_valid_limits_(false),
 	_sorting_index_(0),
	_sorted_(false)
{
	_check_limits();
}

	// Set limits by min, max, and step
pair_binner::pair_binner(CONST_BRG_DISTANCE_REF R_min,
				CONST_BRG_DISTANCE_REF R_max,
				CONST_BRG_DISTANCE_REF R_step,
				CONST_BRG_MASS_REF m_min=-std::numeric_limits<double>::infinity(),
				CONST_BRG_MASS_REF m_max=std::numeric_limits<double>::infinity(),
				CONST_BRG_MASS_REF m_step=std::numeric_limits<double>::infinity(),
				double z_min=-std::numeric_limits<double>::infinity(),
				double z_max=std::numeric_limits<double>::infinity(),
				double z_step=std::numeric_limits<double>::infinity(),
				double mag_min=-std::numeric_limits<double>::infinity(),
				double mag_max=std::numeric_limits<double>::infinity(),
				double mag_step=std::numeric_limits<double>::infinity())
:	_R_bin_limits_(make_limit_vector(R_min,R_max,R_step)),
 	_R_bin_limits_(make_limit_vector(m_min,m_max,m_step)),
 	_R_bin_limits_(make_limit_vector(z_min,z_max,z_step)),
 	_R_bin_limits_(make_limit_vector(mag_min,mag_max,mag_step)),
 	_valid_limits_(false),
 	_sorting_index_(0),
	_sorted_(false)
{
	_check_limits();
}

#endif // Constructors

// Adding and clearing data
#if(1)

void pair_binner::add_pair( const lens_source_pair & new_pair)
{
	_pairs_.push_back(new_pair);
	_sorted_ = false;
}
void pair_binner::clear_pairs()
{
	_pair_bins_.clear();
	_pairs_.clear();
}

#endif

} // end namespace brgastro
