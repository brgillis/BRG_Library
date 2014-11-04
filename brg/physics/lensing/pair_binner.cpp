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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/iterator/zip_iterator.hpp>

#include "brg/global.h"

#include "brg/file_access/ascii_table.h"
#include "brg/physics/units/unit_obj.h"
#include "brg/utility.hpp"
#include "brg/vector/limit_vector.hpp"
#include "brg/vector/make_vector.hpp"
#include "brg/vector/summary_functions.hpp"

#include "pair_binner.h"

namespace brgastro {

double lens_id::unmasked_frac(const BRG_DISTANCE & R_proj) const
{
	try
	{
		return unmasked_frac_bin_limits.interpolate_bins(R_proj,unmasked_fracs);
	}
	catch(const std::logic_error &e)
	{
		return 1;
	}
}

	// Private methods
#if(1)

void pair_binner::_sort() const
{
	if(_sorted_) return;

	if(!valid_limits()) throw std::logic_error("Pairs can't be sorted without valid bin limits.");

	// Check if the pair_bins vector has been generated and is the proper size
	if(_pair_bins_.empty())
	{
		auto bin_limit_func = [&] (size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
		{
			return pair_bin(
					R_limits()[R_i],R_limits()[R_i+1],
					m_limits()[m_i],m_limits()[m_i+1],
					z_limits()[z_i],z_limits()[z_i+1],
					mag_limits()[mag_i],mag_limits()[mag_i+1],
					_z_buffer_);
		};
		make_vector_function(_pair_bins_,bin_limit_func,R_limits().size()-1,m_limits().size()-1,
				z_limits().size()-1,mag_limits().size()-1);
	}
	else
	{
		if(_pair_bins_.size() != R_limits().size()-1)
		{
			_resort();
			return;
		}
		if(_pair_bins_[0].size() != m_limits().size()-1)
		{
			_resort();
			return;
		}
		if(_pair_bins_[0][0].size() != z_limits().size()-1)
		{
			_resort();
			return;
		}
		if(_pair_bins_[0][0][0].size() != mag_limits().size()-1)
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
		// Check bounds first
		BRG_DISTANCE R_proj = _pairs_[_sorting_index_].R_proj();
		if((R_proj < R_limits().front()) || (R_proj > R_limits().back()))
			continue;
		BRG_MASS m = _pairs_[_sorting_index_].m_lens();
		if((m < m_limits().front()) || (m > m_limits().back()))
			continue;
		double z = _pairs_[_sorting_index_].z_lens();
		if((z < z_limits().front()) || (z > z_limits().back()))
			continue;
		double mag = _pairs_[_sorting_index_].mag_lens();
		if((mag < mag_limits().front()) || (mag > mag_limits().back()))
			continue;

		// Check lens-source z separation against buffer
		if(_pairs_[_sorting_index_].z_diff()<_z_buffer_) continue;

		// Find specific position by using a zip iterator to go over limits vectors
		// and the pair bins vector at the same time
		auto R_zip_begin = boost::make_zip_iterator(boost::make_tuple(
				R_limits().begin()+1,_pair_bins_.begin()));
		auto R_zip_end = boost::make_zip_iterator(boost::make_tuple(
				R_limits().end(),_pair_bins_.end()));

		auto R_func = &t1first_lt_v2<decltype(*R_zip_begin),BRG_DISTANCE>;

		auto & pair_R_vec = std::lower_bound(R_zip_begin,R_zip_end,R_proj,R_func)->get<1>();



		auto m_zip_begin = boost::make_zip_iterator(boost::make_tuple(
				m_limits().begin()+1,pair_R_vec.begin()));
		auto m_zip_end = boost::make_zip_iterator(boost::make_tuple(
				m_limits().end(),pair_R_vec.end()));

		auto m_func = &t1first_lt_v2<decltype(*m_zip_begin),BRG_MASS>;

		auto & pair_Rm_vec = std::lower_bound(m_zip_begin,m_zip_end,m,m_func)->get<1>();



		auto z_zip_begin = boost::make_zip_iterator(boost::make_tuple(
				z_limits().begin()+1,pair_Rm_vec.begin()));
		auto z_zip_end = boost::make_zip_iterator(boost::make_tuple(
				z_limits().end(),pair_Rm_vec.end()));

		auto z_func = &t1first_lt_v2<decltype(*z_zip_begin),double>;

		auto & pair_Rmz_vec = std::lower_bound(z_zip_begin,z_zip_end,z,z_func)->get<1>();



		auto mag_zip_begin = boost::make_zip_iterator(boost::make_tuple(
				mag_limits().begin()+1,pair_Rmz_vec.begin()));
		auto mag_zip_end = boost::make_zip_iterator(boost::make_tuple(
				mag_limits().end(),pair_Rmz_vec.end()));

		auto mag_func = &t1first_lt_v2<decltype(*mag_zip_begin),double>;

		auto & pair_Rmzmag = std::lower_bound(mag_zip_begin,mag_zip_end,mag,mag_func)->get<1>();

		// At this point, pair_Rmzmag is a reference to the pair bin we want to add this pair to
		pair_Rmzmag.add_pair(_pairs_[_sorting_index_]);
	}

	// And also add lens ids to each relevant bin
	for(const auto & lens : _lens_ids_)
	{
		// Check bounds first
		const BRG_MASS & m = lens.m;
		if((m < m_limits().front()) || (m > m_limits().back()))
			continue;
		const double & z = lens.z;
		if((z < z_limits().front()) || (z > z_limits().back()))
			continue;
		const double & mag = lens.mag;
		if((mag < mag_limits().front()) || (mag > mag_limits().back()))
			continue;

		// We'll add each id to every radial bin
		for(auto & R_bin : _pair_bins_)
		{
			auto m_zip_begin = boost::make_zip_iterator(boost::make_tuple(
					m_limits().begin()+1,R_bin.begin()));
			auto m_zip_end = boost::make_zip_iterator(boost::make_tuple(
					m_limits().end(),R_bin.end()));

			auto m_func = &t1first_lt_v2<decltype(*m_zip_begin),BRG_MASS>;

			auto & Rm_vec = std::lower_bound(m_zip_begin,m_zip_end,m,m_func)->get<1>();



			auto z_zip_begin = boost::make_zip_iterator(boost::make_tuple(
					z_limits().begin()+1,Rm_vec.begin()));
			auto z_zip_end = boost::make_zip_iterator(boost::make_tuple(
					z_limits().end(),Rm_vec.end()));

			auto z_func = &t1first_lt_v2<decltype(*z_zip_begin),double>;

			auto & Rmz_vec = std::lower_bound(z_zip_begin,z_zip_end,z,z_func)->get<1>();



			auto mag_zip_begin = boost::make_zip_iterator(boost::make_tuple(
					mag_limits().begin()+1,Rmz_vec.begin()));
			auto mag_zip_end = boost::make_zip_iterator(boost::make_tuple(
					mag_limits().end(),Rmz_vec.end()));

			auto mag_func = &t1first_lt_v2<decltype(*mag_zip_begin),double>;

			auto & Rmzmag_bin = std::lower_bound(mag_zip_begin,mag_zip_end,mag,mag_func)->get<1>();

			// At this point, pair_Rmzmag is a reference to the pair bin we want to add this pair to
			Rmzmag_bin.add_lens(lens.id,lens.z,lens.unmasked_frac(Rmzmag_bin.R_mid()));
		}
	}

	// Update the summaries
	make_vector_coerce<4>(_pair_bin_summaries_,_pair_bins_);
}

void pair_binner::_resort() const
{
	_sorting_index_ = 0;
	_sorted_ = false;
	_pair_bins_.clear();
	_sort();
}

#endif // Private functions

// Adding and clearing data
#if(1)

bool pair_binner::binnable( const galaxy & lens ) const
{
	if(!valid_limits()) return false;
	auto m = lens.m();
	if(m<m_limits().front()) return false;
	if(m>m_limits().back()) return false;
	auto z = lens.z();
	if(z<z_limits().front()) return false;
	if(z>z_limits().back()) return false;
	auto mag = lens.mag();
	if(mag<mag_limits().front()) return false;
	if(mag>mag_limits().back()) return false;
	return true;
}
void pair_binner::add_pair( const lens_source_pair & new_pair)
{
	_pairs_.push_back(new_pair);
	_sorted_ = false;
}
void pair_binner::add_lens_id( const size_t & new_lens_id, const BRG_MASS & m, const double & z,
		const double & mag)
{
	_lens_ids_.insert(lens_id(new_lens_id,m,z,mag,_unmasked_frac_bin_limits_,_unmasked_fracs_));
	_sorted_ = false;
}
void pair_binner::set_unmasked_fractions( const std::vector<double> & bin_limits,
		const std::vector<double> & unmasked_fractions)
{
	// Sort at this point, so any lenses already added get the right unmasked fraction
	_sort();

	_unmasked_frac_bin_limits_ = bin_limits;
	_unmasked_fracs_ = unmasked_fractions;
}
void pair_binner::clear()
{
	_sorting_index_ = 0;
	_pair_bin_summaries_.clear();
	_pair_bins_.clear();
	_pairs_.clear();
	_lens_ids_.clear();
}
void pair_binner::empty()
{
	_sorting_index_ = 0;
	_pair_bin_summaries_.clear();
	_pairs_.clear();
	_lens_ids_.clear();

	// Empty each bin in turn, but don't clear the structure
	for(auto R_it=_pair_bins_.begin(); R_it!=_pair_bins_.end(); ++R_it)
	{
		for(auto Rm_it=R_it->begin(); Rm_it!=R_it->end(); ++Rm_it)
		{
			for(auto Rmz_it=Rm_it->begin(); Rmz_it!=Rm_it->end(); ++Rmz_it)
			{
				for(auto Rmzmag_it=Rmz_it->begin(); Rmzmag_it!=Rmz_it->end(); ++Rmzmag_it)
				{
					Rmzmag_it->clear();
				}
			}
		}
	}
}
void pair_binner::sort() const
{
	_sort();
}

#endif // Adding and clearing data

// Accessing summary data for bins
#if(1)

// Access by index (will throw if out of bounds)
#if(1)
BRG_UNITS pair_binner::delta_Sigma_t_mean_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	_sort();
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_mean();
}
BRG_UNITS pair_binner::delta_Sigma_x_mean_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	_sort();
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_x_mean();
}

BRG_UNITS pair_binner::delta_Sigma_t_std_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	_sort();
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_std();
}
BRG_UNITS pair_binner::delta_Sigma_x_std_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	_sort();
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_x_std();
}

BRG_UNITS pair_binner::delta_Sigma_t_stderr_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	_sort();
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_stderr();
}
BRG_UNITS pair_binner::delta_Sigma_x_stderr_for_bin(size_t R_i, size_t m_i, size_t z_i, size_t mag_i)
{
	_sort();
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_stderr();
}
#endif // Access by index

// Access by position
#if(1)
BRG_UNITS pair_binner::delta_Sigma_t_mean_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	_sort();
	size_t R_i = R_limits().get_bin_index(R);
	size_t m_i = m_limits().get_bin_index(m);
	size_t z_i = z_limits().get_bin_index(z);
	size_t mag_i = mag_limits().get_bin_index(mag);
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_mean();
}
BRG_UNITS pair_binner::delta_Sigma_x_mean_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	_sort();
	size_t R_i = R_limits().get_bin_index(R);
	size_t m_i = m_limits().get_bin_index(m);
	size_t z_i = z_limits().get_bin_index(z);
	size_t mag_i = mag_limits().get_bin_index(mag);
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_x_mean();
}

BRG_UNITS pair_binner::delta_Sigma_t_std_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	_sort();
	size_t R_i = R_limits().get_bin_index(R);
	size_t m_i = m_limits().get_bin_index(m);
	size_t z_i = z_limits().get_bin_index(z);
	size_t mag_i = mag_limits().get_bin_index(mag);
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_std();
}
BRG_UNITS pair_binner::delta_Sigma_x_std_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	_sort();
	size_t R_i = R_limits().get_bin_index(R);
	size_t m_i = m_limits().get_bin_index(m);
	size_t z_i = z_limits().get_bin_index(z);
	size_t mag_i = mag_limits().get_bin_index(mag);
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_x_std();
}

BRG_UNITS pair_binner::delta_Sigma_t_stderr_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	_sort();
	size_t R_i = R_limits().get_bin_index(R);
	size_t m_i = m_limits().get_bin_index(m);
	size_t z_i = z_limits().get_bin_index(z);
	size_t mag_i = mag_limits().get_bin_index(mag);
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_t_stderr();
}
BRG_UNITS pair_binner::delta_Sigma_x_stderr_for_bin(CONST_BRG_DISTANCE_REF R, CONST_BRG_MASS_REF m,
		double z, double mag)
{
	_sort();
	size_t R_i = R_limits().get_bin_index(R);
	size_t m_i = m_limits().get_bin_index(m);
	size_t z_i = z_limits().get_bin_index(z);
	size_t mag_i = mag_limits().get_bin_index(mag);
	return _pair_bins_.at(R_i).at(m_i).at(z_i).at(mag_i).delta_Sigma_x_stderr();
}
#endif // Access by index

#endif // Accessing summary data for bins

} // end namespace brgastro
