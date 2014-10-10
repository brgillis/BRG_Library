/**********************************************************************\
 @file shifting_loader.cpp
 ------------------

 This class is used to load in the data tables provided by Chang and
 Jain for the mean relative shift of two points in space, so the data
 can be stored in the same format as my other caches.

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

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "brg/global.h"

#include "brg/file_access/ascii_table_map.hpp"
#include "brg/physics/lensing/magnification/count_fitting_results.hpp"
#include "brg/vector/elementwise_functions.hpp"

#include "expected_count_loader.h"

// Initialisation of static vars
#if (1)
bool brgastro::expected_count_loader::_loaded_(false);
brgastro::table_map_t<double> brgastro::expected_count_loader::_data_map_;
const char * brgastro::expected_count_loader::_class_source_name_ = "expected_count_loader.cpp";
const char * brgastro::expected_count_loader::_data_file_name_ = "count_fitting_results";
#endif

const double num_columns=8;

void brgastro::expected_count_loader::_load()
{
#pragma omp critical(brg_load_expected_count_loader)
	{
		if(!_loaded_)
		{
			std::stringstream ss(count_fitting_result_data);

			_data_map_ = load_table_map<double>(ss);

			assert(_data_map_.size()==num_columns);
			assert(_data_map_.at("z_mid").size()>=2);

			_loaded_ = true;
		}
	}

}
size_t brgastro::expected_count_loader::_lower_z_index(double z)
{
	assert(_data_map_.at("z_mid").size()>=2);

	for(size_t i=1; i<_data_map_.at("z_mid").size(); ++i)
	{
		if(z<_data_map_["z_mid"][i])
			return i-1;
	}
	return _data_map_["z_mid"].size()-2;
}

std::vector<long double> brgastro::expected_count_loader::get(long double z)
{
	if(!_loaded_) _load();

	const size_t zi = _lower_z_index(z);

	const long double zlo = _data_map_["z_mid"][zi];
	const long double zhi = _data_map_["z_mid"][zi+1];

	const long double weight = zhi-zlo;

	std::vector<long double> r_lo, r_hi;

#ifdef _BRG_USE_UNITS_
	r_lo.push_back(unit_obj(_data_map_["N_scale"].at(zi)*(zhi-z),0,0,0,0,-2));
#else
	r_lo.push_back(_data_map_["N_scale"].at(zi)*(zhi-z));
#endif
	r_lo.push_back(_data_map_["m_star_lower"].at(zi)*(zhi-z));
	r_lo.push_back(_data_map_["alpha"].at(zi)*(zhi-z));
	r_lo.push_back(_data_map_["lower_cutoff_sharpness"].at(zi)*(zhi-z));
#ifdef _BRG_USE_UNITS_
	r_lo.push_back(unit_obj(_data_map_["mag23_jump"].at(zi)*(zhi-z),0,0,0,0,-2));
#else
	r_lo.push_back(_data_map_["mag23_jump"].at(zi)*(zhi-z));
#endif
	r_lo.push_back(_data_map_["m_star_upper"].at(zi)*(zhi-z));
	r_lo.push_back(_data_map_["upper_cutoff_sharpness"].at(zi)*(zhi-z));

#ifdef _BRG_USE_UNITS_
	r_hi.push_back(unit_obj(_data_map_["N_scale"].at(zi+1)*(z-zlo),0,0,0,0,-2));
#else
	r_hi.push_back(_data_map_["N_scale"].at(zi+1)*(z-zlo));
#endif
	r_hi.push_back(_data_map_["m_star_lower"].at(zi+1)*(z-zlo));
	r_hi.push_back(_data_map_["alpha"].at(zi+1)*(z-zlo));
	r_hi.push_back(_data_map_["lower_cutoff_sharpness"].at(zi+1)*(z-zlo));
#ifdef _BRG_USE_UNITS_
	r_hi.push_back(unit_obj(_data_map_["mag23_jump"].at(zi+1)*(z-zlo),0,0,0,0,-2));
#else
	r_hi.push_back(_data_map_["mag23_jump"].at(zi+1)*(z-zlo));
#endif
	r_hi.push_back(_data_map_["m_star_upper"].at(zi+1)*(z-zlo));
	r_hi.push_back(_data_map_["upper_cutoff_sharpness"].at(zi+1)*(z-zlo));

	return divide(add(r_lo,r_hi),weight);
}

void brgastro::expected_count_loader::unload()
{
	_data_map_.clear();
	_loaded_ = false;
}
