/**********************************************************************\
 @file expected_count_loader.h
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

// body file: expected_count_loader.cpp

#ifndef _BRG_PHYSICS_LENSING_MAGNIFICATION_EXPECTED_COUNT_LOADER_H_INCLUDED_
#define _BRG_PHYSICS_LENSING_MAGNIFICATION_EXPECTED_COUNT_LOADER_H_INCLUDED_

#include <string>
#include <vector>

#include "IceBRG_main/common.hpp"

#include "IceBRG_main/units/units.hpp"
#include "IceBRG_main/vector/limit_vector.hpp"

namespace IceBRG {

/**
 *
 */
class expected_count_loader {
	static bool _loaded_;
	static IceBRG::limit_vector<flt_t> _z_limits_;
	static std::vector<IceBRG::limit_vector<flt_t>> _mag_limits_;
	static std::vector<std::vector<custom_unit_type<0,0,0,-2,0>>>_smoothed_count_, _smoothed_count_derivative_;
	static std::string _filename_base_, _filename_tail_;

	static void _load();
	static custom_unit_type<0,0,0,-2,0> _get_interp(const flt_t & mag, const flt_t & z,
			const std::vector<std::vector<custom_unit_type<0,0,0,-2,0>>> & table,
			const custom_unit_type<0,0,0,-2,0> & def=units_cast<custom_unit_type<0,0,0,-2,0>>(0));
public:

	// Setting parameters for where the data is stored
#if(1)
	static void set_z_limits(const std::vector<flt_t> & new_limits_vector);
	static void set_z_limits(std::vector<flt_t> && new_limits_vector);
	static void set_filename_base(const std::string & new_filename_base);
	static void set_filename_base(std::string && new_filename_base);
	static void set_filename_tail(const std::string & new_filename_tail);
	static void set_filename_tail(std::string && new_filename_tail);
#endif

	static custom_unit_type<0,0,0,-2,0> get_count(const flt_t & mag, const flt_t & z);
	static custom_unit_type<0,0,0,-2,0> get_derivative(const flt_t & mag, const flt_t & z);

	static void unload();
};

} // end namespace IceBRG

#endif // _BRG_PHYSICS_LENSING_MAGNIFICATION_EXPECTED_COUNT_LOADER_H_INCLUDED_
