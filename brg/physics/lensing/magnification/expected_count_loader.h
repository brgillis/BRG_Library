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

#include "brg/global.h"

namespace brgastro {

/**
 *
 */
class expected_count_loader {
	static bool _loaded_;
	static std::vector<double> _z_limits_;
	static std::vector<std::vector<double>> _mag_limits_, _smoothed_count_, _smoothed_count_derivative_;
	static std::string _filename_base_, _filename_tail_;

	static void _load();
	static double _get_interp(const double & mag, const double & z,
			const std::vector<std::vector<double>> & table,
			const double & def=0);
public:

	// Setting parameters for where the data is stored
#if(1)
	static void set_z_limits(const std::vector<double> & new_limits_vector);
	static void set_z_limits(std::vector<double> && new_limits_vector);
	static void set_filename_base(const std::string & new_filename_base);
	static void set_filename_base(std::string && new_filename_base);
	static void set_filename_tail(const std::string & new_filename_tail);
	static void set_filename_tail(std::string && new_filename_tail);
#endif

	static double get_count(const double & mag, const double & z);
	static double get_derivative(const double & mag, const double & z);

	static void unload();
};

} // end namespace brgastro

#endif // _BRG_PHYSICS_LENSING_MAGNIFICATION_EXPECTED_COUNT_LOADER_H_INCLUDED_
