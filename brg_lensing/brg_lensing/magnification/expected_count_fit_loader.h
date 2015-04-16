/**********************************************************************\
 @file expected_count_loader.h
 ------------------

 This class is used to load in the data tables for fitting results to
 expected count in various magnitude and redshift bins of the CFHTLenS.

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
#ifndef _BRG_EXPECTED_COUNT_FIT_LOADER_H_INCLUDED_
#define _BRG_EXPECTED_COUNT_FIT_LOADER_H_INCLUDED_

#include <cassert>
#include <string>
#include <vector>

#include "brg/global.h"

#include "brg/container/table_typedefs.hpp"

namespace brgastro {

/**
 *
 */
class expected_count_fit_loader {
	static bool _loaded_;
	static brgastro::table_map_t<double> _data_map_;

	static void _load();
	static ssize_t _lower_z_index(double z);
public:

	static std::vector<long double> get(long double z);
	static void unload();
};

} // end namespace brgastro

#endif // _BRG_EXPECTED_COUNT_FIT_LOADER_H_INCLUDED_
