/**********************************************************************\
 @file shifting_loader.h
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

// body file: shifting_loader.cpp
#ifndef _BRG_SHIFTING_LOADER_H_INCLUDED_
#define _BRG_SHIFTING_LOADER_H_INCLUDED_

#include <cassert>
#include <string>
#include <vector>

#include "brg/global.h"


namespace brgastro {

/**
 *
 */
class shifting_loader {
	static bool _loaded_;
	static std::string _full_file_name_;
	static std::vector< std::vector<double> > _data_;

	const size_t _zvals_size_;
	std::vector<double> _zvals_;

	const char * _class_source_name_;
	const char * _corr_alpha_name_;

	void _get_file_name();
	void _load();
	size_t _lower_theta_index(double theta);
	size_t _lower_z_index(double z);
public:

	double get(double theta, double z);
	void unload();

	shifting_loader()
	: _zvals_size_(4), _class_source_name_("shifting_loader.cpp"), _corr_alpha_name_("corr_alpha")
	{
		_zvals_.push_back(0.3);
		_zvals_.push_back(0.5);
		_zvals_.push_back(1.0);
		_zvals_.push_back(1.5);

		assert(_zvals_.size()==_zvals_size_);

	}
};

} // end namespace brgastro

#endif // _BRG_SHIFTING_LOADER_H_INCLUDED_
