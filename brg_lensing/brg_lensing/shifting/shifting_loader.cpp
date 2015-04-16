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

#include "brg/utility.hpp"
#include "brg/vector/manipulations.hpp"

#include "brg_lensing/shifting/corr_alpha.hpp"

#include "shifting_loader.h"

#include <brg/file_access/ascii_table.hpp>

// Initialisation of static vars
#if (1)
bool brgastro::shifting_loader::_loaded_(false);
std::vector< std::vector<double> > brgastro::shifting_loader::_data_;
#endif

void brgastro::shifting_loader::_load()
{
	#ifdef _OPENMP
	#pragma omp critical(brg_load_shifting_loader)
	#endif
	{
		if(!_loaded_)
		{
			std::stringstream ss(corr_alph_data);

			_data_ = load_table<double>(ss);

			_data_ = reverse_vertical(_data_);

			assert(ssize(_data_)==_zvals_size_+1);
			assert(ssize(_data_[0])>=2);

			_loaded_ = true;
		}
	}

}
ssize_t brgastro::shifting_loader::_lower_theta_index(double theta)
{
	if(!_loaded_) _load();

	auto size = ssize(_data_[0]);

	for(ssize_t i=1; i<size; ++i)
	{
		if(theta<_data_[0][i])
			return i-1;
	}
	return size-1;
}
ssize_t brgastro::shifting_loader::_lower_z_index(double z)
{
	assert(_zvals_size_>=2);

	for(ssize_t i=1; i<_zvals_size_; ++i)
	{
		if(z<_zvals_[i])
			return i-1;
	}
	return _zvals_size_-1;
}

double brgastro::shifting_loader::get(double t, double z)
{
	if(!_loaded_) _load();

	const ssize_t ti = _lower_theta_index(t);
	const ssize_t zi = _lower_z_index(z);

	const double tlo = _data_[0][ti];
	const double thi = _data_[0][ti+1];
	const double zlo = _zvals_[zi];
	const double zhi = _zvals_[zi+1];

	const double weight = (thi-tlo)*(zhi-zlo);

	double result = 0;

	result += _data_[zi+1][ti]*(zhi-z)*(thi-t);
	result += _data_[zi+1][ti+1]*(zhi-z)*(t-tlo);
	result += _data_[zi+2][ti]*(z-zlo)*(thi-t);
	result += _data_[zi+2][ti+1]*(z-zlo)*(t-tlo);

	result /= weight;

	return result;
}

void brgastro::shifting_loader::unload()
{
	_data_.clear();
	_loaded_ = false;
}
