/**********************************************************************\
 @file manipulations.hpp
 ------------------

 Functions to manipulate vectors in various ways.

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

#ifndef _BRG_MANIPULATIONS_HPP_INCLUDED_
#define _BRG_MANIPULATIONS_HPP_INCLUDED_

#include <vector>

#include "brg/global.h"

#include "brg/vector/make_vector.hpp"

namespace brgastro {

template<typename T>
std::vector< std::vector<T> > pad(const std::vector< std::vector<T> > & v, const T & default_value=T())
{
	size_t n_cols = v.size();
	size_t n_rows = 0;
	for(size_t i=0; i<n_cols; ++i)
		if(v[i].size() > n_rows) n_rows = v[i].size();

	std::vector< std::vector<T> > result = v;

	for(size_t i=0; i<n_cols; ++i)
		if(result[i].size() < n_rows) result[i].resize(n_rows,default_value);

	return result;
}

template<typename T>
std::vector< std::vector<T> > transpose(const std::vector< std::vector<T> > & v, const T & default_value=T())
{
	size_t n_cols = v.size();
	size_t n_rows = 0;
	for(size_t i=0; i<v.size(); ++i)
		if(v[i].size() > n_rows) n_rows = v[i].size();

	std::vector< std::vector<T> > result;

	auto result_function = [&] (size_t i, size_t j)
		{
			return v[j][i];
		};

	make_vector_function(result,result_function,n_rows,n_cols);

	return result;
}

template<typename T>
std::vector< std::vector<T> > reverse_vertical(const std::vector< std::vector<T> > & v)
{
	size_t n_cols = v.size();
	size_t n_rows = 0;
	for(size_t i=0; i<v.size(); ++i)
		if(v[i].size() > n_rows) n_rows = v[i].size();

	std::vector< std::vector<T> > result;

	auto result_function = [&] (size_t i, size_t j)
		{
			return v[i][n_rows-j-1];
		};

	make_vector_function(result,result_function,n_cols,n_rows);

	return result;
}

} // namespace brgastro

#endif // _BRG_MANIPULATIONS_HPP_INCLUDED_
