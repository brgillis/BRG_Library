/**********************************************************************\
 @file table_utility.hpp
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


// body file: table_utility.cpp


#ifndef _BRG_TABLE_UTILITY_H_INCLUDED_
#define _BRG_TABLE_UTILITY_H_INCLUDED_

#include <string>
#include <vector>

#include "brg/global.h"

#include "brg/file_access/table_map.hpp"
#include "brg/file_access/table_typedefs.hpp"
#include "brg/vector/make_vector.hpp"

namespace brgastro {

// Merge a header and data table into a map
template<typename T>
table_map_t<T> make_table_map(
		const table_t<T> & data,
		const header_t & header,
		const bool silent=false)
{
	table_map_t<T> result;

	size_t h_size = header.size();
	size_t d_size = data.size();

	const header_t *header_to_use = &header;
	header_t new_header;
	const table_t<T> *data_to_use = &data;
	table_t<T> new_data;

	// First, check if the header and data aren't the same size
	if(h_size<d_size)
	{
		// We'll pad the header with default values
		new_header = header;
		new_header.resize(d_size);
		for(size_t i=h_size; i<d_size; ++i)
		{
			std::stringstream ss("col_");
			ss << i+1;
			new_header[i] = ss.str();
		}
		header_to_use = &new_header;
		h_size = d_size;
	}
	else if(d_size<h_size)
	{
		// We'll pad the data with default values
		new_data = data;
		new_data.resize(h_size);
		if(d_size>0)
		{
			// If we have some data, match the size of it in new columns
			for(size_t i=d_size; i<h_size; ++i)
			{
				make_vector_default(new_data[i],new_data[0].size());
			}
		}
		data_to_use = &new_data;
		d_size = h_size;
	}

	for(size_t i=0; i<h_size; ++i)
	{
		result[(*header_to_use)[i]] = (*data_to_use)[i];
	}
	return result;
}

// Splits a string into a vector of "word" strings on whitespace
std::vector< std::string > split_on_whitespace( const std::string & sentence );
header_t convert_to_header( const std::string & line );

} // namespace brgastro

#endif // _BRG_TABLE_UTILITY_H_INCLUDED_
