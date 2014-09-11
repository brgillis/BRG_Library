/**********************************************************************\
  @file file_functions.h

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

// body file: brg/file_access/file_functions.cpp

#ifndef _BRG_FILE_FUNCTIONS_H_INCLUDED_
#define _BRG_FILE_FUNCTIONS_H_INCLUDED_

#include <fstream>
#include <iostream>

#include "brg/global.h"

#include "brg/file_access/table_typedefs.hpp"
#include "brg/utility.hpp"

namespace brgastro
{

/** Global function declarations **/
#if (1)

// Functions to open a file and check that it's been opened successfully. An exception will be thrown
// if the file isn't opened successfully.
void open_file( std::ofstream & stream, const std::string & name,
		const bool silent = false );
void open_file( std::ifstream & stream, const std::string & name,
		const bool silent = false );
void open_file( std::fstream & stream, const std::string & name,
		const bool silent = false );

// Functions to get rid of comments lines (those starting with #) from open fstreams and ifstreams. The "one_line" versions will
// trim the next line if it's commented out, and the "all_at_top" versions will trim all commented lines from the current
// position until they find a line which isn't a comment.
// The "one_line" functions will return 1 if the file is already at the end, and 0 otherwise.
// The "all_at_top" functions will return 1 if they reach the end of the file before the run out of comments, and 0 otherwise.
void trim_comments_one_line( std::istream & stream, const bool silent =
		false );
void trim_comments_one_line( std::fstream & stream, const bool silent =
		false );
void trim_comments_all_at_top( std::istream & stream, const bool silent =
		false );
void trim_comments_all_at_top( std::fstream & stream, const bool silent =
		false );

// Utility functions

// Merge a header and data table into a map
template<typename T>
typename table_map<T>::type make_table_map(
		const typename table<T>::type & data,
		const header::type & header,
		const bool silent=false)
{
	typename table_map<T>::type result;

	size_t h_size = header.size();
	size_t d_size = data.size();

	const header::type *header_to_use = &header;
	header::type new_header;
	const typename table<T>::type *data_to_use = &data;
	typename table<T>::type new_data;

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
				make_array1d(new_data[i],new_data[0].size());
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
inline table_map<std::string>::type make_table_map(
		const table<std::string>::type & data,
		const header::type & header,
		const bool silent=false)
{
	return make_table_map<std::string>(data,header,silent);
}

// Splits a string into a vector of "word" strings on whitespace
std::vector< std::string > split_on_whitespace( const std::string & sentence );
header::type convert_to_header( const std::string & line );

#endif // End global function declarations

}

#endif // __BRG_FILE_FUNCTIONS_H_INCLUDED__
