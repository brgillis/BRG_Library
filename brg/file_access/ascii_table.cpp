/**********************************************************************\
 @file ascii_table.cpp
 ------------------

 Source file for functions defined in ascii_table.h.

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

#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/spirit/home/support/detail/hold_any.hpp>

#include "brg/global.h"

#include "brg/container/table_typedefs.hpp"
#include "brg/file_access/trim_comments.hpp"
#include "brg/vector/manipulations.hpp"

#include "ascii_table.h"


namespace brgastro {

header_t load_header( std::istream & table_stream,
		const bool silent)
{
	std::string temp_line;
	std::vector<header_t> possible_headers;

	// Get all comment lines at the top of the file
	while ( table_stream )
	{
		if ( table_stream.peek() == (int)( *"#" ) )
		{
			getline( table_stream, temp_line );

			header_t temp_header = convert_to_header(temp_line);

			if(temp_header.size() > 0)
				possible_headers.push_back(temp_header);
		}
		else
		{
			break;
		}
	}

	if(possible_headers.size()==1) return possible_headers[0];
	if(possible_headers.size()==0) return header_t();

	// If we get here, more than one line is a possible header candidate. Our next step is to
	// go to the data and count the columns in the first line. If only one possible header has
	// the right length, we know that's the one.

	unsigned short int n_cols = 0;
	do
	{
		getline( table_stream, temp_line );
		std::string junk;
		std::istringstream line_data_stream(temp_line);
		while (line_data_stream >> junk)
		{
			++n_cols;
		}
	} while(n_cols==0 && table_stream);

	if(n_cols==0) // If we can't find any data
	{
		std::cerr << "ERROR: Header line ambiguous; returning null header.\n";
		return header_t();
	}

	// Search through the possible headers, and see if we find exactly one with the right size
	unsigned short int num_right_size = 0;
	size_t i_best = 0;
	for(size_t i=0; i<possible_headers.size(); ++i)
	{
		if(possible_headers[i].size()==n_cols)
		{
			++num_right_size;
			i_best = i;
		}
	}

	if(num_right_size != 1) // If multiple or zero lines are the right size
	{
		std::cerr << "ERROR: Header line ambiguous; returning null header.\n";
		return header_t();
	}

	return possible_headers[i_best];
}

} // namespace brgastro

