/**********************************************************************\
 @file table_utility.cpp
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

#include <sstream>
#include <string>
#include <vector>

#include "brg/global.h"

#include "brg/file_access/table_typedefs.hpp"

#include "table_utility.h"

namespace brgastro {

std::vector< std::string > split_on_whitespace( const std::string & sentence )
{
	std::vector< std::string > result;
	std::istringstream sentence_data_stream(sentence);

	std::string word;
	while (sentence_data_stream >> word)
	{
		result.push_back(word);
	}

	return result;
}
header_t convert_to_header( const std::string & line )
{
	header_t result;
	std::istringstream line_data_stream(line);

	std::string word;

	// Get rid of first word if it's the comment indicator
	if ( line_data_stream.peek() == (int)( *"#" ) )
	{
		line_data_stream >> word;
	}

	while (line_data_stream >> word)
	{

		result.push_back(word);
	}

	return result;
}

}
