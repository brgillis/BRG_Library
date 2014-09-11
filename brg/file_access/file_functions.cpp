/**********************************************************************\
  @file file_functions.cpp

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

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "brg/global.h"

#include "brg/file_access/table_typedefs.hpp"

#include "file_functions.h"

using namespace std;

/** Global function implementations **/
#if (1)

void brgastro::open_file( std::ofstream & stream, const std::string & name,
		const bool silent )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str() );
	if ( !stream )
	{
		throw std::runtime_error("ERROR: Could not open file " + name + "!\n");
	}
}
void brgastro::open_file( std::ifstream & stream, const std::string & name,
		const bool silent )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str() );
	if ( !stream )
	{
		throw std::runtime_error("ERROR: Could not open file " + name + "!\n");
	}
}
void brgastro::open_file( std::fstream & stream, const std::string & name,
		const bool silent )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str() );
	if ( !stream )
	{
		throw std::runtime_error("ERROR: Could not open file " + name + "!\n");
	}
}

void brgastro::trim_comments_one_line( std::istream & stream,
		const bool silent )
{
	std::string file_data;
	if ( stream )
	{
		if ( stream.peek() == (int)( *"#" ) )
			getline( stream, file_data );
	}
}

void brgastro::trim_comments_one_line( std::fstream & stream,
		const bool silent )
{
	std::string file_data;
	if ( stream )
	{
		if ( stream.peek() == (int)( *"#" ) )
			getline( stream, file_data );
	}
}

void brgastro::trim_comments_all_at_top( std::istream & stream,
		const bool silent )
{
	std::string file_data;
	while ( stream )
	{
		if ( stream.peek() == (int)( *"#" ) )
		{
			getline( stream, file_data );
		}
		else
		{
			return;
		}
	}
}

void brgastro::trim_comments_all_at_top( std::fstream & stream,
		const bool silent )
{
	std::string file_data;
	while ( stream )
	{
		if ( stream.peek() == (int)( *"#" ) )
		{
			getline( stream, file_data );
		}
		else
		{
			return;
		}
	}
}

std::vector< std::string > brgastro::split_on_whitespace( const std::string & sentence )
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
brgastro::header::type brgastro::convert_to_header( const std::string & line )
{
	header::type result;
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

#endif // end global function implementations
