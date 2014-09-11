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

// Splits a string into a vector of "word" strings on whitespace
std::vector< std::string > split_on_whitespace( const std::string & sentence );
header::type convert_to_header( const std::string & line );

#endif // End global function declarations

}

#endif // __BRG_FILE_FUNCTIONS_H_INCLUDED__
