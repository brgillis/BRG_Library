/**********************************************************************\
 @file bin_table.cpp
 ------------------

 Source file containing implementations of functions defined in
 bin_table.h.

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

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "brg/global.h"

#include "bin_table.h"

void brgastro::open_bin_file( std::ofstream & stream, const std::string & name,
		const bool silent )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str(), std::ios::out | std::ios::binary );
	if ( !stream )
	{
		throw std::runtime_error("ERROR: Could not open file " + name + "!\n");
	}
}
void brgastro::open_bin_file( std::ifstream & stream, const std::string & name,
		const bool silent )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str(), std::ios::in | std::ios::binary );
	if ( !stream )
	{
		throw std::runtime_error("ERROR: Could not open file " + name + "!\n");
	}
}
void brgastro::open_bin_file( std::fstream & stream, const std::string & name,
		const bool silent )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str(), std::ios::out | std::ios::in | std::ios::binary  );
	if ( !stream )
	{
		throw std::runtime_error("ERROR: Could not open file " + name + "!\n");
	}
}


