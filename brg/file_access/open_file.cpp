/**********************************************************************\
 @file open_file.cpp
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

#include <fstream>
#include <stdexcept>
#include <string>

#include "brg/global.h"

namespace brgastro {

void open_file_input( std::ifstream & stream, const std::string & name )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str(), std::ios::in );
	if ( !stream )
	{
		throw std::runtime_error("Could not open file " + name + ".");
	}
}
void open_file_output( std::ofstream & stream, const std::string & name )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str(), std::ios::out );
	if ( !stream )
	{
		throw std::runtime_error("Could not open file " + name + ".");
	}
}
void open_file_io( std::fstream & stream, const std::string & name )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str(), std::ios::out | std::ios::out );
	if ( !stream )
	{
		throw std::runtime_error("Could not open file " + name + ".");
	}
}
void open_file_append( std::ofstream & stream, const std::string & name )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str(), std::ios::app );
	if ( !stream )
	{
		throw std::runtime_error("Could not open file " + name + ".");
	}
}

void open_bin_file_input( std::ifstream & stream, const std::string & name )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str(), std::ios::in | std::ios::binary );
	if ( !stream )
	{
		throw std::runtime_error("Could not open file " + name + ".");
	}
}
void open_bin_file_io( std::fstream & stream, const std::string & name )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str(), std::ios::out | std::ios::out | std::ios::binary );
	if ( !stream )
	{
		throw std::runtime_error("Could not open file " + name + ".");
	}
}
void open_bin_file_output( std::ofstream & stream, const std::string & name )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str(), std::ios::out | std::ios::binary );
	if ( !stream )
	{
		throw std::runtime_error("Could not open file " + name + ".");
	}
}
void open_bin_file_append( std::ofstream & stream, const std::string & name )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str(), std::ios::app | std::ios::binary );
	if ( !stream )
	{
		throw std::runtime_error("Could not open file " + name + ".");
	}
}


} // namespace brgastro
