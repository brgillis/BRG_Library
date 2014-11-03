/**********************************************************************\
 @file ascii_table_map.hpp
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


#ifndef _BRG_ASCII_TABLE_MAP_HPP_INCLUDED_
#define _BRG_ASCII_TABLE_MAP_HPP_INCLUDED_

#include <iostream>
#include <map>
#include <vector>

#include "brg/global.h"

#include "brg/file_access/ascii_table.h"
#include "brg/file_access/table_typedefs.hpp"
#include "brg/file_access/table_utility.h"

namespace brgastro {

// Print a table from a map
template<typename T>
void print_table_map( std::ostream & out,
		const table_map_t<T> & table_map,
		const bool silent = false)
{
	header_t header;
	table_t<T> data;

	for(auto it=table_map.begin(); it!=table_map.end(); ++it)
	{
		header.push_back(it->first);
		data.push_back(it->second);
	}
	print_table<T>(out,data,header,silent);
}

// And to allow us to print to a file name instead of a stream
template<typename T>
inline void print_table_map( const std::string & file_name,
		const table_map_t<T> & table_map,
		const bool silent = false )
{
	std::ofstream fo;
	open_file_output(fo,file_name);

	print_table_map<T>(fo,table_map,silent);
}

// Directly load a map of the data
template<typename T>
table_map_t<T> load_table_map( std::istream & fi,
		const bool silent=false, const T default_value=T())
{
	header_t header = load_header(fi,silent);
	table_t<T> data = load_table<T>(fi,silent,default_value);
	return make_table_map<T>(data,header,silent);
}

// And to allow us to load from a file name instead of a stream
template<typename T>
table_map_t<T> load_table_map( const std::string & file_name,
		const bool silent = false, const T default_value=T() )
{
	std::ifstream fi;
	open_file_input(fi,file_name);

	return load_table_map<T>(fi,silent,default_value);
}
template<typename T>
void load_table_columns( std::istream & fi,
		std::map< std::string, std::vector<T>* > & column_map,
		const bool case_sensitive=false, const bool silent=false, const T default_value=T())
{
	table_map_t<T> table_map = load_table_map<T>(fi,silent,default_value);

	for(typename std::map< std::string, std::vector<T>* >::iterator it=column_map.begin();
			it!=column_map.end(); ++it)
	{
		*(it->second) = table_map[it->first];

		// Check that we found it
		if(!silent)
		{
			if(it->second->size()==0)
			{
				std::cerr << "WARNING: Column " << it->first << " not found in table.\n";
			}
		}

	}
}

// And to allow us to load from a file name instead of a stream
template<typename T>
void load_table_columns( const std::string & file_name,
		std::map< std::string, std::vector<T>* > & column_map,
		const bool case_sensitive=false, const bool silent=false, const T default_value=T())
{
	std::ifstream fi;
	open_file_input(fi,file_name);

	load_table_columns<T>(fi,column_map,case_sensitive,silent,default_value);
}

} // namespace brgastro

#endif // _BRG_ASCII_TABLE_MAP_HPP_INCLUDED_