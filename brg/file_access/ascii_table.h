/**********************************************************************\
 @file ascii_table.h
 ------------------

 This header file contains definitions of various functions used to
 read from and write to ASCII data tables.

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


// body file: ascii_table.cpp

#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <boost/spirit/home/support/detail/hold_any.hpp>

#include "brg/global.h"

#include "brg/file_access/open_file.h"
#include "brg/file_access/table_typedefs.hpp"
#include "brg/file_access/table_utility.h"
#include "brg/utility.hpp"


#ifndef _BRG_ASCII_TABLE_H_INCLUDED_
#define _BRG_ASCII_TABLE_H_INCLUDED_

namespace brgastro {

// Prints a formatted table in the passed stream. header is a vector of strings representing the labels for each column,
// and data is a 2-d vector of the data to be printed, in the format data[c][r], where c is the column index and r is the row index.
void print_table( std::ostream & out_stream,
		const table_t<boost::spirit::hold_any> & data,
		const header_t & header = header_t(),
		const bool silent = false );

// Some templates to coerce any table of data to be printed out
template<typename T>
void print_table( std::ostream & out_stream,
		const table_t<T> & data,
		const header_t & header = header_t(),
		const bool silent = false )
{
	std::stringstream ss;
	table_t<boost::spirit::hold_any> any_data(0);

	make_array2d(any_data, data);
	for( size_t i=0; i<data.size(); ++i )
	{
		for( size_t j=0; j<data[i].size(); ++j )
		{
			any_data[i][j] = data[i][j];
		}
	}
	print_table(out_stream, any_data, header, silent);
}

// And to allow us to print to a file name instead of a stream
inline void print_table( const std::string & file_name,
		const table_t<boost::spirit::hold_any> & data,
		const header_t & header = header_t(),
		const bool silent = false )
{
	std::ofstream fo;
	open_file_output(fo,file_name);

	print_table(fo,data,header,silent);
}
template<typename T>
void print_table( const std::string & file_name,
		const table_t<T> & data,
		const header_t & header = header_t(),
		const bool silent = false )
{
	std::ofstream fo;
	open_file_output(fo,file_name);

	print_table<T>(fo,data,header,silent);
}

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
inline void print_table_map( std::ostream & out,
		const table_map_t<boost::spirit::hold_any> & table_map,
		const bool silent = false)
{
	header_t header;
	table_t<boost::spirit::hold_any> data;

	for(auto it=table_map.begin(); it!=table_map.end(); ++it)
	{
		header.push_back(it->first);
		data.push_back(it->second);
	}
	print_table(out,data,header,silent);
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
inline void print_table_map( const std::string & file_name,
		const table_map_t<boost::spirit::hold_any> & table_map,
		const bool silent = false )
{
	std::ofstream fo;
	open_file_output(fo,file_name);

	print_table_map(fo,table_map,silent);
}

// Load table, either loading in the entire table, or only loading in certain columns into pointed-to
// variables, found by matching header entries to the strings passed
table_t<boost::spirit::hold_any> load_table( std::istream & table_file_name,
		const bool silent=false);

template<typename T>
table_t<T> load_table( std::istream & table_file_name,
		const bool silent=false)
{
	table_t<boost::spirit::hold_any> any_data(load_table( table_file_name, silent ));
	table_t<T> data;

	make_array2d(data,any_data);

	for(size_t i=0; i<any_data.size(); ++i)
	{
		for(size_t j=0; j<any_data[i].size(); ++j)
		{
			data[i][j] = boost::spirit::any_cast<T>(any_data[i][j]);
		}
	}
	return data;
}

// And to allow us to load from a file name instead of a stream
inline table_t<boost::spirit::hold_any> load_table( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file_input(fi,file_name);

	return load_table(fi,silent);
}
template<typename T>
inline table_t<T> load_table( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file_input(fi,file_name);

	return load_table<T>(fi,silent);
}

// Load a table's header as a vector of strings
header_t load_header( std::istream & table_stream,
		const bool silent=false);

// And to allow us to load from a file name instead of a stream
inline header_t load_header( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file_input(fi,file_name);

	return load_header(fi,silent);
}

// Directly load a map of the data
inline table_map_t<boost::spirit::hold_any> load_table_map( std::istream & fi,
		const bool silent=false)
{
	header_t header = load_header(fi,silent);
	table_t<boost::spirit::hold_any> data = load_table(fi,silent);
	return make_table_map<boost::spirit::hold_any>(data,header,silent);
}

template<typename T>
table_map_t<T> load_table_map( std::istream & fi,
		const bool silent=false)
{
	header_t header = load_header(fi,silent);
	table_t<T> data = load_table<T>(fi,silent);
	return make_table_map<T>(data,header,silent);
}

// And to allow us to load from a file name instead of a stream
inline table_map_t<boost::spirit::hold_any> load_table_map( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file_input(fi,file_name);

	return load_table_map(fi,silent);
}
template<typename T>
table_map_t<T> load_table_map( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file_input(fi,file_name);

	return load_table_map<T>(fi,silent);
}

inline void load_table_columns( std::istream & fi,
		std::map< std::string, std::vector<boost::spirit::hold_any>* > & column_map,
		const bool case_sensitive=false, const bool silent=false)
{
	table_map_t<boost::spirit::hold_any> table_map = load_table_map(fi);

	for(std::map< std::string, std::vector<boost::spirit::hold_any>* >::iterator it=column_map.begin();
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

template<typename T>
void load_table_columns( std::istream & fi,
		std::map< std::string, std::vector<T>* > & column_map,
		const bool case_sensitive=false, const bool silent=false)
{
	table_map_t<T> table_map = load_table_map<T>(fi);

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
inline void load_table_columns( const std::string & file_name,
		std::map< std::string, std::vector<boost::spirit::hold_any>* > & column_map,
		const bool case_sensitive=false, const bool silent=false)
{
	std::ifstream fi;
	open_file_input(fi,file_name);

	load_table_columns(fi,column_map,case_sensitive,silent);
}
template<typename T>
void load_table_columns( const std::string & file_name,
		std::map< std::string, std::vector<T>* > & column_map,
		const bool case_sensitive=false, const bool silent=false)
{
	std::ifstream fi;
	open_file_input(fi,file_name);

	load_table_columns<T>(fi,column_map,case_sensitive,silent);
}

} // namespace brgastro

#endif // _BRG_ASCII_TABLE_H_INCLUDED_
