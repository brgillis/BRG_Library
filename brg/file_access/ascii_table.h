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

#include "brg/global.h"

#include "brg/file_access/file_functions.h"
#include "brg/file_access/table_typedefs.hpp"
#include "brg/utility.hpp"


#ifndef _BRG_ASCII_TABLE_H_INCLUDED_
#define _BRG_ASCII_TABLE_H_INCLUDED_

namespace brgastro {

// Prints a formatted table in the passed stream. header is a vector of strings representing the labels for each column,
// and data is a 2-d vector of the data to be printed, in the format data[c][r], where c is the column index and r is the row index.
void print_table( std::ostream & out_stream,
		const table<std::string>::type & data,
		const header::type & header = header::type(),
		const bool silent = false );

// Some templates to coerce any type of data to be printed out
template<typename T>
void print_table( std::ostream & out_stream,
		const typename table<T>::type & data,
		const header::type & header = header::type(),
		const bool silent = false )
{
	std::stringstream ss;
	table<std::string>::type string_data(0);

	make_array2d(string_data, data);
	for( size_t i=0; i<data.size(); ++i )
	{
		for( size_t j=0; i<data[i].size(); ++j )
		{
			ss.str("");
			ss << data[i][j];
			string_data[i].at(j) = ss.str();
		}
	}
	print_table(out_stream, string_data, header, silent);
}

// And to allow us to print to a file name instead of a stream
inline void print_table( const std::string & file_name,
		const table<std::string>::type & data,
		const header::type & header = header::type(),
		const bool silent = false )
{
	std::ofstream fo;
	open_file(fo,file_name,silent);

	print_table(fo,data,header,silent);
}
template<typename T>
void print_table( const std::string & file_name,
		const typename table<T>::type & data,
		const header::type & header = header::type(),
		const bool silent = false )
{
	std::ofstream of;
	open_file(of,file_name,silent);

	print_table<T>(of,data,header,silent);
}

// Print a table from a map
template<typename T>
void print_table_map( std::ostream & out,
		const typename table_map<T>::type & table_map,
		const bool silent = false)
{
	header::type header;
	typename table<T>::type data;

	for(auto it=table_map.begin(); it!=table_map.end(); ++it)
	{
		header.push_back(it->first);
		data.push_back(it->second);
	}
	print_table<T>(out,data,header,silent);
}
inline void print_table_map( std::ostream & out,
		const table_map<std::string>::type & table_map,
		const bool silent = false)
{
	header::type header;
	table<std::string>::type data;

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
		const typename table_map<T>::type & table_map,
		const bool silent = false )
{
	std::ofstream fo;
	open_file(fo,file_name,silent);

	print_table_map<T>(fo,table_map,silent);
}
void print_table_map( const std::string & file_name,
		const table_map<std::string>::type & table_map,
		const bool silent = false )
{
	std::ofstream of;
	open_file(of,file_name,silent);

	print_table_map(of,table_map,silent);
}

// Load table, either loading in the entire table, or only loading in certain columns into pointed-to
// variables, found by matching header entries to the strings passed
table<std::string>::type load_table( std::istream & table_file_name,
		const bool silent=false);

template<typename T>
typename table<T>::type load_table( std::istream & table_file_name,
		const bool silent=false)
{
	std::stringstream ss;

	table<std::string>::type string_data(load_table( table_file_name, silent ));
	typename table<T>::type data;

	make_array2d(data,string_data);

	for(size_t i=0; i<string_data.size(); ++i)
	{
		for(size_t j=0; j<string_data[i].size(); ++j)
		{
			ss.clear();
			ss.str(string_data[i][j]);
			ss >> data[i].at(j);
		}
	}
	return data;
}

// And to allow us to load from a file name instead of a stream
inline table<std::string>::type load_table( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_table(fi,silent);
}
template<typename T>
inline typename table<T>::type load_table( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_table<T>(fi,silent);
}

// Load a table's header as a vector of strings
header::type load_header( std::istream & table_stream,
		const bool silent=false);

// And to allow us to load from a file name instead of a stream
inline header::type load_header( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_header(fi,silent);
}

// Directly load a map of the data
inline table_map<std::string>::type load_table_map( std::istream & fi,
		const bool silent=false)
{
	header::type header = load_header(fi,silent);
	table<std::string>::type data = load_table(fi,silent);
	return make_table_map(data,header,silent);
}

template<typename T>
typename table_map<T>::type load_table_map( std::istream & fi,
		const bool silent=false)
{
	header::type header = load_header(fi,silent);
	typename table<T>::type data = load_table<T>(fi,silent);
	return make_table_map<T>(data,header,silent);
}

// And to allow us to load from a file name instead of a stream
inline table_map<std::string>::type load_table_map( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_table_map(fi,silent);
}
template<typename T>
typename table_map<T>::type load_table_map( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_table_map<T>(fi,silent);
}

inline void load_table_columns( std::istream & fi,
		std::map< std::string, std::vector<std::string>* > & column_map,
		const bool case_sensitive=false, const bool silent=false)
{
	table_map<std::string>::type table_map = load_table_map(fi);

	for(std::map< std::string, std::vector<std::string>* >::iterator it=column_map.begin();
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
	typename table_map<T>::type table_map = load_table_map<T>(fi);

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
		std::map< std::string, std::vector<std::string>* > & column_map,
		const bool case_sensitive=false, const bool silent=false)
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	load_table_columns(fi,column_map,case_sensitive,silent);
}
template<typename T>
void load_table_columns( const std::string & file_name,
		std::map< std::string, std::vector<T>* > & column_map,
		const bool case_sensitive=false, const bool silent=false)
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	load_table_columns<T>(fi,column_map,case_sensitive,silent);
}

} // namespace brgastro

#endif // _BRG_ASCII_TABLE_H_INCLUDED_
