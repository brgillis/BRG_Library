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

// body file: brg/file_functions.cpp

#ifndef _BRG_FILE_FUNCTIONS_H_INCLUDED_
#define _BRG_FILE_FUNCTIONS_H_INCLUDED_

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

#include "global.h"

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
void open_bin_file( std::ofstream & stream, const std::string & name,
		const bool silent = false );
void open_bin_file( std::ifstream & stream, const std::string & name,
		const bool silent = false );
void open_bin_file( std::fstream & stream, const std::string & name,
		const bool silent = false );

// Prints a formatted table in the passed stream. header is a vector of strings representing the labels for each column,
// and data is a 2-d vector of the data to be printed, in the format data[c][r], where c is the column index and r is the row index.
void print_table( std::ostream & out_stream,
		const std::vector< std::vector< std::string > > & data,
		const std::vector< std::string > & header = std::vector< std::string >(0),
		const bool silent = false );

// Some templates to coerce any type of data to be printed out
template<typename T>
void print_table( std::ostream & out_stream,
		const std::vector< std::vector< T > > & data,
		const std::vector< std::string > & header = std::vector< std::string >(0),
		const bool silent = false )
{
	std::stringstream ss;
	std::vector< std::vector<std::string> > string_data(0);

	make_array(string_data, data.size(), data.at(0).size());
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
		const std::vector< std::vector< std::string > > & data,
		const std::vector< std::string > & header = std::vector< std::string >(0),
		const bool silent = false )
{
	std::ofstream fo;
	open_file(fo,file_name,silent);

	print_table(fo,data,header,silent);
}
template<typename T>
void print_table( const std::string & file_name,
		const std::vector< std::vector< T > > & data,
		const std::vector< std::string > & header = std::vector< std::string >(0),
		const bool silent = false )
{
	std::ofstream of;
	open_file(of,file_name,silent);

	print_table<T>(of,data,header,silent);
}

// Load table, either loading in the entire table, or only loading in certain columns into pointed-to
// variables, found by matching header entries to the strings passed
std::vector<std::vector<std::string> > load_table( std::istream & table_file_name,
		const bool silent=false);

template<typename T>
std::vector<std::vector<T> > load_table( std::istream & table_file_name,
		const bool silent=false)
{
	std::stringstream ss;

	std::vector< std::vector<std::string> > string_data(load_table( table_file_name, silent ));
	std::vector< std::vector<T> > data;

	make_array2d(data,string_data.size(),string_data.at(0).size());

	for(size_t i=0; i<string_data.size(); ++i)
	{
		for(size_t j=0; j<string_data[i].size(); ++j)
		{
			ss.str(string_data[i][j]);
			ss >> data[i].at(j);
		}
	}
	return data;
}

// And to allow us to load from a file name instead of a stream
inline std::vector<std::vector<std::string> > load_table( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_table(fi,silent);
}
template<typename T>
inline std::vector<std::vector<T> > load_table( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_table<T>(fi,silent);
}

void load_table_and_header( std::istream & fi,
		std::vector<std::vector<std::string> > & table_data,
		std::vector<std::string> & header, const bool silent=false);

template<typename T>
void load_table_and_header( std::istream & fi,
		std::vector<std::vector<T> > & table_data,
		std::vector<std::string> & header, const bool silent=false)
{
	std::stringstream ss;

	std::vector< std::vector<std::string> > string_data;
	load_table_and_header( fi, string_data, header, silent );

	make_array2d(table_data,string_data.size(),string_data.at(0).size());

	for(size_t i=0; i<string_data.size(); ++i)
	{
		for(size_t j=0; j<string_data[i].size(); ++j)
		{
			ss.str(string_data[i][j]);
			ss >> table_data[i].at(j);
		}
	}
}

// And to allow us to load from a file name instead of a stream
inline void load_table_and_header( const std::string & file_name,
		std::vector<std::vector< std::string > > & table_data,
		std::vector<std::string> & header, const bool silent = false )
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_table_and_header(fi,table_data,header,silent);
}
template<typename T>
void load_table_and_header( const std::string & file_name,
		std::vector<std::vector<T> > & table_data,
		std::vector<std::string> & header, const bool silent = false )
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_table_and_header<T>(fi,table_data,header,silent);
}

void load_table_columns( std::istream & fi,
		std::vector< std::pair< std::string, std::vector<std::string>* > > & header_links,
		const bool case_sensitive=false, const bool silent=false);

template<typename T>
void load_table_columns( std::istream & fi,
		std::vector< std::pair< std::string, std::vector<T>* > > & header_links,
		const bool case_sensitive=false, const bool silent=false)
{
	std::stringstream ss;

	std::vector< std::pair< std::string, std::vector<std::string>* > > string_header_links;
	std::vector< std::vector< std::string > > string_data(header_links.size());
	for(size_t i=0; i<header_links.size(); ++i)
	{
		string_header_links.push_back(std::make_pair(
				header_links[i].first,&(string_data[i])));
	}

	load_table_columns( fi, string_header_links, case_sensitive, silent);

	for(size_t i=0; i<header_links.size(); ++i)
	{
		header_links[i].second->resize(string_data[i].size());
		for(size_t j=0; j<string_data[i].size(); ++j)
		{
			ss.str("");
			ss.clear();
			ss.str(string_data[i][j]);
			ss >> header_links[i].second->at(j);
		}
	}
}

// And to allow us to load from a file name instead of a stream
inline void load_table_columns( const std::string & file_name,
		std::vector< std::pair< std::string, std::vector<std::string>* > > & header_links,
		const bool case_sensitive=false, const bool silent=false)
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_table_columns(fi,header_links,case_sensitive,silent);
}
template<typename T>
void load_table_columns( const std::string & file_name,
		std::vector< std::pair< std::string, std::vector<T>* > > & header_links,
		const bool case_sensitive=false, const bool silent=false)
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_table_columns<T>(fi,header_links,case_sensitive,silent);
}

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

#endif // End global function declarations

}

#endif // __BRG_FILE_FUNCTIONS_H_INCLUDED__
