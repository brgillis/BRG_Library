/**********************************************************************\
brg_file_functions.h
 -----------

 If this header is used, the source file brg_file_functions.cpp must be
 included and compiled with the project.

 This file includes various functions having to do with file reading and
 writing.

 Everything in this file is declared in the namespace brgastro.

 \**********************************************************************/

#ifndef __BRG_FILE_FUNCTIONS_H_INCLUDED__
#define __BRG_FILE_FUNCTIONS_H_INCLUDED__

#include "brg_global.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>

#include "brg_units.h"
#include "brg_misc_functions.hpp"

namespace brgastro
{

/** Global function declarations **/
#if (1)

// Prints a formatted table in the passed stream. header is a vector of strings representing the labels for each column,
// and data is a 2-d vector of the data to be printed, in the format data[c][r], where c is the column index and r is the row index.
// Prints an error message and returns 1 if an error arises (for instance, vectors are too short).
void _print_table( std::ostream & out_stream,
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
	for( unsigned int i=0; i<data.size(); ++i )
	{
		for( unsigned int j=0; i<data[0].size(); ++j )
		{
			ss.str("");
			ss << data[i][j];
			string_data[i].at(j) = ss.str();
		}
	}
	_print_table(out_stream, string_data, header, silent);
}
void print_table( std::ostream & out_stream,
		const std::vector< std::vector< std::string > > & data,
		const std::vector< std::string > & header = std::vector< std::string >(0),
		const bool silent = false )
{
	_print_table(out_stream, data, header, silent);
}

// Load table, either loading in the entire table, or only loading in certain columns into pointed-to
// variables, found by matching header entries to the strings passed
std::vector<std::vector<std::string> > _load_table( const std::string & table_file_name,
		const bool silent=false);

template<typename T>
std::vector<std::vector<T> > load_table( const std::string & table_file_name,
		const bool silent=false)
{
	std::stringstream ss;

	std::vector< std::vector<std::string> > string_data(_load_table( table_file_name, silent ));
	std::vector< std::vector<T> > data;

	make_array2d(data,string_data.size(),string_data.at(0).size());

	for(unsigned int i=0; i<string_data.size(); ++i)
	{
		for(unsigned int j=0; j<string_data[i].size(); ++j)
		{
			ss.str() = string_data[i][j];
			ss >> data[i].at(j);
		}
	}
	return data;
}
template<> //TODO: Can't specialize functions, overload instead
std::vector<std::vector<std::string> > load_table( const std::string & table_file_name,
		const bool silent)
{
	return _load_table( table_file_name, silent );
}

void _load_table_and_header( const std::string & table_file_name,
		std::vector<std::vector<std::string> > & table_data,
		std::vector<std::string> & header, const bool silent=false);

template<typename T>
void load_table_and_header( const std::string & table_file_name,
		std::vector<std::vector<T> > & table_data,
		std::vector<std::string> & header, const bool silent=false)
{
	std::stringstream ss;

	std::vector< std::vector<std::string> > string_data;
	_load_table_and_header( table_file_name, string_data, header, silent );

	make_array2d(table_data,string_data.size(),string_data.at(0).size());

	for(unsigned int i=0; i<string_data.size(); ++i)
	{
		for(unsigned int j=0; j<string_data[i].size(); ++j)
		{
			ss.str() = string_data[i][j];
			ss >> table_data[i].at(j);
		}
	}
}
template<>
void load_table_and_header( const std::string & table_file_name,
		std::vector<std::vector<std::string> > & table_data,
		std::vector<std::string> & header, const bool silent)
{
	_load_table_and_header( table_file_name, table_data, header, silent );
}

void _load_table_columns( const std::string & table_file_name,
		std::vector< std::pair< std::string, std::vector<std::string>* > > & header_links,
		const bool case_sensitive=false, const bool silent=false);

template<typename T>
void load_table_columns( const std::string & table_file_name,
		std::vector< std::pair< std::string, std::vector<T>* > > & header_links,
		const bool case_sensitive=false, const bool silent=false)
{
	std::stringstream ss;

	std::vector< std::pair< std::string, std::vector<std::string>* > > string_header_links;
	std::vector< std::vector< std::string > > string_data(header_links.size());
	for(unsigned int i=0; i<header_links.size(); ++i)
	{
		string_header_links.push_back(std::make_pair(
				header_links.first(),&(string_data[i])));
	}

	_load_table_columns( table_file_name, string_header_links, case_sensitive, silent);

	for(unsigned int i=0; i<header_links.size(); ++i)
	{
		header_links[i].second()->resize(string_data[i].size());
		for(unsigned int j=0; j<string_data[i].size(); ++j)
		{
			ss.str() = string_data[i][j];
			ss >> header_links[i].second()->at(j);
		}
	}
}
template<>
void load_table_columns( const std::string & table_file_name,
		std::vector< std::pair< std::string, std::vector<std::string>* > > & header_links,
		const bool case_sensitive, const bool silent)
{
	_load_table_columns( table_file_name, header_links, case_sensitive, silent);
}

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

// Functions to get rid of comments lines (those starting with #) from open fstreams and ifstreams. The "one_line" versions will
// trim the next line if it's commented out, and the "all_at_top" versions will trim all commented lines from the current
// position until they find a line which isn't a comment.
// The "one_line" functions will return 1 if the file is already at the end, and 0 otherwise.
// The "all_at_top" functions will return 1 if they reach the end of the file before the run out of comments, and 0 otherwise.
void trim_comments_one_line( std::ifstream & stream, const bool silent =
		false );
void trim_comments_one_line( std::fstream & stream, const bool silent =
		false );
void trim_comments_all_at_top( std::ifstream & stream, const bool silent =
		false );
void trim_comments_all_at_top( std::fstream & stream, const bool silent =
		false );

#endif // End global function declarations

}

#endif // __BRG_FILE_FUNCTIONS_H_INCLUDED__
