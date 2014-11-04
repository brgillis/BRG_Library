/**********************************************************************\
 @file binary_archive.hpp
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

// body file: binary_archive.cpp

#ifndef _BRG_FILE_ACCESS_BINARY_ARCHIVE_HPP_INCLUDED_
#define _BRG_FILE_ACCESS_BINARY_ARCHIVE_HPP_INCLUDED_

#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "brg/file_access/open_file.hpp"

namespace brgastro {

template<typename T>
void binary_save(boost::archive::binary_oarchive &ar, const T & obj)
{
	ar << obj;
}
template<typename T>
void binary_save(std::ostream &out, const T & obj)
{
	boost::archive::binary_oarchive ar(out);
	ar << obj;
}
template<typename T>
void binary_save(const std::string & filename, const T & obj)
{
	std::ofstream out;
	open_file_output(out,filename);
	binary_save(out,obj);
}

// Vector<bool> overload
template<>
void binary_save(boost::archive::binary_oarchive &ar, const std::vector<bool> & obj)
{
	size_t size = obj.size();
	ar << size;

	size_t i;

	// Save most data
	for(i=0; i<obj.size()-8; i+=8)
	{
		unsigned char byte=0;
		for(int b=0;b<8;++b)
		{
			byte |= (obj[i+b] << b);
		}
		ar << byte;
	}

	// Save any remaining data
	unsigned char byte=0;
	for(unsigned char b=0; b < obj.size()-i; ++b)
	{
		byte |= (obj[i+b] << b);
	}
	ar << byte;
}
template<>
void binary_save(std::ostream &out, const std::vector<bool> & obj)
{
	boost::archive::binary_oarchive ar(out);
	return binary_save(ar,obj);
}
template<>
void binary_save(const std::string & filename, const std::vector<bool> & obj)
{
	std::ofstream out;
	open_file_output(out,filename);
	binary_save(out,obj);
}

// Vector overload
template<typename T>
void binary_save(boost::archive::binary_oarchive &ar, const std::vector<T> & obj)
{
	size_t size = obj.size();
	ar << size;
	for(const auto & elem : obj)
	{
		binary_save<T>(ar,elem);
	}
}
template<typename T>
void binary_save(std::ostream &out, const std::vector<T> & obj)
{
	boost::archive::binary_oarchive ar(out);
	binary_save(ar,obj);
}
template<typename T>
void binary_save(const std::string & filename, const std::vector<T> & obj)
{
	std::ofstream out;
	open_file_output(out,filename);
	binary_save(out,obj);
}

template<typename T>
T binary_load(boost::archive::binary_iarchive & ar)
{
	T obj;
	ar >> obj;
	return obj;
}
template<typename T>
T binary_load(std::istream & in)
{
	boost::archive::binary_iarchive ar(in);
	return binary_load<T>(ar);
}
template<typename T>
T binary_load(const std::string & filename)
{
	std::ifstream in;
	open_file_input(in,filename);
	return binary_load<T>(in);
}

// Vector<bool> overload
template<>
std::vector<bool> binary_load(boost::archive::binary_iarchive & ar)
{
	size_t size;
	ar >> size;
	std::vector<bool> result(size);

	size_t i;
	// Load most of the data
	for(i=0; i<size-8; i+=8)
	{
		unsigned char byte;
		ar >> byte;
		for(int b=0; b<8; ++b)
		{
			result[i+b] = byte & (1 << b);
		}
	}

	// Load the remainder of the data
	unsigned char byte;
	ar >> byte;
	for(unsigned char b=0; b<size-i; ++b)
	{
		result[i+b] = byte & (1 << b);
	}

	return result;
}
template<>
std::vector<bool> binary_load(std::istream &in)
{
	boost::archive::binary_iarchive ar(in);
	return binary_load<std::vector<bool>>(ar);
}
template<>
std::vector<bool> binary_load(const std::string & filename)
{
	std::ifstream in;
	open_file_input(in,filename);
	return binary_load<std::vector<bool>>(in);
}

// Vector overload
template<typename T>
T binary_load_vector(boost::archive::binary_iarchive &ar)
{
	size_t size;
	ar >> size;
	T result(size);

	for(size_t i=0; i<size; ++i)
	{
		result[i] = binary_load<typename T::value_type>(ar);
	}

	return result;
}
template<typename T>
T binary_load_vector(std::istream &in)
{
	boost::archive::binary_iarchive ar(in);
	return binary_load_vector<T>(ar);
}
template<typename T>
T binary_load_vector(const std::string & filename)
{
	std::ifstream in;
	open_file_input(in,filename);
	return binary_load_vector<T>(in);
}

} // end namespace brgastro


#endif // _BRG_FILE_ACCESS_BINARY_ARCHIVE_HPP_INCLUDED_
