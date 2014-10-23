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

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "brg/file_access/open_file.hpp"

namespace brgastro {

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
template<typename T>
T binary_load(std::istream & in)
{
	T obj;
	boost::archive::binary_iarchive ar(in);
	ar >> obj;
	return obj;
}
template<typename T>
T binary_load(const std::string & filename)
{
	std::ifstream in;
	open_file_input(in,filename);
	return binary_load<T>(in);
}

} // end namespace brgastro


#endif // _BRG_FILE_ACCESS_BINARY_ARCHIVE_HPP_INCLUDED_
