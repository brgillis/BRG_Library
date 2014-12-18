/**********************************************************************\
 @file print_pixels.hpp
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

#ifndef _BRG_FILE_ACCESS_PRINT_PIXELS_HPP_INCLUDED_
#define _BRG_FILE_ACCESS_PRINT_PIXELS_HPP_INCLUDED_

#include <type_traits>

#include "brg/container/is_container.hpp"
#include "brg/file_access/ascii_table_map.hpp"
#include "brg/file_access/open_file.hpp"

namespace brgastro {

template<typename T, typename T_data=typename T::value_type,
typename std::enable_if<brgastro::is_const_container<T>::value,T>::type* = nullptr,
typename std::enable_if<!brgastro::is_const_container<typename T::value_type>::value,typename T::value_type>::type* = nullptr>
void print_pixels(std::ostream & out_stream, const T & data)
{
	table_map_t<T_data> table;

	for(size_t i=0; i<data.size(); ++i)
	{
		table["row"].push_back(static_cast<T_data>(i));
		table["value"].push_back(static_cast<T_data>(data[i]));
	}

	print_table_map(out_stream,table);
}
template<typename T, typename T_data=typename T::value_type,
typename std::enable_if<brgastro::is_const_container<T>::value,T>::type* = nullptr,
typename std::enable_if<!brgastro::is_const_container<typename T::value_type>::value,typename T::value_type>::type* = nullptr>
void print_pixels(const std::string & file_name, const T & data)
{
	std::ofstream fo;
	open_file_output(fo,file_name);
	print_pixels<T,T_data>(fo,data);
}

template<typename T, typename T_data=typename T::value_type::value_type,
typename std::enable_if<brgastro::is_const_container<T>::value,T>::type* = nullptr,
typename std::enable_if<brgastro::is_const_container<typename T::value_type>::value,typename T::value_type>::type* = nullptr,
typename std::enable_if<!brgastro::is_const_container<typename T::value_type::value_type>::value,typename T::value_type::value_type>::type* = nullptr>
void print_pixels(std::ostream & out_stream, const T & data)
{
	table_map_t<T_data> table;

	for(size_t i=0; i<data.size(); ++i)
	{
		for(size_t j=0; j<data[i].size(); ++j)
		{
			table["row"].push_back(static_cast<T_data>(i));
			table["col"].push_back(static_cast<T_data>(j));
			table["value"].push_back(static_cast<T_data>(data[i][j]));
		}
	}

	print_table_map(out_stream,table);
}
template<typename T, typename T_data=typename T::value_type::value_type,
typename std::enable_if<brgastro::is_const_container<T>::value,T>::type* = nullptr,
typename std::enable_if<brgastro::is_const_container<typename T::value_type>::value,typename T::value_type>::type* = nullptr,
typename std::enable_if<!brgastro::is_const_container<typename T::value_type::value_type>::value,typename T::value_type::value_type>::type* = nullptr>
void print_pixels(const std::string & file_name, const T & data)
{
	std::ofstream fo;
	open_file_output(fo,file_name);
	print_pixels<T,T_data>(fo,data);
}

template<typename T, typename T_data=typename T::value_type::value_type::value_type,
typename std::enable_if<brgastro::is_const_container<T>::value,T>::type* = nullptr,
typename std::enable_if<brgastro::is_const_container<typename T::value_type>::value,typename T::value_type>::type* = nullptr,
typename std::enable_if<brgastro::is_const_container<typename T::value_type::value_type>::value,typename T::value_type::value_type>::type* = nullptr,
typename std::enable_if<!brgastro::is_const_container<typename T::value_type::value_type::value_type>::value,typename T::value_type::value_type::value_type>::type* = nullptr>
void print_pixels(std::ostream & out_stream, const T & data)
{
	table_map_t<T_data> table;

	for(size_t i=0; i<data.size(); ++i)
	{
		for(size_t j=0; j<data[i].size(); ++j)
		{
			for(size_t k=0; j<data[i][j].size(); ++k)
			{
				table["row"].push_back(static_cast<T_data>(i));
				table["col"].push_back(static_cast<T_data>(j));
				table["layer"].push_back(static_cast<T_data>(k));
				table["value"].push_back(static_cast<T_data>(data[i][j][k]));
			}
		}
	}

	print_table_map(out_stream,table);
}
template<typename T, typename T_data=typename T::value_type::value_type::value_type,
typename std::enable_if<brgastro::is_const_container<T>::value,T>::type* = nullptr,
typename std::enable_if<brgastro::is_const_container<typename T::value_type>::value,typename T::value_type>::type* = nullptr,
typename std::enable_if<brgastro::is_const_container<typename T::value_type::value_type>::value,typename T::value_type::value_type>::type* = nullptr,
typename std::enable_if<!brgastro::is_const_container<typename T::value_type::value_type::value_type>::value,typename T::value_type::value_type::value_type>::type* = nullptr>
void print_pixels(const std::string & file_name, const T & data)
{
	std::ofstream fo;
	open_file_output(fo,file_name);
	print_pixels<T,T_data>(fo,data);
}

} // namespace brgastro

#endif // _BRG_FILE_ACCESS_PRINT_PIXELS_HPP_INCLUDED_
