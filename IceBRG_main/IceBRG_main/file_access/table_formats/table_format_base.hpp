/**********************************************************************\
 @file table_format_base.hpp
 ------------------

 TODO <Insert file description here>

 **********************************************************************

 Copyright (C) 2015 brg

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

#ifndef ICEBRG_MAIN_FILE_ACCESS_TABLE_FORMATS_TABLE_FORMAT_BASE_HPP_
#define ICEBRG_MAIN_FILE_ACCESS_TABLE_FORMATS_TABLE_FORMAT_BASE_HPP_

#include <string>
#include <iostream>
#include <Eigen/Core>

#include "IceBRG_main/common.h"

#include "IceBRG_main/container/labeled_array.hpp"

namespace IceBRG {

/**
 * An abstract base class representing a format in which a table can
 * be output.
 */
template<typename T_value>
class table_format_base
{
public:

	/// Destructor (default virtual)
	virtual ~table_format_base() noexcept {}

	/// Name for this format (all in lower case)
	virtual str_type name() const = 0;

	/// Read method to be implemented
	virtual labeled_array<T_value> read(std::istream &) const = 0;

	/// Write method to be implemented
	virtual void write(const labeled_array<T_value> &, std::ostream &) const = 0;

	/// Read from filename - delegates to stream version
	labeled_array<T_value> load(const str_type & file_name)
    {
    	std::ifstream fi;
    	open_file_input(fi,file_name);

    	return read(fi);
    }

	/// Write to filename - delegates to stream version
	void write(const labeled_array<T_value> & table, const str_type & file_name)
    {
    	std::ofstream fi;
    	open_file_output(fi,file_name);

    	write(table,fi);
    }
};

} // namespace IceBRG


#endif // ICEBRG_MAIN_FILE_ACCESS_TABLE_FORMATS_TABLE_FORMAT_BASE_HPP_
