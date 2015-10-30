/**********************************************************************\
 @file ascii_format.hpp
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

#ifndef ICEBRG_MAIN_FILE_ACCESS_TABLE_FORMATS_ASCII_FORMAT_HPP_
#define ICEBRG_MAIN_FILE_ACCESS_TABLE_FORMATS_ASCII_FORMAT_HPP_

#include "IceBRG_main/common.h"
#include "IceBRG_main/file_access/ascii_table.hpp"
#include "IceBRG_main/file_access/table_formats/table_format_base.hpp"

namespace IceBRG {

template<typename T_table>
class ascii_format : public table_format_base<T_table>
{

public:

	/// Destructor (default virtual)
	virtual ~ascii_format() noexcept override {}

	/// Name for this format (all in lower case)
	virtual str_type name() const override
	{
		return "ascii";
	}

	/// Read method
	virtual T_table read(std::istream & fi) const override
	{
		T_table table;
    	table.set_labels(load_header(fi));
    	table.set_rows(load_table<typename T_table::value_type>(fi,Eigen::RowMajor,
    													typename T_table::value_type(),
    													table.num_cols()));
    	return table;
	}

	/// Write method
	virtual void write(const T_table & table, std::ostream & fo) const override
	{
    	header_t header = table.get_labels();

		// First, the header
		fo << "# ";
		for(const auto & label : header)
		{
			fo << label << "\t";
		}
		fo << std::endl;

		// And now the data table
		fo << table.data_table() << std::endl;

    	return;
	}

};

} // namespace IceBRG



#endif // ICEBRG_MAIN_FILE_ACCESS_TABLE_FORMATS_ASCII_FORMAT_HPP_
