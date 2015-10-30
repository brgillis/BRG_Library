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

#ifndef ICEBRG_MAIN_FILE_ACCESS_TABLE_FORMATS_FORMATTED_ASCII_FORMAT_HPP_
#define ICEBRG_MAIN_FILE_ACCESS_TABLE_FORMATS_FORMATTED_ASCII_FORMAT_HPP_

#include "IceBRG_main/common.h"
#include "IceBRG_main/file_access/ascii_table.hpp"
#include "IceBRG_main/file_access/table_formats/table_format_base.hpp"

namespace IceBRG {

template<typename T_value>
class formatted_ascii_format : public table_format_base<T_value>
{

public:

	typedef labeled_array<T_value> table_type;

	/// Destructor (default virtual)
	virtual ~formatted_ascii_format() noexcept override {}

	/// Name for this format (all in lower case)
	virtual str_type name() const override
	{
		return "formatted_ascii";
	}

	/// Read method
	virtual table_type read(std::istream & fi) const override
	{
		table_type table;
    	table.set_labels(load_header(fi));
    	table.set_rows(load_table<T_value>(fi,Eigen::RowMajor,
    			T_value(),
    			table.num_cols()));
    	return table;
	}

	/// Write method
	virtual void write(const labeled_array<T_value> & table, std::ostream & fo) const override
	{
		header_t header = table.get_labels();

		table_t<T_value> data;

		// Fill up the output table
		for( const auto & row : table.rows())
		{
			data.push_back(coerce<std::vector<T_value>>(row.raw()));
		}
		print_table(fo,data,header,Eigen::RowMajor);
	}

};

} // namespace IceBRG



#endif // ICEBRG_MAIN_FILE_ACCESS_TABLE_FORMATS_FORMATTED_ASCII_FORMAT_HPP_
