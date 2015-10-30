/**********************************************************************\
 @file table_formatter.hpp
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

#ifndef ICEBRG_MAIN_FILE_ACCESS_TABLE_FORMATTER_HPP_
#define ICEBRG_MAIN_FILE_ACCESS_TABLE_FORMATTER_HPP_

#include <memory>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "IceBRG_main/file_access/table_formats/ascii_format.hpp"
#include "IceBRG_main/file_access/table_formats/formatted_ascii_format.hpp"
#include "IceBRG_main/file_access/table_formats/table_format_base.hpp"
#include "IceBRG_main/common.h"

namespace IceBRG {

template<typename T_value>
using p_formatter_t = std::unique_ptr<table_format_base<T_value>>;

template<typename T_value>
struct table_formatters
{
	static const std::initializer_list<p_formatter_t<T_value>> formatters;

	static const table_format_base<T_value> & get_formatter(str_type format)
	{
		boost::to_lower(format);

		for(const auto & formatter : formatters)
		{
			if(format==formatter->name())
			{
				return *formatter;
			}
		}

		throw std::runtime_error("Format " + format + " not recognized.");
	}
};

template<typename T_value>
const std::initializer_list<p_formatter_t<T_value>> table_formatters<T_value>::formatters({
	p_formatter_t<T_value>(new ascii_format<T_value>),
	p_formatter_t<T_value>(new formatted_ascii_format<T_value>),
});

} // namespace IceBRG


#endif // ICEBRG_MAIN_FILE_ACCESS_TABLE_FORMATTER_HPP_
