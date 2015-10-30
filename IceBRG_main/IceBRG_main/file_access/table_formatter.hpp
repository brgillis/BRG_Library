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
#include <boost/algorithm/string.hpp>

#include "IceBRG_main/file_access/table_formats/ascii_format.hpp"
#include "IceBRG_main/file_access/table_formats/formatted_ascii_format.hpp"
#include "IceBRG_main/file_access/table_formats/table_format_base.hpp"
#include "IceBRG_main/common.h"

namespace IceBRG {

template<typename T_table>
std::unique_ptr<table_format_base<T_table>> get_formatter(str_type format)
{
	typedef std::unique_ptr<table_format_base<T_table>> p_formatter_t;
	boost::to_lower(format);

	if(format==ascii_format<T_table>().name())
		return p_formatter_t(new ascii_format<T_table>);
	if(format==formatted_ascii_format<T_table>().name())
		return p_formatter_t(new formatted_ascii_format<T_table>);

	throw std::runtime_error("Format " + format + " not recognized.");
}

} // namespace IceBRG


#endif // ICEBRG_MAIN_FILE_ACCESS_TABLE_FORMATTER_HPP_
