/**********************************************************************\
 @file bin_table.h
 ------------------

 Functions for dealing with binary data tables.

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


// body file: bin_table.cpp


#ifndef _BRG_BIN_TABLE_H_INCLUDED_
#define _BRG_BIN_TABLE_H_INCLUDED_

namespace brgastro {

void open_bin_file( std::ofstream & stream, const std::string & name,
		const bool silent = false );
void open_bin_file( std::ifstream & stream, const std::string & name,
		const bool silent = false );
void open_bin_file( std::fstream & stream, const std::string & name,
		const bool silent = false );

} // namespace brgastro

#endif // _BRG_BIN_TABLE_H_INCLUDED_
