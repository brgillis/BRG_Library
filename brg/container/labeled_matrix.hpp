/**********************************************************************\
 @file labeled_matrix.hpp
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

#ifndef _BRG_CONTAINER_LABELED_MATRIX_HPP_INCLUDED_
#define _BRG_CONTAINER_LABELED_MATRIX_HPP_INCLUDED_

#include <string>
#include <unordered_map>
#include <utility>

#include <Eigen/Array>

#include "brg/container/insertion_ordered_map.hpp"
#include "brg/container/labeled_matrix_row.hpp"
#include "brg/container/table_typedefs.hpp"

namespace brgastro {

template<typename data_type, typename key_type=std::string>
class labeled_matrix
{
public:

	// Public typedefs
	typedef typename Eigen::Matrix<data_type,Eigen::Dynamic,Eigen::Dynamic> data_type;
	typedef data_type::Index size_type;
	typedef data_type::ColXpr column_type;
	typedef data_type::ConstColXpr const_column_type;

private:

	// Private typedefs
	typedef typename std::unordered_map<key_type,size_type> map_type;
	typedef typename brgastro::table_map_t<data_type,key_type> buffer_type;

	mutable data_type _data_table_;
	mutable map_type _key_map_;
	mutable buffer_type _buffer_;

	void _add_buffer_to_data_table() const
	{
		// TODO fill in
	}

public:

	// Constructors
#if(1)

	/// Default constructor
	labeled_matrix() {}

	/// Copy from table_map
	labeled_matrix(const buffer_type & init_buffer)
	: _buffer_(init_buffer)
	{
	}

	/// Move from table map
	labeled_matrix(buffer_type && init_buffer)
	: _buffer_(init_buffer)
	{
	}

#endif

	// Data access
	const data_type & data_table() const
	{
		if(!_buffer_.empty())
		{
			_add_buffer_to_data_table();
		}
		return _data_table_;
	}

	// Column access
#if(1)
	column_type & col(const key_type & key)
	{
		return data_table().col(_key_map_.at(key));
	}
	const_column_type & col(const key_type & key) const
	{
		return data_table().col(_key_map_.at(key));
	}
#endif
};

}



#endif // _BRG_CONTAINER_LABELED_MATRIX_HPP_INCLUDED_
