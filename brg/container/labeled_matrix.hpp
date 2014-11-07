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

#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include <Eigen/Array>

#include "brg/container/const_labeled_matrix_row_reference.hpp"
#include "brg/container/insertion_ordered_map.hpp"
#include "brg/container/labeled_matrix_row_reference.hpp"
#include "brg/container/table_typedefs.hpp"

namespace brgastro {

template<typename value_type, typename key_type=std::string>
class labeled_matrix
{
public:

	// Public typedefs
	typedef typename Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic> data_table_type;
	typedef value_type value_type;
	//typedef data_type::Index size_type;
	typedef size_t size_type;
	typedef typename data_table_type::ColXpr column_type;
	typedef typename data_table_type::ConstColXpr const_column_type;
	typedef typename data_table_type::RowXpr row_type;
	typedef typename labeled_matrix_row_reference<data_table_type,key_type> row_reference_type;
	typedef typename data_table_type::ConstRowXpr const_row_type;
	typedef typename const_labeled_matrix_row_reference<data_table_type,key_type> const_row_reference_type;

private:

	// Private typedefs
	typedef typename std::unordered_map<key_type,size_type> map_type;
	typedef typename brgastro::table_map_t<data_table_type,key_type> buffer_type;

	mutable data_table_type _data_table_;
	mutable map_type _key_map_;
	mutable buffer_type _buffer_;

	friend class row_reference_type;
	friend class const_row_reference_type;

	void _add_buffer_to_data_table() const
	{
		// Check that we actually do have columns to add
		if(_buffer_.empty()) return;

		// Get the number of columns and rows
		const size_type num_cols_to_add = _buffer_.size();
		const size_type old_num_cols = _data_table_.cols();
		const size_type new_num_cols = _data_table_.cols() + num_cols_to_add;
		size_type num_rows;

		// If the data table is initially empty, we'll get the number of rows from the first
		// buffered column
		if(_data_table_.empty())
		{
			num_rows = _buffer_.begin()->second.size();
		}
		else // Otherwise use the current number of rows
		{
			num_rows = _data_table_.rows();
		}

		// Check that all columns we'll be adding are the proper size before we change anything

		for( const auto & column : _buffer_)
		{
			// Check that the column is the proper size
			if(column.second.size()!=num_rows)
			{
				throw std::logic_error("All columns added to a labeled_matrix must have the same size.\n");
			}
		}

		// Resize the data table, conserving existing elements
		_data_table_.conservativeResize(num_rows, new_num_cols);

		// Go through the buffer and add each column to the data table. If the key already exists,
		// replace that column instead
		size_type current_column = old_num_cols;
		size_type num_preexisting_cols = 0; // number of keys which we find already exist
		for( const auto & column : _buffer_)
		{
			key_type & key = column.first;
			auto & column_data = column.second;

			// Check if the key already exists
			if(_key_map_.count(key)==1)
			{
				// Instead of adding a new column, we'll replace the existing one
				size_type & column_index = _key_map_.at(key);
				_data_table_.col(column_index) = std::move(column_data);
				++num_preexisting_cols;
			}
			else // New column to add
			{
				// Add this to the key map
				_key_map_[key] = current_column;
				_data_table_.col(current_column) = std::move(column_data);

				// Increment the current column index
				++current_column;
			}
		}

		// If we overwrote any preexisting columns, resize the data table now to the proper size
		if(num_preexisting_cols>0)
		{
			_data_table_.conservativeResize(data_table_type::NoChange_t, new_num_cols-num_preexisting_cols);
		}

		// Clear the buffer, and we're done
		_buffer_.clear();

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
#if(1)
	data_table_type & data_table()
	{
		if(!_buffer_.empty())
		{
			_add_buffer_to_data_table();
		}
		return _data_table_;
	}
	const data_table_type & data_table() const
	{
		if(!_buffer_.empty())
		{
			_add_buffer_to_data_table();
		}
		return _data_table_;
	}
#endif

	// Column access
#if(1)
	column_type & col(const size_type & index)
	{
		return data_table().col(index);
	}
	const_column_type & col(const size_type & index) const
	{
		return data_table().col(index);
	}
	column_type & col_at_key(const key_type & key)
	{
		return data_table().col(_key_map_.at(key));
	}
	const_column_type & col_at_key(const key_type & key) const
	{
		return data_table().col(_key_map_.at(key));
	}
#endif

	// Row access
#if(1)
	row_reference_type & row(const size_type & index)
	{
		return row_reference_type(_key_map_,data_table().row(index));
	}
	const_row_reference_type & row(const size_type & index) const
	{
		return const_row_reference_type(_key_map_,data_table().row(index));
	}
#endif
};

}



#endif // _BRG_CONTAINER_LABELED_MATRIX_HPP_INCLUDED_
