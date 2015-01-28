/**********************************************************************\
 @file labeled_array.hpp
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

#ifndef _BRG_CONTAINER_LABELED_ARRAY_HPP_INCLUDED_
#define _BRG_CONTAINER_LABELED_ARRAY_HPP_INCLUDED_

#include <stdexcept>
#include <string>
#include <sstream>
#include <utility>
#include <type_traits>

#include <boost/bimap.hpp>
#include <boost/iterator.hpp>

#include <Eigen/Core>

#include "brg/container/coerce.hpp"
#include "brg/container/insertion_ordered_map.hpp"
#include "brg/container/is_container.hpp"
#include "brg/container/table_typedefs.hpp"

#include "brg/math/misc_math.hpp"

#include "labeled_array/labeled_array_col_iterator.hpp"
#include "labeled_array/labeled_array_col_reference.hpp"
#include "labeled_array/labeled_array_element_iterator.hpp"
#include "labeled_array/labeled_array_iterator_chooser.hpp"
#include "labeled_array/labeled_array_raw_col_iterator.hpp"
#include "labeled_array/labeled_array_raw_row_iterator.hpp"
#include "labeled_array/labeled_array_row_iterator.hpp"
#include "labeled_array/labeled_array_row_reference.hpp"
#include "labeled_array/labeled_array_vecs.hpp"

namespace brgastro {

template<typename T_value_type=double, char T_major_tag = Eigen::RowMajor, typename T_key_type=std::string>
class labeled_array
{
public:

	// Public typedefs
	typedef typename Eigen::Array<T_value_type,Eigen::Dynamic,Eigen::Dynamic,T_major_tag> data_table_type;
	typedef typename Eigen::Array<const T_key_type,Eigen::Dynamic,Eigen::Dynamic,T_major_tag> const_data_table_type;

	typedef T_value_type value_type;
	typedef const value_type const_value_type;
	typedef value_type & reference;
	typedef const value_type & const_reference;

	typedef T_key_type key_type;
	typedef const key_type const_key_type;

	typedef decltype(data_table_type().size()) size_type;
	typedef ptrdiff_t difference_type;

	// Col typedefs
#if(1)

	typedef typename data_table_type::ColXpr col_type;
	typedef typename data_table_type::ConstColXpr const_col_type;

	typedef labeled_array_col_reference<labeled_array<value_type,T_major_tag,key_type>,col_type> col_reference;
	typedef labeled_array_col_reference<labeled_array<value_type,T_major_tag,key_type>,const_col_type> const_col_reference;

	typedef labeled_array_col_iterator<labeled_array<value_type,T_major_tag,key_type>,
		col_type,col_reference> col_iterator;
	typedef labeled_array_col_iterator<labeled_array<value_type,T_major_tag,key_type>,
		const_col_type,const_col_reference> const_col_iterator;
	typedef typename boost::reverse_iterator<col_iterator> reverse_col_iterator;
	typedef typename boost::reverse_iterator<const_col_iterator> const_reverse_col_iterator;

	typedef labeled_array_vecs<labeled_array<value_type,T_major_tag,key_type>,
		col_type,const_col_type,col_reference,const_col_reference,col_iterator,const_col_iterator> cols_type;
	typedef labeled_array_vecs<labeled_array<value_type,T_major_tag,key_type>,
		const_col_type,const_col_type,const_col_reference,const_col_reference,const_col_iterator,const_col_iterator> const_cols_type;

	typedef labeled_array_raw_col_iterator<labeled_array<value_type,T_major_tag,key_type>,
		col_type,col_reference> raw_col_iterator;
	typedef labeled_array_raw_col_iterator<labeled_array<value_type,T_major_tag,key_type>,
		const_col_type,const_col_reference> const_raw_col_iterator;
	typedef typename boost::reverse_iterator<raw_col_iterator> reverse_raw_col_iterator;
	typedef typename boost::reverse_iterator<const_raw_col_iterator> const_reverse_raw_col_iterator;

	typedef labeled_array_vecs<labeled_array<value_type,T_major_tag,key_type>,
		col_type,const_col_type,col_reference,const_col_reference,
		raw_col_iterator,const_raw_col_iterator> raw_cols_type;
	typedef labeled_array_vecs<labeled_array<value_type,T_major_tag,key_type>,
		const_col_type,const_col_type,const_col_reference,const_col_reference,
		const_raw_col_iterator,const_raw_col_iterator> const_raw_cols_type;

#endif

	// Row typedefs
#if(1)

	typedef typename data_table_type::RowXpr row_type;
	typedef typename data_table_type::ConstRowXpr const_row_type;

	typedef labeled_array_row_reference<labeled_array<value_type,T_major_tag,key_type>,row_type> row_reference;
	typedef labeled_array_row_reference<labeled_array<value_type,T_major_tag,key_type>,const_row_type> const_row_reference;

	typedef labeled_array_row_iterator<labeled_array<value_type,T_major_tag,key_type>,
		row_type,row_reference> row_iterator;
	typedef labeled_array_row_iterator<labeled_array<value_type,T_major_tag,key_type>,
		const_row_type,const_row_reference> const_row_iterator;
	typedef typename boost::reverse_iterator<row_iterator> reverse_row_iterator;
	typedef typename boost::reverse_iterator<const_row_iterator> const_reverse_row_iterator;

	typedef labeled_array_vecs<labeled_array<value_type,T_major_tag,key_type>,
		row_type,const_row_type,row_reference,const_row_reference,row_iterator,const_row_iterator> rows_type;
	typedef labeled_array_vecs<labeled_array<value_type,T_major_tag,key_type>,
		const_row_type,const_row_type,const_row_reference,const_row_reference,const_row_iterator,const_row_iterator> const_rows_type;

	typedef labeled_array_raw_row_iterator<labeled_array<value_type,T_major_tag,key_type>,
		row_type,row_reference> raw_row_iterator;
	typedef labeled_array_raw_row_iterator<labeled_array<value_type,T_major_tag,key_type>,
		const_row_type,const_row_reference> const_raw_row_iterator;
	typedef typename boost::reverse_iterator<raw_row_iterator> reverse_raw_row_iterator;
	typedef typename boost::reverse_iterator<const_raw_row_iterator> const_reverse_raw_row_iterator;

	typedef labeled_array_vecs<labeled_array<value_type,T_major_tag,key_type>,
		row_type,const_row_type,row_reference,const_row_reference,
		raw_row_iterator,const_raw_row_iterator> raw_rows_type;
	typedef labeled_array_vecs<labeled_array<value_type,T_major_tag,key_type>,
		const_row_type,const_row_type,const_row_reference,const_row_reference,
		const_raw_row_iterator,const_raw_row_iterator> const_raw_rows_type;

#endif

	// Iterators
#if(1)

	typedef value_type * iterator;
	typedef const_value_type * const_iterator;
	typedef boost::reverse_iterator<iterator> reverse_iterator;
	typedef boost::reverse_iterator<const_iterator> const_reverse_iterator;

	typedef typename labeled_array_iterator_chooser<value_type,T_major_tag>::col_element_iterator col_element_iterator;
	typedef typename labeled_array_iterator_chooser<const_value_type,T_major_tag>::col_element_iterator const_col_element_iterator;
	typedef typename boost::reverse_iterator<col_element_iterator> reverse_col_element_iterator;
	typedef typename boost::reverse_iterator<const_col_element_iterator> const_reverse_col_element_iterator;

	typedef typename labeled_array_iterator_chooser<value_type,T_major_tag>::row_element_iterator row_element_iterator;
	typedef typename labeled_array_iterator_chooser<const_value_type,T_major_tag>::row_element_iterator const_row_element_iterator;
	typedef typename boost::reverse_iterator<row_element_iterator> reverse_row_element_iterator;
	typedef typename boost::reverse_iterator<const_row_element_iterator> const_reverse_row_element_iterator;

#endif

private:

#if(1)
	// Private typedefs
	typedef typename boost::bimap<key_type,size_type> map_type;

	typedef typename Eigen::Array<value_type,Eigen::Dynamic,1> column_buffer_column_type;
	typedef typename brgastro::insertion_ordered_map<key_type,column_buffer_column_type> column_buffer_type;
	typedef typename column_buffer_type::value_type column_buffer_labeled_column_type;

	typedef typename std::vector<value_type> row_buffer_row_type;
	typedef typename std::vector<row_buffer_row_type> row_buffer_type;
	typedef typename row_buffer_type::value_type row_buffer_labeled_row_type;

	mutable data_table_type _data_table_;
	mutable map_type _key_map_;
	mutable column_buffer_type _column_buffer_;
	mutable row_buffer_type _row_buffer_;

	friend row_reference;
	friend const_row_reference;
	friend cols_type;
	friend const_cols_type;
	friend rows_type;
	friend const_rows_type;

	void _add_buffer_to_data_table() const
	{
		_add_column_buffer_to_data_table();
		_add_row_buffer_to_data_table();
	}

	void _add_column_buffer_to_data_table() const
	{
		// Check that we actually do have columns to add
		if(_column_buffer_.empty()) return;

		// Get the number of columns
		const size_type num_cols_to_add = _column_buffer_.size();
		const size_type old_num_cols = _data_table_.cols();
		const size_type new_num_cols = _data_table_.cols() + num_cols_to_add;
		size_type num_rows;

		// If the data table is initially empty, we'll get the number of rows from the first
		// buffered column
		if(_data_table_.size()==0)
		{
			num_rows = _column_buffer_.begin()->second.size();
		}
		else // Otherwise use the current number of rows
		{
			num_rows = _data_table_.rows();
		}

		// Check that all columns we'll be adding are the proper size before we change anything

		for( const auto & column : _column_buffer_)
		{
			// Check that the column is the proper size
			if(column.second.size()!=num_rows)
			{
				throw std::logic_error("All columns added to a labeled_array must have the same size.\n");
			}
		}

		// Resize the data table, conserving existing elements
		_data_table_.conservativeResize(num_rows, new_num_cols);

		// Go through the buffer and add each column to the data table. If the key already exists,
		// replace that column instead
		size_type current_column = old_num_cols;
		size_type num_preexisting_cols = 0; // number of keys which we find already exist
		for( const auto & column : _column_buffer_)
		{
			const_key_type & key = column.first;
			auto & column_data = column.second;

			// Check if the key already exists
			if(_key_map_.left.count(key)==1)
			{
				// Instead of adding a new column, we'll replace the existing one
				const size_type & column_index = _key_map_.left.at(key);
				for(size_type i=0; i<column_data.size(); ++i)
				{
					_data_table_(i,column_index) = column_data[i];
				}
				++num_preexisting_cols;
			}
			else // New column to add
			{
				// Add this to the key map
				_key_map_.left.insert(std::make_pair(key,current_column));
				for(size_type i=0; i<column_data.size(); ++i)
				{
					_data_table_(i,current_column) = column_data[i];
				}

				// Increment the current column index
				++current_column;
			}
		}

		// If we overwrote any preexisting columns, resize the data table now to the proper size
		if(num_preexisting_cols>0)
		{
			_data_table_.conservativeResize(Eigen::NoChange_t(), new_num_cols-num_preexisting_cols);
		}

		// Clear the buffer, and we're done
		_column_buffer_.clear();

	}

	void _add_row_buffer_to_data_table() const
	{
		// Check that we actually do have rows to add
		if(_row_buffer_.empty()) return;

		// Get the number of rows
		const size_type num_rows_to_add = _row_buffer_.size();
		const size_type old_num_rows = _data_table_.rows();
		const size_type new_num_rows = _data_table_.rows() + num_rows_to_add;
		size_type num_cols;

		// If the data table is initially empty, we'll get the number of columns from the first
		// buffered row
		if(_data_table_.size()==0)
		{
			num_cols = _row_buffer_.front().size();
		}
		else // Otherwise use the current number of columns
		{
			num_cols = _data_table_.cols();
		}

		// Check that all rows we'll be adding are the proper size before we change anything

		for( const auto & row : _row_buffer_)
		{
			// Check that the row is the proper size
			if(static_cast<size_type>(row.size())!=num_cols)
			{
				throw std::logic_error("All rows added to a labeled_array must have the same size.\n");
			}
		}

		// If this is the first row being added, create a dummy key map
		if(_data_table_.size()==0)
		{
			_generate_dummy_key_map(num_cols);
			assert(static_cast<size_type>(_key_map_.left.size())==num_cols);
		}

		// Resize the data table, conserving existing elements
		_data_table_.conservativeResize(new_num_rows, num_cols);

		size_type current_row = old_num_rows;
		// Go through the buffer and add each row to the data table
		for( const auto & row : _row_buffer_)
		{
			for(typename row_buffer_row_type::size_type i=0; i<row.size(); ++i)
			{
				_data_table_(current_row,i) = row[i];
			}
			++current_row;
		}

		// Clear the buffer, and we're done
		_row_buffer_.clear();

	}

	void _generate_dummy_key_map(const size_type & num_cols) const
	{
		_key_map_.left.clear();
		for(size_type i=0; i<num_cols; ++i)
		{
			std::stringstream ss("col_");
			ss << i;
			_key_map_.left.insert(typename map_type::left_value_type(ss.str(),i));
		}
	}

	// Data access after adding in buffer
	data_table_type & _data_table() const
	{
		_add_buffer_to_data_table();
		return _data_table_;
	}
#endif

public:

	// Constructors
#if(1)

	/// Default constructor
	labeled_array() {}

	/// Copy/move from table_map
	template< typename other_table_map_type,
	typename std::enable_if<brgastro::is_const_container<other_table_map_type>::value,other_table_map_type>::type* = nullptr,
	typename std::enable_if<brgastro::is_const_container<typename other_table_map_type::mapped_type>::value,other_table_map_type>::type* = nullptr>
	explicit labeled_array(other_table_map_type && init_column_buffer)
	: _column_buffer_(std::forward<other_table_map_type>(init_column_buffer))
	{
	}

	/// Copy/move from vector of vectors
	template< typename other_table_map_type,
	typename std::enable_if<brgastro::is_const_container<other_table_map_type>::value,other_table_map_type>::type* = nullptr,
	typename std::enable_if<brgastro::is_const_container<typename other_table_map_type::value_type>::value,other_table_map_type>::type* = nullptr>
	explicit labeled_array(other_table_map_type && init_row_buffer)
	: _row_buffer_(std::forward<other_table_map_type>(init_row_buffer))
	{
	}

#endif

	// Data access
#if(1)
	const data_table_type & data_table() const
	{
		_add_buffer_to_data_table();
		return _data_table_;
	}
#endif

	// Column access
#if(1)
	col_reference col(const size_type & index)
	{
		return col_reference(&(get_key_for_column(index)),_data_table().col(index),num_cols());
	}
	const_col_reference col(const size_type & index) const
	{
		return const_col_reference(&(get_key_for_column(index)),_data_table().col(index),num_cols());
	}
	col_type raw_col(const size_type & index)
	{
		return _data_table().col(index);
	}
	const_col_type raw_col(const size_type & index) const
	{
		return _data_table().col(index);
	}
	col_type at(const size_type & index)
	{
		return col(index);
	}
	const_col_type at(const size_type & index) const
	{
		return col(index);
	}
	col_type at_label(const key_type & key)
	{
		return _data_table().col(_key_map_.left.at(key));
	}
	const_col_type at_label(const key_type & key) const
	{
		return _data_table().col(_key_map_.left.at(key));
	}
#endif

	// Iterable table access
#if(1)

	// Get an iterable set of columns
	cols_type cols()
	{
		return cols_type(this,num_cols());
	}
	const_cols_type cols() const
	{
		return const_cols_type(this,num_cols());
	}
	const_cols_type const_cols() const
	{
		return const_cols_type(this,num_cols());
	}

	// Get an iterable set of raw columns
	raw_cols_type raw_cols()
	{
		return raw_cols_type(this,num_cols());
	}
	const_raw_cols_type raw_cols() const
	{
		return const_raw_cols_type(this,num_cols());
	}
	const_raw_cols_type const_raw_cols() const
	{
		return const_raw_cols_type(this,num_cols());
	}

	// Get an iterable set of rows
	rows_type rows()
	{
		return rows_type(this,num_rows());
	}
	const_rows_type rows() const
	{
		return const_rows_type(this,num_rows());
	}
	const_rows_type const_rows() const
	{
		return const_rows_type(this,num_rows());
	}

	// Get an iterable set of raw rows
	raw_rows_type raw_rows()
	{
		return raw_rows_type(this,num_rows());
	}
	const_raw_rows_type raw_rows() const
	{
		return const_raw_rows_type(this,num_rows());
	}
	const_raw_rows_type const_raw_rows() const
	{
		return const_raw_rows_type(this,num_rows());
	}

#endif

    // Column/row insertion
#if(1)
    template< typename new_column_type,
	typename std::enable_if<std::is_convertible<new_column_type,column_buffer_column_type>::value ||
	std::is_convertible<new_column_type,column_buffer_labeled_column_type>::value,char>::type = 0>
    void insert_col(new_column_type && new_column)
    {
    	_add_row_buffer_to_data_table();
    	_column_buffer_.insert(std::forward<new_column_type>(new_column));
    }

    template< typename new_column_type,
	typename std::enable_if<!(std::is_convertible<new_column_type,column_buffer_column_type>::value ||
	std::is_convertible<new_column_type,column_buffer_labeled_column_type>::value),char>::type = 0,
	typename std::enable_if<brgastro::is_const_container<typename std::decay<new_column_type>::type>::value,char>::type = 0>
    void insert_col(new_column_type && new_column)
    {
    	_add_row_buffer_to_data_table();
    	_column_buffer_.insert(brgastro::coerce<column_buffer_column_type>(std::forward<new_column_type>(new_column)));
    }

    template< typename new_row_type,
	typename std::enable_if<std::is_convertible<new_row_type,row_buffer_row_type>::value ||
	std::is_convertible<new_row_type,row_buffer_labeled_row_type>::value,char>::type = 0>
    void insert_row(new_row_type && new_row)
    {
    	_add_column_buffer_to_data_table();
    	_row_buffer_.push_back(std::forward<new_row_type>(new_row));
    }

    template< typename new_row_type,
	typename std::enable_if<!(std::is_convertible<new_row_type,row_buffer_row_type>::value ||
	std::is_convertible<new_row_type,row_buffer_row_type>::value),char>::type = 0,
	typename std::enable_if<brgastro::is_const_container<typename std::decay<new_row_type>::type>::value,char>::type = 0>
    void insert_row(new_row_type && new_column)
    {
    	_add_column_buffer_to_data_table();
    	_row_buffer_.push_back(brgastro::coerce<row_buffer_row_type>(std::forward<new_row_type>(new_column)));
    }
#endif

    // Column access/insertion
#if(1)
    col_type operator[](const key_type & key)
	{
    	// Check if the key is new
    	if(count(key)==0)
    	{
    		// The key is new, so insert a column for it
    		insert_col(std::make_pair(key,column_buffer_column_type(_data_table_.rows())));
    		_data_table();
    	}
    	// Return a reference to this key's column
		return _data_table().col(_key_map_.left.at(key));
	}
#endif

	// Row access
#if(1)
	row_reference row(const size_type & index)
	{
		return row_reference(&_key_map_,_data_table().row(index),num_rows());
	}
	const_row_reference row(const size_type & index) const
	{
		return const_row_reference(&_key_map_,_data_table().row(index),num_rows());
	}
	row_reference operator[](const size_type & index)
	{
		return row(index);
	}
	const_row_reference operator[](const size_type & index) const
	{
		return row(index);
	}
	row_reference operator()(const size_type & index)
	{
		return row(index);
	}
	const_row_reference operator()(const size_type & index) const
	{
		return row(index);
	}
	row_type raw_row(const size_type & index)
	{
		return _data_table().row(index);
	}
	const_row_type raw_row(const size_type & index) const
	{
		return _data_table().row(index);
	}
#endif

	// Element access
#if(1)
	const_reference operator()(const size_type & row_index, const size_type & col_index) const
	{
		return _data_table()(row_index,col_index);
	}
	reference operator()(const size_type & row_index, const size_type & col_index)
	{
		return _data_table()(row_index,col_index);
	}
#endif

	// Size information
#if(1)
	size_type nrows() const
	{
		return num_rows();
	}
	size_type num_rows() const
	{
		return _data_table().rows();
	}
	size_type ncols() const
	{
		return num_cols();
	}
	size_type num_cols() const
	{
		return _data_table().cols();
	}
	size_type size() const
	{
		return _data_table().size();
	}
#endif

	// Key information
#if(1)
	char count (const key_type & k) const
    {
    	_add_buffer_to_data_table();
		return _key_map_.left.count(k);
    }
#endif

    // Key control
#if(1)

    /**
     * Changes a key to another, without altering its mapped data or position.
     *
     * @param init_key
     * @param new_key
     * @return char - 0 if successful
     *                1 if init_key doesn't exist in map
     *                2 if new_key already exists in map
     */
    template< typename new_key_type >
    char change_key(const key_type & init_key, new_key_type && new_key)
    {
    	// Add in the buffer so the key map is fully set up
		_add_buffer_to_data_table();

    	// Get the position of the value we're going to be altering
    	auto it = _key_map_.left.find(init_key);

    	if(it==_key_map_.left.end()) return 1;
    	if(_key_map_.left.count(new_key)>0) return 2;

    	size_t pos = it->second;

    	// Alter the value in the key map by erasing old entry and adding new
    	_key_map_.left.erase(it);
    	_key_map_.left[std::forward<new_key_type>(new_key)] = pos;

    	return 0;
    }

    const key_type & get_key_for_column(const size_type & index) const
    {
    	return _key_map_.right.at(index);
    }

#endif

    // Advanced operations
#if(1)

    // Apply unit conversions
    template< typename unitconvs >
    void apply_unitconvs( const unitconvs & u_map )
    {
    	_add_buffer_to_data_table();
    	for(auto u_it=u_map.begin(); u_it!=u_map.end(); ++u_it)
    	{
    		auto d_it = _key_map_.left.find(u_it->first);
    		// Check if this value is actually in the data table
    		if(d_it!=_key_map_.left.end())
    		{
    			// It's in the table, so multiply its associated vector by the inverse unit conversion
    			double factor = 1./(u_it->second);
    			if(isbad(factor)) factor = 1; // To catch user mistakes
    			col(d_it->second) *= factor;
    		}
    	}
    }

#endif

};

}



#endif // _BRG_CONTAINER_LABELED_ARRAY_HPP_INCLUDED_
