/**********************************************************************\
 @file labeled_array_row_reference.hpp
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

#ifndef _BRG_CONTAINER_LABELED_ARRAY_ROW_REFERENCE_HPP_INCLUDED_
#define _BRG_CONTAINER_LABELED_ARRAY_ROW_REFERENCE_HPP_INCLUDED_

#include <cstdlib>
#include <utility>

#include <boost/type_traits/is_convertible.hpp>

#include "brg/container/labeled_array.hpp"

namespace brgastro {

template<typename labeled_array_type, typename T_row_type>
class labeled_array_row_reference
{
public:

	typedef typename labeled_array_type::key_type key_type;
	typedef typename labeled_array_type::value_type value_type;
	typedef typename labeled_array_type::size_type size_type;

	typedef T_row_type row_type;
	typedef typename labeled_array_type::const_row_type const_row_type;

	typedef typename labeled_array_type::reference reference;
	typedef typename labeled_array_type::const_reference const_reference;

	typedef typename labeled_array_type::row_element_iterator iterator;
	typedef typename labeled_array_type::const_row_element_iterator const_iterator;
	typedef typename labeled_array_type::reverse_row_element_iterator reverse_iterator;
	typedef typename labeled_array_type::const_reverse_row_element_iterator const_reverse_iterator;

	typedef size_t difference_type;

private:

	// Private typedefs
	typedef typename labeled_array_type::map_type map_type;

	// Members
	map_type * _key_map_;
	row_type _row_;
	size_type _num_cols_;

public:

	/// Constructor. Requires a pointer to a labeled_array's key map and row
	template <typename T_init_row_type>
	labeled_array_row_reference(map_type * key_map, T_init_row_type && row, const size_type & num_cols)
	: _key_map_(key_map),
	  _row_(std::forward<T_init_row_type>(row)),
	  _num_cols_(num_cols)
	{
	}

	/// Virtual destructor
	virtual ~labeled_array_row_reference() {}

	// Iterator methods
#if(1)
	/// begin
	const_iterator begin() const noexcept
	{
		return const_iterator(_row_.data(),_num_cols_);
	}
	/// begin
	iterator begin() noexcept
	{
		return iterator(_row_.data(),_num_cols_);
	}
	/// end
	const_iterator end() const noexcept
	{
		return const_iterator(_row_.data()+_row_.size()*_num_cols_,_num_cols_);
	}
	/// end
	iterator end() noexcept
	{
		return iterator(_row_.data()+_row_.size()*_num_cols_,_num_cols_);
	}
	/// rbegin
	const_reverse_iterator rbegin() const noexcept
	{
		return const_reverse_iterator(_row_.data()+_row_.size()*_num_cols_-_num_cols_,_num_cols_);
	}
	/// rbegin
	reverse_iterator rbegin() noexcept
	{
		return reverse_iterator(_row_.data()+_row_.size()*_num_cols_-_num_cols_,_num_cols_);
	}
	/// rend
	const_reverse_iterator rend() const noexcept
	{
		return const_reverse_iterator(_row_.data()-_num_cols_,_num_cols_);
	}
	/// rend
	reverse_iterator rend() noexcept
	{
		return reverse_iterator(_row_.data()-_num_cols_,_num_cols_);
	}
	/// cbegin
	const_iterator cbegin() const noexcept
	{
		return const_iterator(_row_.data(),_num_cols_);
	}
	/// cend
	const_iterator cend() const noexcept
	{
		return const_iterator(_row_.data()+_row_.size()*_num_cols_,_num_cols_);
	}
	/// crbegin
	const_reverse_iterator crbegin() const noexcept
	{
		return const_reverse_iterator(_row_.data()+_row_.size()*_num_cols_-_num_cols_,_num_cols_);
	}
	/// crend
	const_reverse_iterator crend() const noexcept
	{
		return const_reverse_iterator(_row_.data()-_num_cols_,_num_cols_);
	}

#endif // Iterator methods

	// Size methods
#if(1)

	/// size
	size_type size() const noexcept
	{
		return _row_.size();
	}

	/// empty
	bool empty() const noexcept
	{
		return _row_.empty();
	}

#endif // Capacity methods

	// Element access
#if(1)

	/// Element access
	const_reference operator[] (const size_type & n) const
	{
		return _row_[n];
	}

	/// Element access
	reference operator[] (const size_type & n)
	{
		return _row_[n];
	}

	/// Range-checked element access
	const_reference operator()( const size_type & n ) const
	{
		return _row_(n);
	}

	/// Range-checked element access
	reference operator()( const size_type & n )
	{
		return _row_(n);
	}

	/// Range-checked element access
	const_reference at( const size_type & n ) const
	{
		return _row_(n);
	}

	/// Range-checked element access
	reference at( const size_type & n )
	{
		return _row_(n);
	}

	/// Range-checked element access
	const_reference at_label( const key_type & key ) const
	{
		return _row_(_key_map_->left.at(key));
	}

	/// Range-checked element access
	reference at_label( const key_type & key )
	{
		return _row_(_key_map_->left.at(key));
	}

	/// Access first element
	const_reference front() const
	{
		return _row_(0);
	}

	/// Access first element
	reference front()
	{
		return _row_(0);
	}

	/// Access last element
	const_reference back() const
	{
		return _row_(_row_.size()-1);
	}

	/// Access last element
	reference back()
	{
		return _row_(_row_.size()-1);
	}

	/// Access data
	const value_type* data() const noexcept
	{
		return _row_.data();
	}

	/// Access data
	value_type* data() noexcept
	{
		return _row_.data();
	}

#endif // Element access

	// Label access
#if(1)

	key_type label(const size_type & n) const
	{
		return _key_map_->right.at(n);
	}

#endif

	// Casting
#if(1)

	/// Cast to row_type
	operator row_type &() const
	{
		return _row_;
	}

	/// Cast non-const version to const version
	template <typename other_row_type,
	typename std::enable_if<std::is_convertible<other_row_type,row_type>::value, other_row_type>::type* = nullptr>
	labeled_array_row_reference( const labeled_array_row_reference<labeled_array_type,other_row_type> & other)
	: _key_map_(other._key_map_), _row_(other._row_), _num_cols_(other._num_cols_) {}

#endif

};

}



#endif // _BRG_CONTAINER_LABELED_ARRAY_ROW_REFERENCE_HPP_INCLUDED_
