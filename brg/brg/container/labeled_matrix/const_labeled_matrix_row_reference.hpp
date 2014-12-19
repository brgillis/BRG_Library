/**********************************************************************\
 @file const_labeled_matrix_row_reference.hpp
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

#ifndef _BRG_CONTAINER_CONST_LABELED_MATRIX_ROW_REFERENCE_HPP_INCLUDED_
#define _BRG_CONTAINER_CONST_LABELED_MATRIX_ROW_REFERENCE_HPP_INCLUDED_

#include <utility>

#include "brg/container/labeled_matrix.hpp"

namespace brgastro {

template<typename data_type, typename key_type=std::string>
class const_labeled_matrix_row_reference
{
private:

	// Private typedefs
	typedef labeled_matrix<data_type,key_type> labeled_matrix_type;
	typedef typename labeled_matrix_type::map_type map_type;
	typedef typename labeled_matrix_type::const_row_type row_type;

	// Members
	const map_type & _key_map_;
	row_type & _row_;

public:

	// Public typedefs
	typedef key_type key_type;
	typedef data_type value_type;
	typedef typename labeled_matrix_type::size_type size_type;
	typedef typename labeled_matrix_type::reference reference;
	typedef typename labeled_matrix_type::const_reference const_reference;
	typedef typename labeled_matrix_type::pointer pointer;
	typedef typename labeled_matrix_type::const_pointer const_pointer;
	typedef typename row_type::const_iterator const_iterator;
	typedef typename row_type::const_reverse_iterator const_reverse_iterator;
	typedef typename row_type::difference_type difference_type;

	/// Constructor. Requires a non-const reference to a labeled_matrix's key map and row
	const_labeled_matrix_row_reference(const map_type & key_map, row_type & row)
	: _key_map_(key_map),
	  _row_(row)
	{
	}

	/// Construct from non-const labeled matrix row reference
	const_labeled_matrix_row_reference(const labeled_matrix_row_reference<data_type,key_type> & other )
	: _key_map_(other._key_map_),
	  _row_(other._row_)
	{
	}


	/// Virtual destructor
	virtual ~const_labeled_matrix_row_reference() {}

	/// Delete default copy and move constructors and assignments

	/// Delete copy constructors and assignments
#if(1)
	const_labeled_matrix_row_reference(const const_labeled_matrix_row_reference & other) = delete;
	const_labeled_matrix_row_reference & operator=(const const_labeled_matrix_row_reference & other) = delete;
#endif

	/// Move constructors and assignments
#if(1)
	const_labeled_matrix_row_reference(const_labeled_matrix_row_reference && other)
	: _key_map_(std::move(other._key_map_)),
	  _row_(std::move(other._row_))
	{
	}
	const_labeled_matrix_row_reference & operator=(const_labeled_matrix_row_reference && other)
	{
		_key_map_ = std::move(other._key_map_);
		_row_ = std::move(other._row_);
		return *this;
	}
#endif

	// Iterator methods
#if(1)
	/// begin
	const_iterator begin() const noexcept
	{
		return _row_.begin();
	}
	/// end
	const_iterator end() const noexcept
	{
		return _row_.end();
	}
	/// rbegin
	const_iterator rbegin() const noexcept
	{
		return _row_.rbegin();
	}
	/// rend
	const_iterator rend() const noexcept
	{
		return _row_.rend();
	}
	/// cbegin
	const_iterator cbegin() const noexcept
	{
		return _row_.cbegin();
	}
	/// cend
	const_iterator cend() const noexcept
	{
		return _row_.cend();
	}
	/// crbegin
	const_iterator crbegin() const noexcept
	{
		return _row_.crbegin();
	}
	/// crend
	const_iterator crend() const noexcept
	{
		return _row_.crend();
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

	/// Range-checked element access
	const_reference at( const size_type & n ) const
	{
		return _row_.at(n);
	}

	/// Range-checked element access
	const_reference at_label( const key_type & key ) const
	{
		return _row_[_key_map_.at(key)];
	}

	/// Access first element
	const_reference front() const
	{
		return _row_.front();
	}

	/// Access last element
	const_reference back() const
	{
		return _row_.back();
	}

	/// Access data
	const value_type* data() const noexcept
	{
		return _row_.data();
	}

#endif // Element access

	// Casting
#if(1)

	/// Cast to row_type
	operator row_type &() const
	{
		return _row_;
	}

#endif

};

}



#endif // _BRG_CONTAINER_LABELED_MATRIX_ROW_REFERENCE_HPP_INCLUDED_
