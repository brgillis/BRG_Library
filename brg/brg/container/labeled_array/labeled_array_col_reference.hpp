/**********************************************************************\
 @file labeled_array_col_reference.hpp
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

#ifndef _BRG_CONTAINER_LABELED_ARRAY_COL_REFERENCE_HPP_INCLUDED_
#define _BRG_CONTAINER_LABELED_ARRAY_COL_REFERENCE_HPP_INCLUDED_

#include <utility>

#include <boost/type_traits/is_convertible.hpp>

#include "brg/container/labeled_array.hpp"

namespace brgastro {

template<typename labeled_array_type, typename col_type>
class labeled_array_col_reference
{

public:

	// Public typedefs

	typedef typename labeled_array_type::key_type key_type;
	typedef typename labeled_array_type::data_type value_type;
	typedef typename labeled_array_type::size_type size_type;

	typedef col_type col_type;
	typedef typename labeled_array_type::const_col_type const_col_type;

	typedef typename col_type::reference reference;
	typedef typename col_type::const_reference const_reference;

	typedef typename col_type::iterator iterator;
	typedef typename col_type::const_iterator const_iterator;
	typedef typename col_type::reverse_iterator reverse_iterator;
	typedef typename col_type::const_reverse_iterator const_reverse_iterator;

	typedef typename col_type::difference_type difference_type;

private:

	// Private Members
	key_type * _key_;
	col_type * _col_;

public:

	/// Constructor. Requires a pointer to a labeled_array's key and column
	labeled_array_col_reference(key_type * key, col_type * col)
	: _key_(key),
	  _col_(col)
	{
	}

	/// Virtual destructor
	virtual ~labeled_array_col_reference() {}

	// Iterator methods
#if(1)
	/// begin
	const_iterator begin() const noexcept
	{
		return _col_->begin();
	}
	/// begin
	iterator begin() noexcept
	{
		return _col_->begin();
	}
	/// end
	const_iterator end() const noexcept
	{
		return _col_->end();
	}
	/// end
	iterator end() noexcept
	{
		return _col_->end();
	}
	/// rbegin
	const_reverse_iterator rbegin() const noexcept
	{
		return _col_->rbegin();
	}
	/// rbegin
	reverse_iterator rbegin() noexcept
	{
		return _col_->rbegin();
	}
	/// rend
	const_reverse_iterator rend() const noexcept
	{
		return _col_->rend();
	}
	/// rend
	reverse_iterator rend() noexcept
	{
		return _col_->rend();
	}
	/// cbegin
	const_iterator cbegin() const noexcept
	{
		return _col_->cbegin();
	}
	/// cend
	const_iterator cend() const noexcept
	{
		return _col_->cend();
	}
	/// crbegin
	const_reverse_iterator crbegin() const noexcept
	{
		return _col_->crbegin();
	}
	/// crend
	const_reverse_iterator crend() const noexcept
	{
		return _col_->crend();
	}

#endif // Iterator methods

	// Size methods
#if(1)

	/// size
	size_type size() const noexcept
	{
		return _col_->size();
	}

	/// empty
	bool empty() const noexcept
	{
		return _col_->empty();
	}

#endif // Capacity methods

	// Element access
#if(1)

	/// Element access
	const_reference operator[] (const size_type & n) const
	{
		return _col_[n];
	}

	/// Element access
	reference operator[] (const size_type & n)
	{
		return _col_[n];
	}

	/// Range-checked element access
	const_reference at( const size_type & n ) const
	{
		return _col_->at(n);
	}

	/// Range-checked element access
	reference at( const size_type & n )
	{
		return _col_->at(n);
	}

	/// Access first element
	const_reference front() const
	{
		return _col_->front();
	}

	/// Access first element
	reference front()
	{
		return _col_->front();
	}

	/// Access last element
	const_reference back() const
	{
		return _col_->back();
	}

	/// Access last element
	reference back()
	{
		return _col_->back();
	}

	/// Access data
	const value_type* data() const noexcept
	{
		return _col_->data();
	}

	/// Access data
	value_type* data() noexcept
	{
		return _col_->data();
	}

#endif // Element access

	// Label access
#if(1)
	/// Get the label
	const key_type & label() const noexcept
	{
		return *_key_;
	}
#endif

	// Casting
#if(1)

	/// Cast to col_type
	operator col_type &() const
	{
		return *_col_;
	}

	/// Cast non-const version to const version
	template <typename other_col_type,
	typename std::enable_if<boost::is_convertible<other_col_type,col_type>, other_col_type>::type* = nullptr>
	labeled_array_col_reference( const labeled_array_col_reference<labeled_array_type,other_col_type> & other)
	: _key_(other._key_), _col_(other._col_) {}

#endif

};

}



#endif // _BRG_CONTAINER_LABELED_ARRAY_COL_REFERENCE_HPP_INCLUDED_
