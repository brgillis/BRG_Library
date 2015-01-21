/**********************************************************************\
 @file labeled_array_vecs.hpp
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

#ifndef _BRG_CONTAINER_LABELED_ARRAY_VECS_HPP_INCLUDED_
#define _BRG_CONTAINER_LABELED_ARRAY_VECS_HPP_INCLUDED_

#include <utility>

#include "brg/container/labeled_array.hpp"

namespace brgastro {

template<typename labeled_array_type,
typename vec_type, typename const_vec_type,
typename reference, typename const_reference,
typename iterator, typename const_iterator>
class labeled_array_vecs
{
private:

	// Private typedefs
	typedef typename labeled_array_type::map_type map_type;

	// Members
	labeled_array_type * const _array_;
	typename labeled_array_type::size_type _size_;

public:

	// Public typedefs
	typedef typename labeled_array_type::key_type key_type;
	typedef typename labeled_array_type::value_type value_type;
	typedef typename labeled_array_type::size_type size_type;

	typedef vec_type vec_type;
	typedef const_vec_type const_vec_type;

	typedef reference reference;
	typedef const_reference const_reference;

	typedef iterator iterator;
	typedef const_iterator const_iterator;
	typedef iterator reverse_iterator; // FIXME
	typedef const_iterator const_reverse_iterator; // FIXME

	typedef typename labeled_array_type::difference_type difference_type;

	/// Constructor. Requires a non-const pointer to a labeled array
	labeled_array_vecs(labeled_array_type * init_array, const size_type & init_size)
	: _array_(init_array), _size_(init_size)
	{
	}

	/// Virtual destructor
	virtual ~labeled_array_vecs() {}

	// Iterator methods
#if(1)
	/// begin
	const_iterator begin() const noexcept
	{
		return const_iterator(_array_,0);
	}
	/// begin
	iterator begin() noexcept
	{
		return iterator(_array_,0);
	}
	/// end
	const_iterator end() const noexcept
	{
		return const_iterator(_array_,_size_);
	}
	/// end
	iterator end() noexcept
	{
		return iterator(_array_,_size_);
	}
	/// begin
	const_reverse_iterator rbegin() const noexcept
	{
		return const_reverse_iterator(_array_,0);
	}
	/// begin
	reverse_iterator rbegin() noexcept
	{
		return reverse_iterator(_array_,0);
	}
	/// end
	const_reverse_iterator rend() const noexcept
	{
		return const_reverse_iterator(_array_,_array_->num_rows());
	}
	/// end
	reverse_iterator rend() noexcept
	{
		return reverse_iterator(_array_,_array_->num_rows());
	}
	/// begin
	const_iterator cbegin() const noexcept
	{
		return const_iterator(_array_,0);
	}
	/// end
	const_iterator cend() const noexcept
	{
		return const_iterator(_array_,_array_->num_rows());
	}
	/// begin
	const_reverse_iterator crbegin() const noexcept
	{
		return const_reverse_iterator(_array_,0);
	}
	/// end
	const_reverse_iterator crend() const noexcept
	{
		return const_reverse_iterator(_array_,_array_->num_rows());
	}

#endif // Iterator methods

	// Size methods
#if(1)

	/// size
	size_type size() const noexcept
	{
		return _size_;
	}

	/// empty
	bool empty() const noexcept
	{
		return _size_==0;
	}

#endif // Capacity methods

	// Row access
#if(1)

	/// Element access
	const_reference operator[] (const size_type & n) const
	{
		return const_reference(_array_->_key_map_,n);
	}

	/// Element access
	reference operator[] (const size_type & n)
	{
		return reference(_array_->_key_map_,n);
	}

	/// Range-checked element access
	const_reference at( const size_type & n ) const
	{
		if((n<0)||n>_size_) throw std::out_of_range;
		return const_reference(_array_->_key_map_,n);
	}

	/// Range-checked element access
	reference at( const size_type & n )
	{
		if((n<0)||n>_size_) throw std::out_of_range;
		return reference(_array_->_key_map_,n);
	}

	/// Access first element
	const_reference front() const
	{
		return const_reference(_array_->_key_map_,0);
	}

	/// Access first element
	reference front()
	{
		return reference(_array_->_key_map_,0);
	}

	/// Access last element
	const_reference back() const
	{
		return const_reference(_array_->_key_map_,_size_-1);
	}

	/// Access last element
	reference back()
	{
		return reference(_array_->_key_map_,_size_-1);
	}

	/// Access data
	const value_type* data() const noexcept
	{
		return _array_->data();
	}

	/// Access data
	value_type* data() noexcept
	{
		return _array_->data();
	}

#endif // Row access

	// Casting
#if (1)

	/// Cast non-const version to const version
	template <typename other_vec_type, typename other_reference, typename other_iterator,
	typename std::enable_if<boost::is_convertible<other_vec_type,vec_type>, other_vec_type>::type* = nullptr,
	typename std::enable_if<boost::is_convertible<other_reference,reference>, other_reference>::type* = nullptr,
	typename std::enable_if<boost::is_convertible<other_iterator,iterator>, other_iterator>::type* = nullptr>
	labeled_array_vecs( const labeled_array_vecs<labeled_array_type,other_vec_type,const_vec_type,other_reference,const_reference,
						other_iterator,const_iterator> & other)
	: _array_(other._array_) {}

#endif

};

}



#endif // _BRG_CONTAINER_LABELED_ARRAY_VECS_HPP_INCLUDED_
