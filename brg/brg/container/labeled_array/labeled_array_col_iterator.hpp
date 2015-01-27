/**********************************************************************\
 @file labeled_array_col_iterator.hpp
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

#ifndef _BRG_CONTAINER_LABELED_ARRAY_LABELED_ARRAY_COL_ITERATOR_HPP_INCLUDED_
#define _BRG_CONTAINER_LABELED_ARRAY_LABELED_ARRAY_COL_ITERATOR_HPP_INCLUDED_

#include <boost/iterator/iterator_facade.hpp>

namespace brgastro {

template< typename labeled_array_type, typename T, typename T_reference>
class labeled_array_col_iterator :
	public boost::iterator_facade<
		labeled_array_col_iterator<labeled_array_type,T,T_reference>, // CRTP
		T, // Value type
		boost::random_access_traversal_tag, // Traversal type
		T_reference, // Reference type
		size_t> // Difference type
{
private:
	// Private typedefs
	typedef typename boost::iterator_facade<labeled_array_col_iterator<labeled_array_type,T,T_reference>,
		T, boost::random_access_traversal_tag, T_reference, size_t> base;
public:
	// Public typedefs
	typedef typename base::difference_type difference_type;

private:
	labeled_array_type * _base_;
	size_t _col_index_;

	friend class boost::iterator_core_access;

	// Implementation of necessary operations
#if(1)

    void increment() { ++_col_index_; }
    void decrement() { --_col_index_; }
    void advance(const difference_type & n) { _col_index_ += n; }

    bool equal(const labeled_array_col_iterator<labeled_array_type,T,T_reference> & other) const
    {
        return (this->_base_ == other._base_) &&
        	(this->_col_index_ == other._col_index_);
    }

    difference_type distance_to(
    	const labeled_array_col_iterator<labeled_array_type,T,T_reference> & other) const
    {
    	return other._col_index_-_col_index_;
    }

    T_reference dereference() const { return _base_->col(_col_index_); }

#endif

public:
	labeled_array_col_iterator()
	: _base_(), _col_index_()
	{
	}

	labeled_array_col_iterator(labeled_array_type * base, const size_t & col_index)
	: _base_(base), _col_index_(col_index)
	{
	}

	template <class T_o, typename T_o_reference>
	labeled_array_col_iterator(
		const labeled_array_raw_col_iterator<labeled_array_type,T_o,T_o_reference> & other,
		typename std::enable_if<std::is_convertible<T_o_reference,T_reference>::value, T_o_reference *>::type = nullptr
	)
	: _base_(other._base_), _col_index_(other._col_index_) {}

};

}


#endif // _BRG_CONTAINER_LABELED_MATRIX_ARRAY_ARRAY_COL_ITERATOR_HPP_INCLUDED_
