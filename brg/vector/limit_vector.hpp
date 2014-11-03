/*
 * limit_vector.hpp
 *
 *  Created on: Nov 3, 2014
 *      Author: user
 */

#ifndef _BRG_VECTOR_LIMIT_VECTOR_HPP_
#define _BRG_VECTOR_LIMIT_VECTOR_HPP_

#include <cassert>
#include <limits>
#include <vector>

#include "brg/vector/summary_functions.hpp"

namespace brgastro {

template < class T, class A = std::allocator<T> >
class limit_vector: protected std::vector<T,A>
{
public:
	// Some typedefs that we can allow users to access
	using std::vector<T,A>::value_type;
	using std::vector<T,A>::allocator_type;
	using std::vector<T,A>::reference;
	using std::vector<T,A>::const_reference;
	using std::vector<T,A>::pointer;
	using std::vector<T,A>::const_pointer;
	using std::vector<T,A>::iterator;
	using std::vector<T,A>::const_iterator;
	using std::vector<T,A>::reverse_iterator;
	using std::vector<T,A>::const_reverse_iterator;
	using std::vector<T,A>::difference_type;
	using std::vector<T,A>::size_type;

	// Enum to describe type of limit vector this is
	enum class limit_type {GENERAL, LINEAR, LOG};

private:
	std::vector<T,A> base;

	limit_type type;

	// Private functions
#if(1)
	/// Init function - sets limits to cover full numeric range. Must not be called when base is not empty
	void _init()
	{
		assert(base.empty());

		type = limit_type::GENERAL;

		base.push_back(std::numeric_limits<T>::lowest());
		base.push_back(std::numeric_limits<T>::max());
	} // void _init()

	/// Check if is valid, and throw if it isn't
	void _check_if_valid_construction()
	{
		if(!is_monotonically_increasing(base))
			throw std::runtime_error("Limit vector cannot be constructed from vector which isn't monotonic increasing.\n");
	}

	/// Check if is valid, and throw if it isn't
	void _check_if_valid_setting(const std::vector<T,A> & other)
	{
		if(!is_monotonically_increasing(other))
			throw std::runtime_error("Limit vector cannot be set to vector which isn't monotonic increasing.\n");
	}
#endif // Private functions

public:

	/// Clear function - leaves it in initial state
	void clear()
	{
		base.clear();
		_init();
	}

	// Constructors
#if(1)

	/// Default constructor
	explicit limit_vector(const allocator_type& alloc = allocator_type())
	: base(alloc)
	{
		_init();
	}

	/// Construct as linear or log format

	/// Range constructor
	template <class InputIterator>
	limit_vector(InputIterator first, InputIterator last,
			const allocator_type& alloc = allocator_type())
	: base(first,last,alloc),
	  type(limit_type::GENERAL)
	{
		_check_if_valid_construction();
	}

	/// Copy constructor
	limit_vector(const limit_vector<T,A> & other) = default;

	/// Copy and set allocator
	limit_vector(const limit_vector<T,A> & other, const allocator_type & alloc)
	: base(other.base, alloc),
	  type(other.type)
	{
	}

	/// Construct from vector
	limit_vector(const vector<T,A> & other)
	: base(other),
	  type(limit_type::GENERAL)
	{
		_check_if_valid_construction();
	}
	/// Construct from vector and set allocator
	limit_vector(const vector<T,A> & other, const allocator_type & alloc)
	: base(other, alloc),
	  type(limit_type::GENERAL)
	{
		_check_if_valid_construction();
	}

	/// Coerce from vector
	template<typename To, typename Ao>
	limit_vector(const vector<To,Ao> & other)
	: base(other.begin(),other.end(),other.get_allocator()),
	  type(limit_type::GENERAL)
	{
		_check_if_valid_construction();
	}
	/// Coerce from vector and set allocator
	template<typename To, typename Ao>
	limit_vector(const vector<To,Ao> & other, const allocator_type& alloc)
	: base(other.begin(),other.end(),alloc),
	  type(limit_type::GENERAL)
	{
		_check_if_valid_construction();
	}

	/// Move constructor
	limit_vector(limit_vector<T,A> && other) = default;

	/// Move and set allocator
	limit_vector(limit_vector<T,A> && other, const allocator_type& alloc)
	: base(other.base, alloc),
	  type(other.type)
	{
	}

	/// Move from vector
	limit_vector(vector<T,A> && other, const allocator_type & alloc = other.get_allocator())
	: base(alloc),
	  type(other.type)
	{
		if(is_monotonically_increasing(other))
		{
			using std::swap;
			swap(base,other);
		}
		else
		{
			throw std::runtime_error("Limit vector cannot be constructed from vector which isn't monotonic increasing.\n");
		}
	}


	/// Construct from initializer list
	limit_vector(initializer_list<value_type> il,
	       const allocator_type& alloc = allocator_type())
	: base(il,alloc),
	  type(limit_type::GENERAL)
	{
		_check_if_valid_construction();
	}

#endif // Constructors

	/// Virtual destructor
	virtual ~limit_vector() {}


};

} /* namespace brgastro */

#endif /* _BRG_VECTOR_LIMIT_VECTOR_HPP_ */
