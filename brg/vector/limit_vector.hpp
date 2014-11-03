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
#include <utility>
#include <vector>

#include "brg/vector/limit_vector_operations.hpp"
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

		base.reserve(2);

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

	// Constructors
#if(1)

	/// Default constructor
	explicit limit_vector(const allocator_type& alloc = allocator_type())
	: base(alloc)
	{
		_init();
	}

	/// Construct as linear or log format
	limit_vector(const limit_type & init_type, const T & min, const T & max, const size_t & num_bins,
			const allocator_type& alloc = allocator_type())
	{
		if(init_type==limit_type::LOG)
		{
			type = init_type;
			base = make_log_limit_vector<T>(min,max,bins);
		}
		else
		{
			type = limit_type::LINEAR;
			base = make_linear_limit_vector<T>(min,max,bins);
		}
		base.shrink_to_fit();
	}

	/// Range constructor
	template <class InputIterator>
	limit_vector(InputIterator first, InputIterator last,
			const allocator_type& alloc = allocator_type())
	: base(first,last,alloc),
	  type(limit_type::GENERAL)
	{
		_check_if_valid_construction();
		base.shrink_to_fit();
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
		base.shrink_to_fit();
	}
	/// Coerce from vector and set allocator
	template<typename To, typename Ao>
	limit_vector(const vector<To,Ao> & other, const allocator_type& alloc)
	: base(other.begin(),other.end(),alloc),
	  type(limit_type::GENERAL)
	{
		_check_if_valid_construction();
		base.shrink_to_fit();
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

	// Assignment
#if(1)
	/// Copy assignment
	limit_vector<T,A> & operator=(limit_vector<T,A> other)
	{
		swap(other);
		return *this;
	}

	/// Copy from vector
	limit_vector<T,A> & operator=(std::vector<T,A> new_base)
	{
		_check_if_valid_setting(new_base);

		using std::swap;
		swap(base,test_base);

		type = limit_type::GENERAL;

		return *this;
	}

	/// Coerce copy from container
	template <typename ContainerType>
	limit_vector<T,A> & operator=(const ContainerType & other)
	{
		return *this = std::vector<T,A>(other.begin(),other.end());
	}

	/// Move assignment
	limit_vector<T,A> & operator=(limit_vector<T,A> && other)
	{
		swap(other);
		return *this;
	}

	/// Move from vector
	limit_vector<T,A> & operator=(std::vector<T,A> && new_base)
	{
		_check_if_valid_setting(new_base);

		using std::swap;
		swap(base,test_base);

		type = limit_type::GENERAL;

		return *this;
	}

	/// Assign from initializer_list
	limit_vector<T,A> & operator= (const initializer_list<value_type> & il)
	{
		return *this = std::vector<T,A>(il);
	}
#endif

	// Iterator methods
#if(1)
	/// begin (only const version allowed)
	const_iterator begin() const noexcept
	{
		return base.begin();
	}
	/// end (only const version allowed)
	const_iterator end() const noexcept
	{
		return base.end();
	}
	/// rbegin (only const version allowed)
	const_iterator rbegin() const noexcept
	{
		return base.rbegin();
	}
	/// rend (only const version allowed)
	const_iterator rend() const noexcept
	{
		return base.rend();
	}
	/// cbegin
	const_iterator cbegin() const noexcept
	{
		return base.cbegin();
	}
	/// cend
	const_iterator cend() const noexcept
	{
		return base.cend();
	}
	/// crbegin
	const_iterator crbegin() const noexcept
	{
		return base.crbegin();
	}
	/// crend
	const_iterator crend() const noexcept
	{
		return base.crend();
	}

#endif // Iterator methods

	// Capacity methods
#if(1)

	/// Get size of the base vector
	size_type size() const noexcept
	{
		return base.size();
	}

	/// Get number of bins
	size_type num_bins() const noexcept
	{
		return base.size()-1;
	}

	/// Get max size of the base vector
	size_type max_size() const noexcept
	{
		return base.max_size();
	}

	/// Get capacity of the base vector
	size_type capacity() const noexcept
	{
		return base.capacity();
	}

	/// Empty test - will always be false (included here for compatibility)
	bool empty() const noexcept
	{
		return false;
	}

	/// Request a change in capacity of the base vector
	void reserve(const size_type & n)
	{
		base.reserve(n);
	}

	/// Reduce base capacity to fit its size
	void shrink_to_fit() const
	{
		base.shrink_to_fit();
	}

#endif // Capacity methods

	// Element access
#if(1)

	/// Element access (const only)
	const_reference operator[] (const size_type & n) const
	{
		return base[n];
	}

	/// Range-checked element access (const only)
	const_reference at( const size_type & n ) const
	{
		return base.at(n);
	}

	/// Access first element (const only)
	const_reference front() const
	{
		return base.front();
	}

	/// Access last element (const only)
	const_reference back() const
	{
		return base.back();
	}

	/// Access data (const only)
	const value_type* data() const noexcept
	{
		return base.data();
	}

#endif // Element access

	// Modifiers
#if(1)

	/// Swap with another limit_vector
	void swap(limit_vector<T,A> & other)
	{
		using std::swap;
		swap(base,other.base);
		swap(type,other.type);
	}

	/// Clear function - leaves it in initial state
	void clear()
	{
		base.clear();
		_init();
	}

#endif // Modifiers

	// Limit_vector-specific functions
#if(1)

	/// Returns minimum limit
	const_reference min() const noexcept
	{
		return base.front();
	}
	/// Returns maximum limit
	const_reference max() const noexcept
	{
		return base.back();
	}

	/// TODO Finish this


#endif

};

} /* namespace brgastro */

// Non-member function overloads for limit_vector
#if(1)

template <typename T, typename A = std::allocator<T>>
void swap(brgastro::limit_vector<T,A> & same, brgastro::limit_vector<T,A> & other)
{
	same.swap(other);
}

namespace std
{

template <typename T, typename A = std::allocator<T>>
void swap(brgastro::limit_vector<T,A> & same, brgastro::limit_vector<T,A> & other)
{
	same.swap(other);
}

}

#endif // Non-member function overloads for limit_vector

#endif /* _BRG_VECTOR_LIMIT_VECTOR_HPP_ */
