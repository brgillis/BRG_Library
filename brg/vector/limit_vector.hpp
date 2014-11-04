/*
 * limit_vector.hpp
 *
 *  Created on: Nov 3, 2014
 *      Author: user
 */

#ifndef _BRG_VECTOR_LIMIT_VECTOR_HPP_
#define _BRG_VECTOR_LIMIT_VECTOR_HPP_

#include <cassert>
#include <cstdlib>
#include <limits>
#include <utility>
#include <valarray>
#include <vector>

#include "brg/vector/limit_vector_operations.hpp"
#include "brg/vector/make_vector.hpp"
#include "brg/vector/summary_functions.hpp"

namespace brgastro {

template < class T, class A = std::allocator<T> >
class limit_vector
{
public:
	// Some typedefs that we can allow users to access
	typedef typename std::vector<T,A>::value_type value_type;
	typedef typename std::vector<T,A>::allocator_type allocator_type;
	typedef typename std::vector<T,A>::reference reference;
	typedef typename std::vector<T,A>::const_reference const_reference;
	typedef typename std::vector<T,A>::pointer pointer;
	typedef typename std::vector<T,A>::const_pointer const_pointer;
	typedef typename std::vector<T,A>::iterator iterator;
	typedef typename std::vector<T,A>::const_iterator const_iterator;
	typedef typename std::vector<T,A>::reverse_iterator reverse_iterator;
	typedef typename std::vector<T,A>::const_reverse_iterator const_reverse_iterator;
	typedef typename std::vector<T,A>::difference_type difference_type;
	typedef typename std::vector<T,A>::size_type size_type;

	// Enum to describe type of limit vector this is
	enum class limit_type {GENERAL, LINEAR, LOG};

private:
	/// The base vector storing the data
	std::vector<T,A> base;

	/// What type of limit vector this is (linear, log, or general)
	limit_type type;

	/// Stored step size (to improve speed at the cost of a bit of memory)
	T step;

	/// Stored log of minimum value (to improve speed at the cost of a bit of memory)
	double lmin;

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

		step = 1;
		lmin = 0;
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
	limit_vector(const limit_type & init_type, const T & min, const T & max, const size_type & num_bins,
			const allocator_type& alloc = allocator_type())
	{
		if(num_bins<1)
			throw std::logic_error("Limit vector cannot be constructed with zero bins.\n");

		if(init_type==limit_type::LOG)
		{
			type = init_type;
			base = make_log_limit_vector<T>(min,max,num_bins);

			using std::pow;
			step = pow(max/min,1./num_bins);

			using std::log;
			lmin = log(min);
		}
		else
		{
			type = limit_type::LINEAR;
			base = make_limit_vector<T>(min,max,num_bins);

			step = (max-min)/num_bins;
			lmin = 0;
		}
		base.shrink_to_fit();
	}

	/// Range constructor
	template <class InputIterator>
	limit_vector(InputIterator first, InputIterator last,
			const allocator_type& alloc = allocator_type())
	: base(first,last,alloc),
	  type(limit_type::GENERAL),
	  step(1),
	  lmin(0)
	{
		_check_if_valid_construction();
		base.shrink_to_fit();
	}

	/// Copy constructor
	limit_vector(const limit_vector<T,A> & other) = default;

	/// Copy and set allocator
	limit_vector(const limit_vector<T,A> & other, const allocator_type & alloc)
	: base(other.base, alloc),
	  type(other.type),
	  step(other.step),
	  lmin(other.lmin)
	{
	}

	/// Construct from vector
	limit_vector(const std::vector<T,A> & other)
	: base(other),
	  type(limit_type::GENERAL),
	  step(1),
	  lmin(0)
	{
		_check_if_valid_construction();
	}
	/// Construct from vector and set allocator
	limit_vector(const std::vector<T,A> & other, const allocator_type & alloc)
	: base(other, alloc),
	  type(limit_type::GENERAL),
	  step(1),
	  lmin(0)
	{
		_check_if_valid_construction();
	}

	/// Coerce from vector
	template<typename To, typename Ao>
	limit_vector(const std::vector<To,Ao> & other)
	: base(other.begin(),other.end(),other.get_allocator()),
	  type(limit_type::GENERAL),
	  step(1),
	  lmin(0)
	{
		_check_if_valid_construction();
		base.shrink_to_fit();
	}
	/// Coerce from vector and set allocator
	template<typename To, typename Ao>
	limit_vector(const std::vector<To,Ao> & other, const allocator_type& alloc)
	: base(other.begin(),other.end(),alloc),
	  type(limit_type::GENERAL),
	  step(1),
	  lmin(0)
	{
		_check_if_valid_construction();
		base.shrink_to_fit();
	}

	/// Move constructor
	limit_vector(limit_vector<T,A> && other) = default;

	/// Move and set allocator
	limit_vector(limit_vector<T,A> && other, const allocator_type& alloc)
	: base(other.base, alloc),
	  type(other.type),
	  step(other.step),
	  lmin(other.lmin)
	{
	}

	/// Move from vector
	limit_vector(std::vector<T,A> && other)
	: base(other.get_allocator()),
	  type(limit_type::GENERAL),
	  step(1),
	  lmin(0)
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

	/// Move from vector and set allocator
	limit_vector(std::vector<T,A> && other, const allocator_type& alloc)
	: base(alloc),
	  type(limit_type::GENERAL),
	  step(1),
	  lmin(0)
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
	limit_vector(const std::initializer_list<value_type> & il,
	       const allocator_type& alloc = allocator_type())
	: base(il,alloc),
	  type(limit_type::GENERAL),
	  step(1),
	  lmin(0)
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
		swap(base,new_base);

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
		swap(base,new_base);

		type = limit_type::GENERAL;

		return *this;
	}

	/// Assign from initializer_list
	limit_vector<T,A> & operator= (const std::initializer_list<value_type> & il)
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

	/// Returns true if val is equal to or greater than the maximum limit of this vector
	bool above_limits(const T & val)
	{
		return val>=max();
	}

	/// Returns true if val is less than the minimum limit of this vector
	bool under_limits(const T & val)
	{
		return val<min();
	}

	/// Returns true if val is outside the limits of this vector. Equivalent to
	/// above_limits OR below_limits.
	bool outside_limits(const T & val)
	{
		return (above_limits(val) or under_limits(val));
	}

	/// Returns true if val is inside the limits of the vector (inclusive the lower
	/// limit, exclusive the upper limit).
	bool inside_limits(const T & val)
	{
		return !outside_limits(val);
	}

	/// Gets the bin number val falls within. 0 implies lowest valid bin, num_bins()-1 implies
	/// highest valid bin, num_bins() implies below limits, num_bins()+1 implies above limits.
	size_type get_bin_index(const T & val)
	{
		if(under_limits(val)) return num_bins();
		if(above_limits(val)) return num_bins()+1;

		switch(type)
		{
		case limit_type::LINEAR:

			T step = (max()-min())/num_bins();

			return static_cast<size_type>((val-min())/step);

			break;

		case limit_type::LOG:

			using std::log;
			double lmax = log(static_cast<double>(max()));
			double lmin = log(static_cast<double>(min()));
			double lval = log(static_cast<double>(val));

			double lstep = (lmax-lmin)/num_bins();

			return static_cast<size_type>((lval-lmin)/lstep);

			break;

		default: // limit_type::GENERAL

			for(size_t i=1; i<size(); ++i)
			{
				if(base[i]>=val) return i-1;
			}

			assert(false); // If we reach this path, there's an error
			return 0; // Just to suppress editor errors

			break;
		}
	}


#endif

	// Casting
#if(1)

	/// Cast to vector
	operator std::vector<T,A>()
	{
		return base;
	}

	/// Coerce to vector
	template<typename To, typename Ao>
	std::vector<To,Ao> to_vector()
	{
		std::vector<To,Ao> result;
		make_vector_coerce<1>(result,base);
		return result;
	}

	/// Coerce to valarray
	template<typename To>
	std::valarray<To> to_valarray()
	{
		std::vector<To> result;
		make_vector_coerce<1>(result,base);
		return result;
	}

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
