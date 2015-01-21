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
#include <type_traits>
#include <utility>
#include <valarray>
#include <vector>

#include <boost/serialization/vector.hpp>

#include "brg/global.h"

#include "brg/vector/make_limit_vector_base.hpp"
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
	enum class type {GENERAL, LINEAR, LOG, DISJOINT}; // TODO Add DISJOINT type and support for it

private:
	/// The base vector storing the data
	std::vector<T,A> _base_;

	/// What type of limit vector this is (linear, log, or general)
	type _type_;

	/// Stored step size (to improve speed at the cost of a bit of memory)
	T _step_;

	/// Stored log of minimum value (to improve speed at the cost of a bit of memory)
	double _lmin_;

	// Private functions
#if(1)
	/// Init function - sets limits to cover full numeric range. Must not be called when base is not empty
	void _init()
	{
		_type_ = type::GENERAL;

		_base_.reserve(2);

		_base_.push_back(std::numeric_limits<T>::lowest());
		_base_.push_back(std::numeric_limits<T>::max());

		_step_ = 1;
		_lmin_ = 0;
	} // void _init()

	/// Check if is valid, and throw if it isn't
	void _check_if_valid_construction()
	{
		if(!is_monotonically_increasing(_base_))
			throw std::runtime_error("Limit vector cannot be constructed from vector which isn't monotonic increasing.\n");
	}

	/// Check if is valid, and throw if it isn't
	void _check_if_valid_setting(const std::vector<T,A> & other)
	{
		if(!is_monotonically_increasing(other))
			throw std::runtime_error("Limit vector cannot be set to vector which isn't monotonic increasing.\n");
	}

	/// Check if is valid, and throw if it isn't
	void _check_if_valid_disjoint_setting(const std::vector<T,A> & other)
	{
		if((!is_monotonically_increasing(other)) || (other.size() % 2 != 0))
			throw std::runtime_error("Limit vector cannot be made disjoint from vector which isn't monotonic increasing with even size.\n");
	}

	/// Check if the base is even size (which is necessary for it to be a valid disjoint limit vector)
	void _check_if_even_size()
	{
		if(_base_.size() % 2 != 0)
			throw std::runtime_error("Limit vector with odd size cannot be made disjoint.\n");
	}
#endif // Private functions

public:

	// Constructors
#if(1)

	/// Default constructor
	explicit limit_vector(const allocator_type& alloc = allocator_type())
	: _base_(alloc)
	{
		_init();
	}

	/// Construct as linear format
	limit_vector(const T & min, const T & max, const size_type & num_bins,
			const allocator_type& alloc = allocator_type())
	: _base_(alloc)
	{
		reconstruct(type::LINEAR,min,max,num_bins);
	}

	/// Construct as linear or log format
	limit_vector(const type & init_type, const T & min, const T & max, const size_type & num_bins,
			const allocator_type& alloc = allocator_type())
	: _base_(alloc)
	{
		reconstruct(init_type,min,max,num_bins);
	}

	/// Range constructor
	template <class InputIterator>
	limit_vector(InputIterator first, InputIterator last,
			const allocator_type& alloc = allocator_type())
	: _base_(first,last,alloc),
	  _type_(type::GENERAL),
	  _step_(1),
	  _lmin_(0)
	{
		_check_if_valid_construction();
		_base_.shrink_to_fit();
	}

	/// Copy constructor
	limit_vector(const limit_vector<T,A> & other) = default;

	/// Copy and set allocator
	limit_vector(const limit_vector<T,A> & other, const allocator_type & alloc)
	: _base_(other._base_, alloc),
	  _type_(other._type_),
	  _step_(other._step_),
	  _lmin_(other._lmin_)
	{
	}

	/// Construct from vector
	explicit limit_vector(const std::vector<T,A> & other)
	: _base_(other),
	  _type_(type::GENERAL),
	  _step_(1),
	  _lmin_(0)
	{
		_check_if_valid_construction();
	}
	/// Construct from vector and set allocator
	limit_vector(const std::vector<T,A> & other, const allocator_type & alloc)
	: _base_(other, alloc),
	  _type_(type::GENERAL),
	  _step_(1),
	  _lmin_(0)
	{
		_check_if_valid_construction();
	}

	/// Coerce from vector
	template<typename To, typename Ao>
	explicit limit_vector(const std::vector<To,Ao> & other, const type & init_type=type::GENERAL)
	: _base_(other.begin(),other.end(),other.get_allocator()),
	  _type_(init_type),
	  _step_(1),
	  _lmin_(0)
	{
		// If Linear or log type was passed, switch it to general, since we can't trust it'll be safe
		if(_type_==type::LINEAR) _type_ = type::GENERAL;
		if(_type_==type::LOG) _type_ = type::GENERAL;

		// If Disjoint type was passed, check that it's even size
		if(_type_==type::DISJOINT) _check_if_even_size();

		_check_if_valid_construction();
		_base_.shrink_to_fit();
	}
	/// Coerce from vector and set allocator
	template<typename To, typename Ao>
	limit_vector(const std::vector<To,Ao> & other, const allocator_type& alloc, const type & init_type=type::GENERAL)
	: _base_(other.begin(),other.end(),alloc),
	  _type_(init_type),
	  _step_(1),
	  _lmin_(0)
	{
		// If Linear or log type was passed, switch it to general, since we can't trust it'll be safe
		if(_type_==type::LINEAR) _type_ = type::GENERAL;
		if(_type_==type::LOG) _type_ = type::GENERAL;

		// If Disjoint type was passed, check that it's even size
		if(_type_==type::DISJOINT) _check_if_even_size();

		_check_if_valid_construction();
		_base_.shrink_to_fit();
	}

	/// Move constructor
	limit_vector(limit_vector<T,A> && other) = default;

	/// Move and set allocator
	limit_vector(limit_vector<T,A> && other, const allocator_type& alloc)
	: _base_(std::move(other._base_), alloc),
	  _type_(std::move(other._type_)),
	  _step_(std::move(other._step_)),
	  _lmin_(std::move(other._lmin_))
	{
	}

	/// Move from vector
	explicit limit_vector(std::vector<T,A> && other)
	: _base_(other.get_allocator()),
	  _type_(type::GENERAL),
	  _step_(1),
	  _lmin_(0)
	{
		if(is_monotonically_increasing(other))
		{
			using std::swap;
			swap(_base_,other);
		}
		else
		{
			throw std::runtime_error("Limit vector cannot be constructed from vector which isn't monotonic increasing.\n");
		}
	}

	/// Move from vector and set allocator
	limit_vector(std::vector<T,A> && other, const allocator_type& alloc)
	: _base_(alloc),
	  _type_(type::GENERAL),
	  _step_(1),
	  _lmin_(0)
	{
		if(is_monotonically_increasing(other))
		{
			using std::swap;
			swap(_base_,other);
		}
		else
		{
			throw std::runtime_error("Limit vector cannot be constructed from vector which isn't monotonic increasing.\n");
		}
	}


	/// Construct from initializer list
	explicit limit_vector(const std::initializer_list<value_type> & il,
	       const allocator_type& alloc = allocator_type())
	: _base_(il,alloc),
	  _type_(type::GENERAL),
	  _step_(1),
	  _lmin_(0)
	{
		_check_if_valid_construction();
	}

#endif // Constructors

	/// Virtual destructor
	virtual ~limit_vector() {}

	// Reconstruction
#if(1)
	/// Reconstruct as linear or log format
	void reconstruct(const type & init_type, const T & min, const T & max,
			const size_type & num_bins)
	{
		if(num_bins<1)
			throw std::logic_error("Limit vector cannot be constructed with zero bins.\n");

		if(init_type==type::LOG)
		{
			_type_ = init_type;
			_base_ = make_log_limit_vector_base<T,A>(min,max,num_bins);

			using std::pow;
			using std::log;
			_lmin_ = log(min);

			_step_ = (log(max)-_lmin_)/num_bins;

		}
		else
		{
			_type_ = type::LINEAR;
			_base_ = make_linear_limit_vector_base<T,A>(min,max,num_bins);

			_step_ = (max-min)/num_bins;
			_lmin_ = 0;
		}
		_base_.shrink_to_fit();
	}

	/// Reconstruct from a vector of the bin middles. Since some information is lost here,
	/// the method has to take some guesses.
	template<typename To, typename Ao>
	void reconstruct_from_bin_mids(const std::vector<To,Ao> & vec)
	{
		if(vec.size()==1)
		{
			clear();
			return;
		}

		if(!is_monotonically_increasing(vec))
			throw std::logic_error("Cannot reconstruct from a mids vector which isn't monotonically increasing.\n");

		_base_.resize(vec.size()+1);

		size_t i;

		for(i=0; i<vec.size(); ++i)
		{
			_base_[i] = (3*vec[i]-vec[i+1])/2;
		}

		// Special handling for the final element
		_base_[i] = (3*vec[i-1]-vec[i-2])/2;
		_type_ = type::GENERAL;
		return;
	}
	/// Reconstruct from a vector of the bin middles (move version). Since some information is lost here,
	/// the method has to take some guesses.
	template<typename To, typename Ao>
	void reconstruct_from_bin_mids(std::vector<To,Ao> && vec)
	{
		if(vec.size()==1)
		{
			clear();
			return;
		}

		if(!is_monotonically_increasing(vec))
			throw std::logic_error("Cannot reconstruct from a mids vector which isn't monotonically increasing.\n");

		size_t i;

		for(i=0; i<vec.size()-1; ++i)
		{
			vec[i] -= (vec[i+1]-vec[i])/2;
		}

		// Special handling for the final element
		To d_last = (vec[i]-vec[i-1])/3;
		vec.push_back(vec[i]+d_last);
		vec[i]-=d_last;

		_base_ = std::move(vec);
		_type_ = type::GENERAL;
	}
	/// Reconstruct as disjoint
	template<typename To>
	void reconstruct_as_disjoint(To && vec)
	{
		_check_if_valid_disjoint_setting(vec);

		_base_ = std::forward<To>(vec);
		_type_ = type::DISJOINT;
	}
#endif

	// Change type to general or disjoint
#if(1)

	void make_general()
	{
		_type_ = type::GENERAL;
	}

	void make_disjoint()
	{
		_check_if_even_size();
		_type_ = type::DISJOINT;
	}

#endif

	// Assignment
#if(1)
	/// Copy assignment
	limit_vector<T,A> & operator=(limit_vector<T,A> other)
	{
		swap(other);
		return *this;
	}

	/// Copy from vector
	limit_vector<T,A> & operator=(const std::vector<T,A> & new_base)
	{
		_check_if_valid_setting(new_base);

		_base_ = new_base;
		_type_ = type::GENERAL;

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
		swap(_base_,new_base);

		_type_ = type::GENERAL;

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
		return _base_.begin();
	}
	/// end (only const version allowed)
	const_iterator end() const noexcept
	{
		return _base_.end();
	}
	/// rbegin (only const version allowed)
	const_iterator rbegin() const noexcept
	{
		return _base_.rbegin();
	}
	/// rend (only const version allowed)
	const_iterator rend() const noexcept
	{
		return _base_.rend();
	}
	/// cbegin
	const_iterator cbegin() const noexcept
	{
		return _base_.cbegin();
	}
	/// cend
	const_iterator cend() const noexcept
	{
		return _base_.cend();
	}
	/// crbegin
	const_iterator crbegin() const noexcept
	{
		return _base_.crbegin();
	}
	/// crend
	const_iterator crend() const noexcept
	{
		return _base_.crend();
	}

#endif // Iterator methods

	// Capacity methods
#if(1)

	/// Get size of the base vector
	size_type size() const noexcept
	{
		return _base_.size();
	}

	/// Get number of bins
	size_type num_bins() const noexcept
	{
		if(_type_==type::DISJOINT)
			return _base_.size()/2;
		else
			return _base_.size()-1;
	}

	/// Get max size of the base vector
	size_type max_size() const noexcept
	{
		return _base_.max_size();
	}

	/// Get capacity of the base vector
	size_type capacity() const noexcept
	{
		return _base_.capacity();
	}

	/// Empty test - will always be false (included here for compatibility)
	bool empty() const noexcept
	{
		return false;
	}

	/// Request a change in capacity of the base vector
	void reserve(const size_type & n)
	{
		_base_.reserve(n);
	}

	/// Reduce base capacity to fit its size
	void shrink_to_fit() const
	{
		_base_.shrink_to_fit();
	}

#endif // Capacity methods

	// Element access
#if(1)

	/// Element access (const only)
	const_reference operator[] (const size_type & n) const
	{
		return _base_[n];
	}

	/// Range-checked element access (const only)
	const_reference at( const size_type & n ) const
	{
		return _base_.at(n);
	}

	/// Access first element (const only)
	const_reference front() const
	{
		return _base_.front();
	}

	/// Access last element (const only)
	const_reference back() const
	{
		return _base_.back();
	}

	/// Access data (const only)
	const value_type* data() const noexcept
	{
		return _base_.data();
	}

#endif // Element access

	// Modifiers
#if(1)

	/// Swap with another limit_vector
	void swap(limit_vector<T,A> & other)
	{
		using std::swap;
		swap(_base_,other._base_);
		swap(_type_,other._type_);
		swap(_lmin_,other._lmin_);
		swap(_step_,other._step_);
	}

	/// Clear function - leaves it in initial state
	void clear()
	{
		_base_.clear();
		_init();
	}

	void fixbad()
	{
		for(auto & val : _base_)
			brgastro::fixbad(val);
	}
#endif // Modifiers

	// Limit_vector-specific functions
#if(1)

	/// Returns minimum limit
	const_reference min() const noexcept
	{
		return _base_.front();
	}
	/// Returns maximum limit
	const_reference max() const noexcept
	{
		return _base_.back();
	}

	/// Get lower limit of specific bin
	const_reference lower_limit( const size_type & i ) const noexcept
	{
		if(_type_ != type::DISJOINT)
		{
			return _base_[i];
		}
		else
		{
			return _base_[2*i];
		}
	}

	/// Get upper limit of specific bin
	const_reference upper_limit( const size_type & i ) const noexcept
	{
		if(_type_ != type::DISJOINT)
		{
			return _base_[i+1];
		}
		else
		{
			return _base_[2*i+1];
		}
	}

	/// Returns true if val is equal to or greater than the maximum limit of this vector
	bool above_limits(const T & val) const noexcept
	{
		return val>=max();
	}

	/// Returns true if val is less than the minimum limit of this vector
	bool under_limits(const T & val) const noexcept
	{
		return val<min();
	}

	/// Returns true if val is outside the limits of this vector. Equivalent to
	/// above_limits OR below_limits.
	bool outside_limits(const T & val) const noexcept
	{
		if (above_limits(val) or under_limits(val)) return true;
		if(_type_ != type::DISJOINT) return false;

		// If we get here, it's a disjoint vector, and it's between the upper and lower limits
		// We'll check which "general" bin it falls in. If it's in a gap, it's outside the limits,
		// otherwise it's inside

		for(size_t i=1; i<size(); ++i)
		{
			if(_base_[i]>=val)
			{
				if ( (i-1)%2 == 0 )
					return false;
				else
					return true;
			}
		}
		return true; // We shouldn't reach this path, but if we do, it's outside
	}

	/// Returns true if val is inside the limits of the vector
	bool inside_limits(const T & val) const noexcept
	{
		return !outside_limits(val);
	}

	/// Gets the bin number val falls within. 0 implies lowest valid bin, num_bins()-1 implies
	/// highest valid bin, num_bins() implies below limits, num_bins()+1 implies above limits.
	size_type get_bin_index(const T & val) const
	{
		if(under_limits(val)) return num_bins();
		if(above_limits(val)) return num_bins()+1;

		switch(_type_)
		{
		case type::LINEAR:

			{
				if(std::is_integral<T>::value)
				{
					return static_cast<size_type>((val-min())/_step_);
				}
				else
				{
					// For floating-point types, we add an epsilon at this step to prevent values at the exact bin boundaries
					// from being assigned the wrong bin
					return static_cast<size_type>((val-min()+std::numeric_limits<T>::epsilon())/_step_);
				}
			}

			break;

		case type::LOG:

			{
				using std::log;
				T lval = log(val);
				return static_cast<size_type>((lval-_lmin_+std::numeric_limits<T>::epsilon())/_step_);

			}

			break;

		case type::DISJOINT:

			for(size_t i=1; i<size(); ++i)
			{
				if(_base_[i]>=val)
				{
					if ( (i-1)%2 == 0 )
						return (i-1)/2;
					else
						return num_bins()+2;
				}
			}

			break;

		default: // limit_type::GENERAL

			for(size_t i=1; i<size(); ++i)
			{
				if(_base_[i]>=val) return i-1;
			}

			break;
		}

		assert(false); // If we reach this path, there's an error
		return 0; // Just to suppress editor errors
	}

	/// Interpolate between values for successive bins
	template<typename Tv, typename To>
	Tv interpolate_bins(const Tv & val, const To & val_vec) const
	{
		if(val_vec.size()!=num_bins())
			throw std::logic_error("Value vector's size must equal num_bins() in interpolate_bins.\n");

		if(num_bins()==1)
		{
			return val_vec[0];
		}

		size_t bin_i=size()-3;

		const char step_length = _type_==type::DISJOINT ? 2 : 1;

		for(size_t i=0; i<size()-2; i += step_length)
		{
			if((_base_[i]+_base_[i+1])/2>=val)
			{
				if(i==0)
					bin_i=0;
				else
					bin_i=i-1;
				break;
			}
		}

		T xlo = (_base_[bin_i]+_base_[bin_i+1])/2;
		T xhi = (_base_[bin_i+step_length]+_base_[bin_i+step_length+1])/2;
		const Tv & ylo = val_vec[bin_i];
		const Tv & yhi = val_vec[bin_i+1];

		return ylo + (yhi-ylo)/(xhi-xlo) * (val-xlo);
	}

	std::vector<T,A> get_bin_mids() const
	{
		if(_type_!=type::DISJOINT)
		{
			std::vector<T,A> result(_base_);
			for(size_t i=0; i<num_bins(); ++i)
			{
				result[i] += (result[i+1]-result[i])/2;
			}

			result.pop_back();
			return result;
		}
		else
		{
			std::vector<T,A> result(num_bins());
			for(size_t i=0; i<size(); i+=2)
			{
				result.push_back((_base_[i]+_base_[i+1])/2);
			}
			return result;
		}
	}

#endif

	// Casting
#if(1)

	/// Cast to vector
	operator std::vector<T,A>() const
	{
		return _base_;
	}

	/// Coerce to vector
	template<typename To, typename Ao>
	std::vector<To,Ao> to_vector() const
	{
		std::vector<To,Ao> result;
		make_vector_coerce<1>(result,_base_);
		return result;
	}

	/// Coerce to valarray
	template<typename To>
	std::valarray<To> to_valarray() const
	{
		std::valarray<To> result;
		make_vector_coerce<1>(result,_base_);
		return result;
	}

#endif

	// Operator overloads
#if(1)

	bool operator==(const brgastro::limit_vector<T,A> & other) const
	{
		if((_type_==type::DISJOINT)||(other._type_==type::DISJOINT))
			return (_base_==other._base_)&&(_type_==other._type_);
		else
			return _base_==other._base_;
	}

#endif

	// Serialization (to allow it to be saved)
#if(1)
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    	ar & _base_;
    	ar & _type_;
    	ar & _step_;
    	ar & _lmin_;
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

namespace brgastro {
template<typename T, typename A=std::allocator<T>>
bool is_monotonically_increasing(const brgastro::limit_vector<T,A> &v)
{
	return true;
}
}

#endif // Non-member function overloads for limit_vector

#endif /* _BRG_VECTOR_LIMIT_VECTOR_HPP_ */
