/**********************************************************************\
 @file tuple.hpp
 ------------------

 Element-wise operations for tuples of arbitrary length.

 **********************************************************************

 Copyright (C) 2015 brg

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

#ifndef BRG_CONTAINER_TUPLE_HPP_
#define BRG_CONTAINER_TUPLE_HPP_

#include <type_traits>
#include <utility>

#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/is_convertible.hpp>

#include "brg/container/is_boost_tuple.hpp"
#include "brg/vector/elementwise_functions.hpp" // So we have primary definition of math funcs

namespace brgastro {

using boost::tuple;

namespace tuples { // Various helper functions the end-user won't need to worry about

/// Helper function to construct a cons list with argument deduction
template<typename Th, typename Tt>
boost::tuples::cons<Th,Tt> make_cons(Th && h, Tt && t)
{
	return boost::tuples::cons<Th,Tt>(std::forward<Th>(h),std::forward<Tt>(t));
}

// add_typeof_helper
#if(1)

/**
 * Helper structure to determine the type of adding two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2, class Enable = void>
struct add_typeof_helper
{
};

/**
 * Helper structure to determine the type of adding two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2>
struct add_typeof_helper<T1,T2,typename std::enable_if<
	is_boost_tuple<T1>::value && is_boost_tuple<T2>::value>::type>
{
	typedef typename std::decay<T1>::type T1d;
	typedef typename std::decay<T2>::type T2d;

	typedef boost::tuples::cons<decltype(add(T1d().get_head(),T2d().get_head())),
			typename add_typeof_helper<decltype(T1d().get_tail()),decltype(T2d().get_tail())>::type> type;
};

/**
 * Helper structure to determine the type of adding two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2>
struct add_typeof_helper<T1,T2,typename std::enable_if<
	is_boost_tuple<T1>::value && !is_boost_tuple<T2>::value>::type>
{
	typedef typename std::decay<T1>::type T1d;
	typedef typename std::decay<T2>::type T2d;

	typedef boost::tuples::cons<decltype(add(T1d().get_head(),T2d())),
			typename add_typeof_helper<decltype(T1d().get_tail()),T2>::type> type;
};

/**
 * Helper structure to determine the type of adding two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2>
struct add_typeof_helper<T1,T2,typename std::enable_if<
	!is_boost_tuple<T1>::value && is_boost_tuple<T2>::value>::type>
{
	typedef typename std::decay<T1>::type T1d;
	typedef typename std::decay<T2>::type T2d;

	typedef boost::tuples::cons<decltype(add(T1d(),T2d().get_head())),
			typename add_typeof_helper<T1d,decltype(T2d().get_tail())>::type> type;
};

/**
 * Null type overload to end recursion.
 */
template<typename T>
struct add_typeof_helper<boost::tuples::null_type,T,void>
{
	typedef boost::tuples::null_type type;
};
/**
 * Null type overload to end recursion.
 */
template<typename T>
struct add_typeof_helper<T,boost::tuples::null_type,void>
{
	typedef boost::tuples::null_type type;
};
/**
 * Null type overload to end recursion.
 */
template<>
struct add_typeof_helper<boost::tuples::null_type,boost::tuples::null_type,void>
{
	typedef boost::tuples::null_type type;
};

#endif // add_typeof_helper

// subtract_typeof_helper
#if(1)

/**
 * Helper structure to determine the type of subtracting two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2, class Enable = void>
struct subtract_typeof_helper
{
};

/**
 * Helper structure to determine the type of subtracting two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2>
struct subtract_typeof_helper<T1,T2,typename std::enable_if<
	is_boost_tuple<T1>::value && is_boost_tuple<T2>::value>::type>
{
	typedef typename std::decay<T1>::type T1d;
	typedef typename std::decay<T2>::type T2d;

	typedef boost::tuples::cons<decltype(subtract(T1d().get_head(),T2d().get_head())),
			typename subtract_typeof_helper<decltype(T1d().get_tail()),decltype(T2d().get_tail())>::type> type;
};

/**
 * Helper structure to determine the type of subtracting two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2>
struct subtract_typeof_helper<T1,T2,typename std::enable_if<
	is_boost_tuple<T1>::value && !is_boost_tuple<T2>::value>::type>
{
	typedef typename std::decay<T1>::type T1d;
	typedef typename std::decay<T2>::type T2d;

	typedef boost::tuples::cons<decltype(subtract(T1d().get_head(),T2d())),
			typename subtract_typeof_helper<decltype(T1d().get_tail()),T2>::type> type;
};

/**
 * Helper structure to determine the type of subtracting two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2>
struct subtract_typeof_helper<T1,T2,typename std::enable_if<
	!is_boost_tuple<T1>::value && is_boost_tuple<T2>::value>::type>
{
	typedef typename std::decay<T1>::type T1d;
	typedef typename std::decay<T2>::type T2d;

	typedef boost::tuples::cons<decltype(subtract(T1d(),T2d().get_head())),
			typename subtract_typeof_helper<T1d,decltype(T2d().get_tail())>::type> type;
};

/**
 * Null type overload to end recursion.
 */
template<typename T>
struct subtract_typeof_helper<boost::tuples::null_type,T,void>
{
	typedef boost::tuples::null_type type;
};
/**
 * Null type overload to end recursion.
 */
template<typename T>
struct subtract_typeof_helper<T,boost::tuples::null_type,void>
{
	typedef boost::tuples::null_type type;
};
/**
 * Null type overload to end recursion.
 */
template<>
struct subtract_typeof_helper<boost::tuples::null_type,boost::tuples::null_type,void>
{
	typedef boost::tuples::null_type type;
};

#endif // subtract_typeof_helper

// multiply_typeof_helper
#if(1)

/**
 * Helper structure to determine the type of multiplying two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2, class Enable = void>
struct multiply_typeof_helper
{
};

/**
 * Helper structure to determine the type of multiplying two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2>
struct multiply_typeof_helper<T1,T2,typename std::enable_if<
	is_boost_tuple<T1>::value && is_boost_tuple<T2>::value>::type>
{
	typedef typename std::decay<T1>::type T1d;
	typedef typename std::decay<T2>::type T2d;

	typedef boost::tuples::cons<decltype(multiply(T1d().get_head(),T2d().get_head())),
			typename multiply_typeof_helper<decltype(T1d().get_tail()),decltype(T2d().get_tail())>::type> type;
};

/**
 * Helper structure to determine the type of multiplying two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2>
struct multiply_typeof_helper<T1,T2,typename std::enable_if<
	is_boost_tuple<T1>::value && !is_boost_tuple<T2>::value>::type>
{
	typedef typename std::decay<T1>::type T1d;
	typedef typename std::decay<T2>::type T2d;

	typedef boost::tuples::cons<decltype(multiply(T1d().get_head(),T2d())),
			typename multiply_typeof_helper<decltype(T1d().get_tail()),T2>::type> type;
};

/**
 * Helper structure to determine the type of multiplying two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2>
struct multiply_typeof_helper<T1,T2,typename std::enable_if<
	!is_boost_tuple<T1>::value && is_boost_tuple<T2>::value>::type>
{
	typedef typename std::decay<T1>::type T1d;
	typedef typename std::decay<T2>::type T2d;

	typedef boost::tuples::cons<decltype(multiply(T1d(),T2d().get_head())),
			typename multiply_typeof_helper<T1d,decltype(T2d().get_tail())>::type> type;
};

/**
 * Null type overload to end recursion.
 */
template<typename T>
struct multiply_typeof_helper<boost::tuples::null_type,T,void>
{
	typedef boost::tuples::null_type type;
};
/**
 * Null type overload to end recursion.
 */
template<typename T>
struct multiply_typeof_helper<T,boost::tuples::null_type,void>
{
	typedef boost::tuples::null_type type;
};
/**
 * Null type overload to end recursion.
 */
template<>
struct multiply_typeof_helper<boost::tuples::null_type,boost::tuples::null_type,void>
{
	typedef boost::tuples::null_type type;
};

#endif // multiply_typeof_helper

// divide_typeof_helper
#if(1)

/**
 * Helper structure to determine the type of divideing two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2, class Enable = void>
struct divide_typeof_helper
{
};

/**
 * Helper structure to determine the type of divideing two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2>
struct divide_typeof_helper<T1,T2,typename std::enable_if<
	is_boost_tuple<T1>::value && is_boost_tuple<T2>::value>::type>
{
	typedef typename std::decay<T1>::type T1d;
	typedef typename std::decay<T2>::type T2d;

	typedef boost::tuples::cons<decltype(divide(T1d().get_head(),T2d().get_head())),
			typename divide_typeof_helper<decltype(T1d().get_tail()),decltype(T2d().get_tail())>::type> type;
};

/**
 * Helper structure to determine the type of divideing two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2>
struct divide_typeof_helper<T1,T2,typename std::enable_if<
	is_boost_tuple<T1>::value && !is_boost_tuple<T2>::value>::type>
{
	typedef typename std::decay<T1>::type T1d;
	typedef typename std::decay<T2>::type T2d;

	typedef boost::tuples::cons<decltype(divide(T1d().get_head(),T2d())),
			typename divide_typeof_helper<decltype(T1d().get_tail()),T2>::type> type;
};

/**
 * Helper structure to determine the type of divideing two values or tuples together.
 * This is needed since the compiler won't fully recurse function decltypes, but it
 * will fully recurse a structure.
 */
template< class T1, class T2>
struct divide_typeof_helper<T1,T2,typename std::enable_if<
	!is_boost_tuple<T1>::value && is_boost_tuple<T2>::value>::type>
{
	typedef typename std::decay<T1>::type T1d;
	typedef typename std::decay<T2>::type T2d;

	typedef boost::tuples::cons<decltype(divide(T1d(),T2d().get_head())),
			typename divide_typeof_helper<T1d,decltype(T2d().get_tail())>::type> type;
};

/**
 * Null type overload to end recursion.
 */
template<typename T>
struct divide_typeof_helper<boost::tuples::null_type,T,void>
{
	typedef boost::tuples::null_type type;
};
/**
 * Null type overload to end recursion.
 */
template<typename T>
struct divide_typeof_helper<T,boost::tuples::null_type,void>
{
	typedef boost::tuples::null_type type;
};
/**
 * Null type overload to end recursion.
 */
template<>
struct divide_typeof_helper<boost::tuples::null_type,boost::tuples::null_type,void>
{
	typedef boost::tuples::null_type type;
};

#endif // divide_typeof_helper

} // namespace tuples

// for_each
#if(1)

// Unary for_each
#if(1)

template<typename f>
inline void unary_for_each( const f & func, const boost::tuples::null_type &)
{
	return;
}

template <typename f, class Th, class Tt>
inline void unary_for_each( const f & func, const boost::tuples::cons<Th, Tt> & t1)
{
	func(t1.get_head()); // Apply to head member

	unary_for_each(func,t1.get_tail()); // Recurse to next
}

template <typename f, class Th, class Tt>
inline void unary_for_each( const f & func, boost::tuples::cons<Th, Tt> & t1)
{
	func(t1.get_head()); // Apply to head member

	unary_for_each(func,t1.get_tail()); // Recurse to next
}

#endif // Unary for_each

// Binary for_each
#if(1)

template<typename f>
inline void binary_for_each( const f & func, const boost::tuples::null_type &,
		const boost::tuples::null_type &)
{
	return;
}

template <typename f, class Th, class Tt, typename T2>
inline void binary_for_each( const f & func, const boost::tuples::cons<Th, Tt> & t1,
		T2 && t2)
{
	func(t1.get_head(),t2.get_head()); // Apply to head members

	binary_for_each(func,t1.get_tail(),t2.get_tail()); // Recurse to next
}

template <typename f, class Th, class Tt, typename T2>
inline void binary_for_each( const f & func, boost::tuples::cons<Th, Tt> & t1,
		T2 && t2)
{
	func(t1.get_head(),t2.get_head()); // Apply to head members

	binary_for_each(func,t1.get_tail(),t2.get_tail()); // Recurse to next
}

#endif // Binary for_each

// Trinary for_each
#if(1)

template<typename f>
inline void trinary_for_each( const f & func, const boost::tuples::null_type &,
		const boost::tuples::null_type &)
{
	return;
}

template <typename f, class Th, class Tt, typename T2, typename T3>
inline void trinary_for_each( const f & func, const boost::tuples::cons<Th, Tt> & t1,
		T2 && t2, T3 && t3)
{
	func(t1.get_head(),t2.get_head(),t3.get_head()); // Apply to head members

	trinary_for_each(func,t1.get_tail(),t2.get_tail(),t3.get_tail()); // Recurse to next
}

template <typename f, class Th, class Tt, typename T2, typename T3>
inline void trinary_for_each( const f & func, boost::tuples::cons<Th, Tt> & t1,
		T2 && t2, T3 && t3)
{
	func(t1.get_head(),t2.get_head(),t3.get_head()); // Apply to head members

	trinary_for_each(func,t1.get_tail(),t2.get_tail(),t3.get_tail()); // Recurse to next
}

#endif // Trinary for_each

#endif // for_each

// Basic arithmetic
#if(1)

// Addition
#if(1)

template< class T1, class T2,
	typename std::enable_if<std::is_same<typename std::decay<T1>::type,boost::tuples::null_type>::value
	|| std::is_same<typename std::decay<T2>::type,boost::tuples::null_type>::value,
		char>::type = 0>
inline boost::tuples::null_type add( T1 &&,
		T2 &&)
{
	return boost::tuples::null_type();
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<is_boost_tuple<T2>::value,char>::type = 0>
inline typename tuples::add_typeof_helper<T1,T2>::type add( const T1 & t1,
		const T2 & t2)
{
	return tuples::make_cons(add(t1.get_head(),t2.get_head()),add(t1.get_tail(),t2.get_tail()));
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<!is_boost_tuple<T2>::value,char>::type = 0>
inline typename tuples::add_typeof_helper<T1,T2>::type add( const T1 & t1,
		const T2 & t2)
{
	return tuples::make_cons(add(t1.get_head(),t2),add(t1.get_tail(),t2));
}

template <class T1, class T2,
typename std::enable_if<!is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<is_boost_tuple<T2>::value,char>::type = 0>
inline typename tuples::add_typeof_helper<T1,T2>::type add( const T1 & t1,
		const T2 & t2)
{
	return tuples::make_cons(add(t1,t2.get_head()),add(t1,t2.get_tail()));
}

#endif // Addition

// Subtraction
#if(1)

template< class T1, class T2,
	typename std::enable_if<std::is_same<typename std::decay<T1>::type,boost::tuples::null_type>::value
	|| std::is_same<typename std::decay<T2>::type,boost::tuples::null_type>::value,
		char>::type = 0>
inline boost::tuples::null_type subtract( T1 &&,
		T2 &&)
{
	return boost::tuples::null_type();
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<is_boost_tuple<T2>::value,char>::type = 0>
inline typename tuples::subtract_typeof_helper<T1,T2>::type subtract( const T1 & t1,
		const T2 & t2)
{
	return tuples::make_cons(subtract(t1.get_head(),t2.get_head()),subtract(t1.get_tail(),t2.get_tail()));
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<!is_boost_tuple<T2>::value,char>::type = 0>
inline typename tuples::subtract_typeof_helper<T1,T2>::type subtract( const T1 & t1,
		const T2 & t2)
{
	return tuples::make_cons(subtract(t1.get_head(),t2),subtract(t1.get_tail(),t2));
}

template <class T1, class T2,
typename std::enable_if<!is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<is_boost_tuple<T2>::value,char>::type = 0>
inline typename tuples::subtract_typeof_helper<T1,T2>::type subtract( const T1 & t1,
		const T2 & t2)
{
	return tuples::make_cons(subtract(t1,t2.get_head()),subtract(t1,t2.get_tail()));
}

#endif // Subtraction

// Multiplication
#if(1)

template< class T1, class T2,
	typename std::enable_if<std::is_same<typename std::decay<T1>::type,boost::tuples::null_type>::value
	|| std::is_same<typename std::decay<T2>::type,boost::tuples::null_type>::value,
		char>::type = 0>
inline boost::tuples::null_type multiply( T1 &&,
		T2 &&)
{
	return boost::tuples::null_type();
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<is_boost_tuple<T2>::value,char>::type = 0>
inline typename tuples::multiply_typeof_helper<T1,T2>::type multiply( const T1 & t1,
		const T2 & t2)
{
	return tuples::make_cons(multiply(t1.get_head(),t2.get_head()),multiply(t1.get_tail(),t2.get_tail()));
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<!is_boost_tuple<T2>::value,char>::type = 0>
inline typename tuples::multiply_typeof_helper<T1,T2>::type multiply( const T1 & t1,
		const T2 & t2)
{
	return tuples::make_cons(multiply(t1.get_head(),t2),multiply(t1.get_tail(),t2));
}

template <class T1, class T2,
typename std::enable_if<!is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<is_boost_tuple<T2>::value,char>::type = 0>
inline typename tuples::multiply_typeof_helper<T1,T2>::type multiply( const T1 & t1,
		const T2 & t2)
{
	return tuples::make_cons(multiply(t1,t2.get_head()),multiply(t1,t2.get_tail()));
}

#endif // Multiplication

// Division
#if(1)

template< class T1, class T2,
	typename std::enable_if<std::is_same<typename std::decay<T1>::type,boost::tuples::null_type>::value
	|| std::is_same<typename std::decay<T2>::type,boost::tuples::null_type>::value,
		char>::type = 0>
inline boost::tuples::null_type divide( T1 &&,
		T2 &&)
{
	return boost::tuples::null_type();
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<is_boost_tuple<T2>::value,char>::type = 0>
inline typename tuples::divide_typeof_helper<T1,T2>::type divide( const T1 & t1,
		const T2 & t2)
{
	return tuples::make_cons(divide(t1.get_head(),t2.get_head()),divide(t1.get_tail(),t2.get_tail()));
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<!is_boost_tuple<T2>::value,char>::type = 0>
inline typename tuples::divide_typeof_helper<T1,T2>::type divide( const T1 & t1,
		const T2 & t2)
{
	return tuples::make_cons(divide(t1.get_head(),t2),divide(t1.get_tail(),t2));
}

template <class T1, class T2,
typename std::enable_if<!is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<is_boost_tuple<T2>::value,char>::type = 0>
inline typename tuples::divide_typeof_helper<T1,T2>::type divide( const T1 & t1,
		const T2 & t2)
{
	return tuples::make_cons(divide(t1,t2.get_head()),divide(t1,t2.get_tail()));
}

#endif // Division

#endif // Basic arithmetic

// Compound assignment-arithmetic
#if(1)

// Addition
#if(1)

template< class T1, class T2,
	typename std::enable_if<std::is_same<typename std::decay<T1>::type,boost::tuples::null_type>::value,
		char>::type = 0>
inline const T1 & add_equals( const T1 & t1,
		T2 &&)
{
	return t1;
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<is_boost_tuple<T2>::value,char>::type = 0>
inline T1 & add_equals( T1 & t1,
		const T2 & t2)
{
	static_assert(boost::is_convertible<typename tuples::add_typeof_helper<T1,T2>::type,
			typename std::decay<T1>::type>::value,
			"add_equals cannot be compiled due to incompatible types.");

	add_equals(t1.get_head(),t2.get_head()); // Add heads
	add_equals(t1.get_tail(),t2.get_tail()); // Recursively add tails

	return t1;
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<!is_boost_tuple<T2>::value,char>::type = 0>
inline T1 & add_equals( T1 & t1,
		const T2 & t2)
{
	static_assert(boost::is_convertible<typename tuples::add_typeof_helper<T1,T2>::type,
			typename std::decay<T1>::type>::value,
			"add_equals cannot be compiled due to incompatible types.");

	add_equals(t1.get_head(),t2); // Add to head
	add_equals(t1.get_tail(),t2); // Recursively add to tail

	return t1;
}

#endif // Addition

// Subtraction
#if(1)

template< class T1, class T2,
	typename std::enable_if<std::is_same<typename std::decay<T1>::type,boost::tuples::null_type>::value,
		char>::type = 0>
inline const T1 & subtract_equals( const T1 & t1,
		T2 &&)
{
	return t1;
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<is_boost_tuple<T2>::value,char>::type = 0>
inline T1 & subtract_equals( T1 & t1,
		const T2 & t2)
{
	static_assert(boost::is_convertible<typename tuples::subtract_typeof_helper<T1,T2>::type,
			typename std::decay<T1>::type>::value,
			"subtract_equals cannot be compiled due to incompatible types.");

	subtract_equals(t1.get_head(),t2.get_head()); // subtract heads
	subtract_equals(t1.get_tail(),t2.get_tail()); // Recursively subtract tails

	return t1;
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<!is_boost_tuple<T2>::value,char>::type = 0>
inline T1 & subtract_equals( T1 & t1,
		const T2 & t2)
{
	static_assert(boost::is_convertible<typename tuples::subtract_typeof_helper<T1,T2>::type,
			typename std::decay<T1>::type>::value,
			"subtract_equals cannot be compiled due to incompatible types.");

	subtract_equals(t1.get_head(),t2); // subtract from head
	subtract_equals(t1.get_tail(),t2); // Recursively subtract from tail

	return t1;
}

#endif // Subtraction

// Multiplication
#if(1)

template< class T1, class T2,
	typename std::enable_if<std::is_same<typename std::decay<T1>::type,boost::tuples::null_type>::value,
		char>::type = 0>
inline const T1 & multiply_equals( const T1 & t1,
		T2 &&)
{
	return t1;
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<is_boost_tuple<T2>::value,char>::type = 0>
inline T1 & multiply_equals( T1 & t1,
		const T2 & t2)
{
	static_assert(boost::is_convertible<typename tuples::multiply_typeof_helper<T1,T2>::type,
			typename std::decay<T1>::type>::value,
			"multiply_equals cannot be compiled due to incompatible types.");

	multiply_equals(t1.get_head(),t2.get_head()); // multiply heads
	multiply_equals(t1.get_tail(),t2.get_tail()); // Recursively multiply tails

	return t1;
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<!is_boost_tuple<T2>::value,char>::type = 0>
inline T1 & multiply_equals( T1 & t1,
		const T2 & t2)
{
	static_assert(boost::is_convertible<typename tuples::multiply_typeof_helper<T1,T2>::type,
			typename std::decay<T1>::type>::value,
			"multiply_equals cannot be compiled due to incompatible types.");

	multiply_equals(t1.get_head(),t2); // multiply with head
	multiply_equals(t1.get_tail(),t2); // Recursively multiply with tail

	return t1;
}

#endif // Multiplication

// Division
#if(1)

template< class T1, class T2,
	typename std::enable_if<std::is_same<typename std::decay<T1>::type,boost::tuples::null_type>::value,
		char>::type = 0>
inline const T1 & divide_equals( const T1 & t1,
		T2 &&)
{
	return t1;
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<is_boost_tuple<T2>::value,char>::type = 0>
inline T1 & divide_equals( T1 & t1,
		const T2 & t2)
{
	static_assert(boost::is_convertible<typename tuples::divide_typeof_helper<T1,T2>::type,
			typename std::decay<T1>::type>::value,
			"divide_equals cannot be compiled due to incompatible types.");

	divide_equals(t1.get_head(),t2.get_head()); // divide heads
	divide_equals(t1.get_tail(),t2.get_tail()); // Recursively divide tails

	return t1;
}

template <class T1, class T2,
typename std::enable_if<is_boost_tuple<T1>::value,char>::type = 0,
typename std::enable_if<!is_boost_tuple<T2>::value,char>::type = 0>
inline T1 & divide_equals( T1 & t1,
		const T2 & t2)
{
	static_assert(boost::is_convertible<typename tuples::divide_typeof_helper<T1,T2>::type,
			typename std::decay<T1>::type>::value,
			"divide_equals cannot be compiled due to incompatible types.");

	divide_equals(t1.get_head(),t2); // divide from head
	divide_equals(t1.get_tail(),t2); // Recursively divide from tail

	return t1;
}

#endif // Division

#endif // Compound assignment-arithmetic

}

#endif // BRG_CONTAINER_TUPLE_HPP_
