/**********************************************************************\
brg_functpr.hpp
 -----------



 Everything in this file is declared in the namespace brgastro.

 \**********************************************************************/

#ifndef __BRG_FUNCTOR_HPP_INCLUDED__
#define __BRG_FUNCTOR_HPP_INCLUDED__

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "brg_global.h"

#include "brg_vector_functions.hpp"

namespace brgastro
{

/** Class definitions **/
#if (1)

template< typename f1, typename f2, typename T=BRG_UNITS >
class functor_product
{
	/*****************************************************
	 functor_product
	 -----------------------------------------------------

	 An example of a functor which returns the
	 product of two other functors (see the
	 integrate_weighted* functions for
	 how this is used).

	 ****************************************************/
private:

	const f1 *_f1_ptr_;
	const f2 *_f2_ptr_;

public:
	// Constructors

	functor_product()
	{
		_f1_ptr_ = _f2_ptr_ = NULL;
	}

	functor_product( const f1 *new_f1, const f2 *new_f2 )
	{
		_f1_ptr_ = new_f1;
		_f2_ptr_ = new_f2;
	}
	functor_product( const f1 *new_f1, const f2 &new_f2 )
	{
		_f1_ptr_ = &new_f1;
		_f2_ptr_ = &new_f2;
	}
	functor_product( const f1 &new_f1, const f2 *new_f2 )
	{
		_f1_ptr_ = new_f1;
		_f2_ptr_ = new_f2;
	}
	functor_product( const f1 &new_f1, const f2 &new_f2 )
	{
		_f1_ptr_ = &new_f1;
		_f2_ptr_ = &new_f2;
	}

	// Virtual destructor
	virtual ~functor_product()
	{
	}

	// Set methods

	void set_f1_ptr( const f1* new_f1_ptr )
	{
		_f1_ptr_ = new_f1_ptr;
	}
	void set_f1_ptr( const f1 &new_f1_ptr )
	{
		_f1_ptr_ = &new_f1_ptr;
	}

	void set_f2_ptr( const f2 *new_f2_ptr )
	{
		_f2_ptr_ = new_f2_ptr;
	}
	void set_f2_ptr( const f2 &new_f2_ptr )
	{
		_f2_ptr_ = &new_f2_ptr;
	}

	void set_f1_f2_ptrs( const f1 *new_f1_ptr, const f2 *new_f2_ptr )
	{
		set_f1_ptr( new_f1_ptr );
		set_f2_ptr( new_f2_ptr );
	}
	void set_f1_f2_ptrs( const f1 *new_f1_ptr, const f2 &new_f2_ptr )
	{
		set_f1_ptr( new_f1_ptr );
		set_f2_ptr( new_f2_ptr );
	}
	void set_f1_f2_ptrs( const f1 &new_f1_ptr, const f2 *new_f2_ptr )
	{
		set_f1_ptr( new_f1_ptr );
		set_f2_ptr( new_f2_ptr );
	}
	void set_f1_f2_ptrs( const f1 &new_f1_ptr, const f2 &new_f2_ptr )
	{
		set_f1_ptr( new_f1_ptr );
		set_f2_ptr( new_f2_ptr );
	}

	// Function method

	const T operator()( const T & in_param,	const bool silent = false ) const
	{
		std::assert((_f1_ptr_!=NULL) && (_f2_ptr_!=NULL));

		return multiply(_f1_ptr_( in_param, silent ),_f2_ptr_( in_param, silent ));
	}
}; // class functor_product< f1, f2, T>

#endif // end class declarations

} // end namespace brgastro

#endif // __BRG_FUNCTOR_HPP_INCLUDED__
