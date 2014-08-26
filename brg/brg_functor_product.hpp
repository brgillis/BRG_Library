/**********************************************************************\
brg_functpr.hpp
 -----------



 Everything in this file is declared in the namespace brgastro.

 \**********************************************************************/

#ifndef __BRG_FUNCTOR_HPP_INCLUDED__
#define __BRG_FUNCTOR_HPP_INCLUDED__

#include <iostream>
#include <vector>

#include "brg_global.h"

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
	bool _f1_set_up_, _f2_set_up_;

public:
	// Constructors

	functor_product()
	{
		_f1_ptr_ = _f2_ptr_ = 0;
		_f1_set_up_ = _f2_set_up_ = false;
	}

	functor_product( const f1 *new_f1, const f2 *new_f2 )
	{
		_f1_ptr_ = new_f1;
		_f2_ptr_ = new_f2;
		_f1_set_up_ = _f2_set_up_ = true;
	}

	functor_product( const f1 &new_f1, const f2 &new_f2 )
	{
		_f1_ptr_ = &new_f1;
		_f2_ptr_ = &new_f2;
		_f1_set_up_ = _f2_set_up_ = true;
	}

	// Virtual destructor
	virtual ~functor_product()
	{
	}

	// Set methods

	const int set_f1_ptr( const f1* new_f1_ptr )
	{
		_f1_ptr_ = new_f1_ptr;
		_f1_set_up_ = true;
		return 0;
	}
	const int set_f1_ptr( const f1 &new_f1_ptr )
	{
		_f1_ptr_ = &new_f1_ptr;
		_f1_set_up_ = true;
		return 0;
	}

	const int set_f2_ptr( const f2 *new_f2_ptr )
	{
		_f2_ptr_ = new_f2_ptr;
		_f2_set_up_ = true;
		return 0;
	}
	const int set_f2_ptr( const f2 &new_f2_ptr )
	{
		_f2_ptr_ = &new_f2_ptr;
		_f2_set_up_ = true;
		return 0;
	}

	const int set_f1_f2_ptrs( const f1 *new_f1_ptr, const f2 *new_f2_ptr )
	{
		if ( set_f1_ptr( new_f1_ptr ) )
			return 1;
		return set_f2_ptr( new_f2_ptr );
	}
	const int set_f1_f2_ptrs( const f1 &new_f1_ptr, const f2 &new_f2_ptr )
	{
		if ( set_f1_ptr( new_f1_ptr ) )
			return 1;
		return set_f2_ptr( new_f2_ptr );
	}

	// Function method

	const int operator()( const T & in_param, T & out_param,
			const bool silent = false ) const
	{
		if ( ( !_f1_set_up_ ) || ( !_f2_set_up_ ) )
		{
			return NOT_SET_UP_ERROR;
		}

		T f1_out_param, f2_out_param;

		if ( ( *_f1_ptr_ )( in_param, f1_out_param, silent ) )
			return 1;
		if ( ( *_f2_ptr_ )( in_param, f2_out_param, silent ) )
			return 1;

		out_param = f1_out_param * f2_out_param;

		return 0;
	}
};

template< typename f1, typename f2, typename T >
class functor_product< f1, f2, std::vector< T > >
{
	/**********************************************
	 functor_product
	 -------------------------

	 A specialization of functor_profuct for vector
	 input/output.

	 **********************************************/
private:

	const f1 *_f1_ptr_;
	const f2 *_f2_ptr_;
	bool _f1_set_up_, _f2_set_up_;

public:
	// Constructors

	functor_product()
	{
		_f1_ptr_ = _f2_ptr_ = 0;
		_f1_set_up_ = _f2_set_up_ = false;
	}

	functor_product( const f1 *new_f1, const f2 *new_f2 )
	{
		_f1_ptr_ = new_f1;
		_f2_ptr_ = new_f2;
		_f1_set_up_ = _f2_set_up_ = true;
	}
	functor_product( const f1 &new_f1, const f2 &new_f2 )
	{
		_f1_ptr_ = &new_f1;
		_f2_ptr_ = &new_f2;
		_f1_set_up_ = _f2_set_up_ = true;
	}

	// Virtual destructor
	virtual ~functor_product()
	{
	}

	// Set methods

	const int set_f1_ptr( const f1 *new_f1_ptr )
	{
		_f1_ptr_ = new_f1_ptr;
		_f1_set_up_ = true;
		return 0;
	}
	const int set_f1_ptr( const f1 &new_f1_ptr )
	{
		_f1_ptr_ = &new_f1_ptr;
		_f1_set_up_ = true;
		return 0;
	}

	const int set_f2_ptr( const f2 *new_f2_ptr )
	{
		_f2_ptr_ = new_f2_ptr;
		_f2_set_up_ = true;
		return 0;
	}
	const int set_f2_ptr( const f2 &new_f2_ptr )
	{
		_f2_ptr_ = &new_f2_ptr;
		_f2_set_up_ = true;
		return 0;
	}

	const int set_f1_f2_ptrs( const f1 *new_f1_ptr, const f2 *new_f2_ptr )
	{
		if ( set_f1_ptr( new_f1_ptr ) )
			return 1;
		return set_f2_ptr( new_f2_ptr );
	}
	const int set_f1_f2_ptrs( const f1 &new_f1_ptr, const f2 &new_f2_ptr )
	{
		if ( set_f1_ptr( new_f1_ptr ) )
			return 1;
		return set_f2_ptr( new_f2_ptr );
	}

	// Function method

	const int operator()( const std::vector< T > & in_params,
			std::vector< T > & out_params, const bool silent = false ) const
	{
		if ( ( !_f1_set_up_ ) || ( !_f2_set_up_ ) )
		{
			return NOT_SET_UP_ERROR;
		}

		std::vector< T > f1_out_params( 0 ), f2_out_params( 0 );

		if ( _f1_ptr_( in_params, f1_out_params, silent ) )
			return 1;
		if ( _f2_ptr_( in_params, f2_out_params, silent ) )
			return 1;

		if ( f1_out_params.size() != f2_out_params.size() )
		{
			if ( !silent )
				std::cerr
						<< "ERROR: Functions assigned to function_product_function have\n"
						<< "different numbers of output parameters.\n";
			return UNSPECIFIED_ERROR;
		}

		unsigned int num_out_params = f1_out_params.size();
		out_params.resize( num_out_params );

		for ( unsigned int i = 0; i < num_out_params; i++ )
		{
			out_params.at( i ) = f1_out_params.at( i ) * f2_out_params.at( i );
		}

		return 0;
	}
};

#endif // end class declarations

} // end namespace brgastro

#endif // __BRG_FUNCTOR_HPP_INCLUDED__
