/**       @file interpolator_derivative_functors.h
 *
 *     Project: brg
 *        Path: /brg/interpolator/interpolator_derivative_functors.h
 *
 *  Created on: 29 Aug 2014
 *      Author: brg
 */

#ifndef _INTERPOLATOR_DERIVATIVE_FUNCTORS_H_INCLUDED_
#define _INTERPOLATOR_DERIVATIVE_FUNCTORS_H_INCLUDED_

#include "brg/brg_global.h"

#include "brg/brg_units.h"
#include "brg_interpolator.h"

namespace brgastro {

class interpolator_functor
{
	/************************************************************
	 interpolator_functor
	 ---------------

	 This class is used to provide a functor for getting
	 points along a spline (which is used here to differentiate the
	 spline).

	 Use of this class is handled by the spline_derivative class.
	 No need for the end-user to worry too much about how this
	 works.

	 \************************************************************/

private:

	const brgastro::interpolator *_interpolator_ptr_;
	bool _interpolator_ptr_set_up_;

public:

	// Swap functions
	void swap(interpolator_functor & other);
	friend void swap(interpolator_functor & same, interpolator_functor & other) {same.swap(other);}

	// Constructors
	interpolator_functor();
	interpolator_functor(const interpolator_functor& other);
	interpolator_functor(const brgastro::interpolator *init_interpolator_ptr );

	// Destructor
	virtual ~interpolator_functor()
	{
	}

	// Operator=
	interpolator_functor& operator=(interpolator_functor other);

	// Set functions
	const int set_interpolator_ptr( const brgastro::interpolator *new_interpolator_ptr );

	// Function method
	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param, const bool silent = false ) const;

};
// class interpolator_functor

class interpolator_derivative_functor
{
	/************************************************************
	 interpolator_derivative_functor
	 --------------------------

	 This class is used to provide a functor * for getting
	 the derivative of an interpolator at a given point.

	 Use of this class is handled by the interpolator_derivative class.
	 No need for the end-user to worry too much about how this
	 works.

	 \************************************************************/

private:

	interpolator_functor _interpolator_functor_;
	bool _interpolator_functor_set_up_;

public:

	// Swap functions
	void swap(interpolator_derivative_functor& other);
	friend void swap(interpolator_derivative_functor& same, interpolator_derivative_functor& other)
		{same.swap(other);}

	// Constructors
	interpolator_derivative_functor();
	interpolator_derivative_functor(const interpolator_derivative_functor& other);
	interpolator_derivative_functor( brgastro::interpolator *init_interpolator_ptr );

	// Destructor
	virtual ~interpolator_derivative_functor()
	{
	}

	// Operator=
	interpolator_derivative_functor& operator=(interpolator_derivative_functor other);

	// Set functions
	const int set_interpolator_ptr( const brgastro::interpolator *new_interpolator_ptr );

	// Function method
	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param, const bool silent = false ) const;

};
// class interpolator_derivative_functor

class interpolator_derivative_weight_functor
{
	/************************************************************
	 interpolator_derivative_weight_functor
	 ---------------------------------

	 Child of functor

	 This class is used to provide a functor * for getting
	 the weight of various points in the smoothing kernel for
	 calculating the derivative of an interpolator.

	 Use of this class is handled by the interpolator_derivative class.
	 No need for the end-user to worry too much about how this
	 works.

	 \************************************************************/

private:

	double _sample_scale_, _sample_max_width_;
	double _t_min_, _t_max_, _centre_point_;

public:

	// Swap functions
	void swap(interpolator_derivative_weight_functor &other);
	friend void swap(interpolator_derivative_weight_functor &same,
			interpolator_derivative_weight_functor &other)
				{same.swap(other);}

	// Constructors
	interpolator_derivative_weight_functor();
	interpolator_derivative_weight_functor(const interpolator_derivative_weight_functor &other);

	// Destructor
	virtual ~interpolator_derivative_weight_functor()
	{
	}

	// Operator=
	interpolator_derivative_weight_functor & operator=(interpolator_derivative_weight_functor other);

	// Set functions
	const int set_sample_scale( double new_sample_scale );

	const int set_sample_max_width( double new_sample_max_width );

	const int set_center_point( double new_center_point );

	const int set_t_min( double new_t_min );

	const int set_t_max( double new_t_max );

	// Function method
	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param, const bool silent = false ) const;
};
// class interpolator_derivative_sample_functor

} // end namespace brgastro

#endif /* _INTERPOLATOR_DERIVATIVE_FUNCTORS_H_INCLUDED_ */
