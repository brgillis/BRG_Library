/**********************************************************************\
 @file correlation_function_estimator.h
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

// body file: correlation_function_estimator.cpp
#ifndef _BRG_PHYSICS_CORRELATION_FUNCTION_ESTIMATOR_H_INCLUDED_
#define _BRG_PHYSICS_CORRELATION_FUNCTION_ESTIMATOR_H_INCLUDED_

#include <functional>
#include <utility>
#include <vector>

#include "brg/global.h"

#include "brg/vector/limit_vector.hpp"

namespace brgastro {

/**
 *
 */
class correlation_function_estimator {
private:

	typedef std::vector<std::pair<double,double>> position_list;

	position_list _D1_pos_list_, _D2_pos_list_, _R1_pos_list_, _R2_pos_list_;

	brgastro::limit_vector<double> _r_bin_limits_;

	bool _set_up();

public:

	// Constructors
#if(1)

	/// Default constructor. Will not be set up for calculation
	correlation_function_estimator()
	{
	}

	/// Construct with only bin limits. Will not set up for calculation
	template<typename T>
	correlation_function_estimator(T && D_bin_limits)
			: _r_bin_limits_(std::forward<T>(D_bin_limits))
	{
	}

	/// Construct with only position lists
	template<typename TR1, typename TR2, typename TM1, typename TM2>
	correlation_function_estimator(TR1 && D1_pos_list, TR2 && D2_pos_list,
			TM1 && R1_pos_list, TM2 && R2_pos_list)
			: _D1_pos_list_(std::forward<TR1>(D1_pos_list)),
			  _D2_pos_list_(std::forward<TR2>(D2_pos_list)),
			  _R1_pos_list_(std::forward<TM1>(R1_pos_list)),
			  _R2_pos_list_(std::forward<TM2>(R2_pos_list))
	{
	}

	/// Construct with bin limits and position lists
	template<typename T, typename TR1, typename TR2, typename TM1, typename TM2>
	correlation_function_estimator(T && D_bin_limits,
			TR1 && D1_pos_list, TR2 && D2_pos_list,
			TM1 && R1_pos_list, TM2 && R2_pos_list)
			: _D1_pos_list_(std::forward<TR1>(D1_pos_list)),
			  _D2_pos_list_(std::forward<TR2>(D2_pos_list)),
			  _R1_pos_list_(std::forward<TM1>(R1_pos_list)),
			  _R2_pos_list_(std::forward<TM2>(R2_pos_list)),
			  _r_bin_limits_(std::forward<T>(D_bin_limits))
	{
	}

#endif

	// Virtual destructor
	virtual ~correlation_function_estimator()
	{
	}

	// Calculation function
	std::valarray<double> calculate_weighted(const std::function<double(double)> &
			weight_function = [] (const double & theta) {return 1.;});
	std::valarray<double> calculate();
};

} // end namespace brgastro

#endif // _BRG_PHYSICS_CORRELATION_FUNCTION_ESTIMATOR_H_INCLUDED_
