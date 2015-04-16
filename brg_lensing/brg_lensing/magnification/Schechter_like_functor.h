/**********************************************************************\
 @file Schechter_like_functor.h
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

// body file: Schechter_like_functor.cpp

#ifndef _BRG_SCHECHTER_LIKE_FUNCTOR_H_INCLUDED_
#define _BRG_SCHECHTER_LIKE_FUNCTOR_H_INCLUDED_

#include <cassert>
#include <vector>

#include "brg/math/functor/functor.hpp"
#include "brg/utility.hpp"

namespace brgastro {

/**
 *
 */
class Schechter_like_functor: public brgastro::functor<long double,std::vector<long double>> {
private:
	constexpr static ssize_t _num_params_ = 7;

public:
	Schechter_like_functor()
	{
	}
	Schechter_like_functor(const std::vector<long double> & init_params )
	: functor(init_params)
	{
	}
	virtual ~Schechter_like_functor()
	{
	}

	// Params accessors
#if (1)

	long double N_scale() const
	{
		assert(ssize(params())==_num_params_);
		return params()[0];
	}
	long double m_star() const
	{
		assert(ssize(params())==_num_params_);
		return params()[1];
	}
	long double alpha() const
	{
		assert(ssize(params())==_num_params_);
		return params()[2];
	}
	long double mag_lower_lim_sharpness() const
	{
		assert(ssize(params())==_num_params_);
		return params()[3];
	}
	long double mag23_jump() const
	{
		assert(ssize(params())==_num_params_);
		return params()[4];
	}
	long double mag_upper_lim() const
	{
		assert(ssize(params())==_num_params_);
		return params()[5];
	}
	long double mag_upper_lim_sharpness() const
	{
		assert(ssize(params())==_num_params_);
		return params()[6];
	}

#endif

	// Function method

	long double operator()( const long double & in_params,
			const bool silent = false ) const;

};

} // namespace brgastro

#endif // _BRG_SCHECHTER_LIKE_FUNCTOR_H_INCLUDED_
