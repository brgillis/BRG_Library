/**********************************************************************\
 @file functor.hpp
 ------------------

 Abstract base class for a functor. Inherit from this to ensure
 compatibility with all functions in the library.

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


#ifndef _BRG_FUNCTOR_HPP_INCLUDED_
#define _BRG_FUNCTOR_HPP_INCLUDED_

#include <vector>

#include <boost/call_traits.hpp>

#include "brg/global.h"

namespace brgastro {

template<typename T, typename param_struct=bool> // Defaults to using the minimum size for the param structure
class functor
{
private:

	param_struct _params_;

public:

	functor()
	{
	}

	functor(param_struct init_params)
	: _params_(init_params)
	{
	}

	virtual ~functor()
	{
	}

	virtual void set_params(typename boost::call_traits<param_struct>::param_type new_params)
	{
		_params_ = new_params;
	}

	typename boost::call_traits<param_struct>::const_reference params() const
	{
		return _params_;
	}

	virtual typename boost::call_traits<T>::value_type
		operator()(const typename boost::call_traits<T>::param_type in_param, const bool silent=false) const =0;

};

// Vector specialisation
template<typename T, typename param_struct>
class functor<std::vector<T>,param_struct>
{
private:

	param_struct _params_;

public:

	functor()
	{
	}

	functor(param_struct init_params)
	: _params_(init_params)
	{
	}

	virtual ~functor()
	{
	}

	virtual void set_params(typename boost::call_traits<param_struct>::param_type new_params)
	{
		_params_ = new_params;
	}

	typename boost::call_traits<param_struct>::const_reference params() const
	{
		return _params_;
	}

	virtual typename boost::call_traits<std::vector<T>>::value_type
		operator()(const typename boost::call_traits<std::vector<T>>::param_type in_params, const bool silent=false) const =0;

};

} // namespace brgastro

#endif // _BRG_FUNCTOR_HPP_INCLUDED_
