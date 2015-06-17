/**********************************************************************\
 @file is_Eigen_container.hpp
 ------------------

 Credit to Nawaz on StackOverflow at:
 http://stackoverflow.com/a/9407521

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


#ifndef _BRG_CONTAINER_IS_EIGEN_CONTAINER_HPP_INCLUDED_
#define _BRG_CONTAINER_IS_EIGEN_CONTAINER_HPP_INCLUDED_

#include <type_traits>

#include "brg/Eigen.hpp"
#include "brg/container/is_container.hpp"

namespace brgastro
{

template<typename T>
struct is_eigen_container : std::integral_constant<bool, std::is_base_of<Eigen::DenseBase<typename std::decay<T>::type>,
	typename std::decay<T>::type>::value>
{ };

template<typename T>
struct is_stl_or_eigen_container : std::integral_constant<bool,
	is_stl_container<T>::value || is_eigen_container<T>::value>
{ };

} // end namespace brgastro



#endif // _BRG_CONTAINER_IS_EIGEN_CONTAINER_HPP_INCLUDED_
