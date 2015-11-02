/**********************************************************************\
 @file Eigen.hpp
 ------------------

 A file to include the Eigen library, plus set up some basic extensions
 for it - most notably, begin() and end() support

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

#ifndef _BRG_EIGEN_HPP_
#define _BRG_EIGEN_HPP_

#include <boost/serialization/split_free.hpp>
#include <Eigen/Core>

#include "common.h"

// Free begin and end definitions
#if(1)
namespace Eigen {

template <typename Derived, int Nrow, int Ncol, int Options, int Maxrow, int Maxcol>
auto begin(const Eigen::Array<Derived,Nrow,Ncol,Options,Maxrow,Maxcol> & vec) -> decltype(vec.data())
{
	return vec.data();
}

template <typename Derived, int Nrow, int Ncol, int Options, int Maxrow, int Maxcol>
auto begin(Eigen::Array<Derived,Nrow,Ncol,Options,Maxrow,Maxcol> & vec) -> decltype(vec.data())
{
	return vec.data();
}

template <typename Derived, int Nrow, int Ncol, int Options, int Maxrow, int Maxcol>
auto end(const Eigen::Array<Derived,Nrow,Ncol,Options,Maxrow,Maxcol> & vec) -> decltype(vec.data())
{
	return vec.data()+vec.size();
}

template <typename Derived, int Nrow, int Ncol, int Options, int Maxrow, int Maxcol>
auto end(Eigen::Array<Derived,Nrow,Ncol,Options,Maxrow,Maxcol> & vec) -> decltype(vec.data())
{
	return vec.data()+vec.size();
}

} // namespace Eigen

#endif // Global namespace begin and end definitions

// Array typedefs
#if(1)

namespace IceBRG {

typedef Eigen::Array<flt_type,Eigen::Dynamic,1> flt_array_type;
typedef Eigen::Array<long_flt_type,Eigen::Dynamic,1> long_flt_array_type;

typedef Eigen::Array<int_type,Eigen::Dynamic,1> int_array_type;
typedef Eigen::Array<long_int_type,Eigen::Dynamic,1> long_int_array_type;
typedef Eigen::Array<short_int_type,Eigen::Dynamic,1> short_int_array_type;

typedef Eigen::Array<unsigned_int_type,Eigen::Dynamic,1> unsigned_int_array_type;
typedef Eigen::Array<long_unsigned_int_type,Eigen::Dynamic,1> long_unsigned_int_array_type;
typedef Eigen::Array<short_unsigned_int_type,Eigen::Dynamic,1> short_unsigned_int_array_type;

typedef Eigen::Array<complex_type,Eigen::Dynamic,1> complex_array_type;
typedef Eigen::Array<long_complex_type,Eigen::Dynamic,1> long_complex_array_type;

}
#endif // Array typedefs

// Serialization
#if(1)

namespace boost {
namespace serialization {

template<typename Archive, typename T_value>
void save(Archive & ar, const Eigen::Array<T_value,-1,-1,1,-1,-1> & v, const unsigned int version)
{
	int_type nrows(v.rows());
	int_type ncols(v.cols());

	ar << nrows << ncols;

	for(const auto & val : v)
	{
		ar << val;
	}
}
template<typename Archive, typename T_value>
void load(Archive & ar, Eigen::Array<T_value,-1,-1,1,-1,-1> & v, const unsigned int version)
{
	int_type nrows;
	int_type ncols;

	ar >> nrows >> ncols;

	v.resize(nrows,ncols);

	for(auto & val : v)
	{
		ar >> val;
	}
}

template<typename Archive, typename T_value>
inline void serialize(
        Archive & ar,
		Eigen::Array<T_value,-1,-1,1,-1,-1> & t,
        const unsigned int file_version
)
{
        split_free(ar, t, file_version);
}

} // namespace serialization
} // namespace boost

#endif

#endif // _BRG_EIGEN_HPP_
