/**********************************************************************\
 @file rm_conversions.hpp
 ------------------

 TODO <Insert file description here>

 **********************************************************************

 Copyright (C) 2018 brg

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

#ifndef ICEBRG_STRIPPING_RM_CONVERSIONS_HPP_
#define ICEBRG_STRIPPING_RM_CONVERSIONS_HPP_

#include <tuple>

#include <Eigen/Core>
#include "IceBRG_main/common.hpp"
#include "IceBRG_main/math/cache/cache_2d.hpp"

namespace IceBRG
{

DECLARE_BRG_CACHE_2D(rratio_cache,r_ratio,flt_t,flt_t,flt_t);
DECLARE_BRG_CACHE_2D(mratio_cache,m_ratio,flt_t,flt_t,flt_t);

flt_t get_delta_c(flt_t const & z);

Eigen::Array2f calculate_RM_ratio(flt_t const & z, flt_t const & c);

flt_t get_R_ratio(flt_t const & z, flt_t const & c);
flt_t get_M_ratio(flt_t const & z, flt_t const & c);

} // namespace IceBRG


#endif // ICEBRG_STRIPPING_RM_CONVERSIONS_HPP_
