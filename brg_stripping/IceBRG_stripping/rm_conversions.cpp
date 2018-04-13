/**********************************************************************\
 @file rm_conversions.cpp
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "IceBRG_main/common.hpp"
#include "IceBRG_main/math/solvers/solvers.hpp"
#include "IceBRG_physics/constants.hpp"
#include "rm_conversions.hpp"

using namespace IceBRG;

flt_t IceBRG::get_delta_c(flt_t const & z)
{
    flt_t x = -1.0/(Omega_m/Omega_l*(1+z)*(1+z)*(1+z)+1.0);
    flt_t delta_c = 18*pi*pi+82*x-39*x*x;

    return delta_c;
}

Eigen::Array2f IceBRG::calculate_RM_ratio(flt_t const & z, flt_t const & c)
{
    /* compute the ratio R_a/R_b, so for instance to get R_200/R_100 set Delta_a = 200, Delta_b = 100
     * Delta_a and Delta_b can be relative to critical or background density; may need to provide Omega_M
     * need to provide a concentration parameter, c; MW is about 10-15?, cluster is about 4-10?
     * also compute M_a/M_b
     * should be valid for any redshift, but careful with some definitions of Rvir, Delta is a fn of z
     * */

    flt_t a_c = get_delta_c(z);

    flt_t b_c = 200; // Compare to R200

    flt_t a_b = a_c/b_c;

    auto f = [&] (flt_t const & ra_rb)
    {
        return square(std::pow(ra_rb,3) - (1.0/a_b) * (std::log(1+c)-c/(1+c))/(std::log(1+(1/ra_rb)*c)-c/(ra_rb+c)));
    };

    flt_t R_ratio = solve_sd(f,1);

    flt_t M_ratio = a_b * std::pow(R_ratio,3);

    Eigen::Array2f result;

    result << R_ratio, M_ratio;

    return result;
}

// Initialisation for caches
DEFINE_BRG_CACHE_2D( rratio_cache, flt_t, flt_t, flt_t,
        0, 5, 0.05, 3, 14, 0.1,
            return calculate_RM_ratio(in_param_1, in_param_2)[0];
        ,

        ,
)
DEFINE_BRG_CACHE_2D( mratio_cache, flt_t, flt_t, flt_t,
        0, 5, 0.05, 3, 14, 0.1,
            return calculate_RM_ratio(in_param_1, in_param_2)[1];
        ,

        ,
)

flt_t IceBRG::get_R_ratio(flt_t const & z, flt_t const & c)
{
    return rratio_cache().get(z,c);
}
flt_t IceBRG::get_M_ratio(flt_t const & z, flt_t const & c)
{
    return mratio_cache().get(z,c);
}

