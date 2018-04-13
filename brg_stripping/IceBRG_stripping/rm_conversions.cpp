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
#include "IceBRG_stripping/rm_conversions.hpp"

using namespace IceBRG;

flt_t get_delta_c(flt_t const & z)
{
    flt_t x = -1.0/(Omega_m/Omega_l*(1+z)*(1+z)*(1+z)+1.0);
    flt_t delta_c = 18*pi*pi+82*x-39*x*x;

    return delta_c;
}

std::tuple<flt_t,flt_t> calculate_RM_ratio(flt_t const & z, flt_t const & c)
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

    return std::tie(R_ratio,M_ratio);
}

