/**********************************************************************\
 @file mass_func_test.cpp
 ------------------

 TODO <Insert file description here>

 **********************************************************************

 Copyright (C) 2016 brg

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

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "IceBRG_main/math/calculus/integrate.hpp"
#include "IceBRG_main/units/units.hpp"
#include "IceBRG_main/units/unit_conversions.hpp"
#include "IceBRG_physics/astro.h"

using namespace IceBRG;

BOOST_AUTO_TEST_SUITE (Mass_Func_Test)

BOOST_AUTO_TEST_CASE( mass_func_test )
{
	flt_t z = 0;

	volume_type Gpc_cubed = cube(1e9*unitconv::pctom*m);

	flt_t log10msun_m1 = 8;
	flt_t log10msun_m2 = 10;
	flt_t log10msun_m3 = 12;
	flt_t log10msun_m4 = 14;
	flt_t log10msun_m5 = 16;

	// Try for z=0

	flt_t n1 = log10_mass_function(log10msun_m1,z)*Gpc_cubed;
	flt_t n2 = log10_mass_function(log10msun_m2,z)*Gpc_cubed;
	flt_t n3 = log10_mass_function(log10msun_m3,z)*Gpc_cubed;
	flt_t n4 = log10_mass_function(log10msun_m4,z)*Gpc_cubed;
	flt_t n5 = log10_mass_function(log10msun_m5,z)*Gpc_cubed;

	BOOST_CHECK_LE(value_of(n2),value_of(n1));
	BOOST_CHECK_LE(value_of(n3),value_of(n2));
	BOOST_CHECK_LE(value_of(n4),value_of(n3));
	BOOST_CHECK_LE(value_of(n5),value_of(n4));

	auto density_at_scale_0 = [&] (flt_t const & log10_mass)
	{
		return log10_mass_function(log10_mass,z)*std::pow(10.,log10_mass)*unitconv::Msuntokg*kg;
	};

	density_type dens = integrate_Romberg(density_at_scale_0,0.,16.);

	density_type matter_dens = rho_bar(z);

	BOOST_CHECK_CLOSE(value_of(dens),value_of(matter_dens),10);

	// Try for z=0.9

	z = 0.9;

	n1 = log10_mass_function(log10msun_m1,z)*Gpc_cubed;
	n2 = log10_mass_function(log10msun_m2,z)*Gpc_cubed;
	n3 = log10_mass_function(log10msun_m3,z)*Gpc_cubed;
	n4 = log10_mass_function(log10msun_m4,z)*Gpc_cubed;
	n5 = log10_mass_function(log10msun_m5,z)*Gpc_cubed;

	BOOST_CHECK_LE(value_of(n2),value_of(n1));
	BOOST_CHECK_LE(value_of(n3),value_of(n2));
	BOOST_CHECK_LE(value_of(n4),value_of(n3));
	BOOST_CHECK_LE(value_of(n5),value_of(n4));

	auto density_at_scale_0_9 = [&] (flt_t const & log10_mass)
	{
		return log10_mass_function(log10_mass,z)*std::pow(10.,log10_mass)*unitconv::Msuntokg*kg;
	};

	dens = integrate_Romberg(density_at_scale_0_9,0.,16.);

	matter_dens = rho_bar(z);

	BOOST_CHECK_CLOSE(value_of(dens),value_of(matter_dens),10);
}

BOOST_AUTO_TEST_SUITE_END()


