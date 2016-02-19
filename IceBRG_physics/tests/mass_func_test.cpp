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

#include "IceBRG_main/units/units.hpp"
#include "IceBRG_main/units/unit_conversions.hpp"
#include "IceBRG_physics/astro.h"

using namespace IceBRG;

BOOST_AUTO_TEST_SUITE (Mass_Func_Test)

BOOST_AUTO_TEST_CASE( mass_func_test )
{
	mass_type m1 = 1e10*unitconv::Msuntokg*kg;
	mass_type m2 = 1e12*unitconv::Msuntokg*kg;
	mass_type m3 = 1e14*unitconv::Msuntokg*kg;
	mass_type m4 = 1e16*unitconv::Msuntokg*kg;

	flt_t n1 = mass_function(m1);
	flt_t n2 = mass_function(m2);
	flt_t n3 = mass_function(m3);
	flt_t n4 = mass_function(m4);

	BOOST_CHECK_LE(n2,n1);
	BOOST_CHECK_LE(n3,n2);
	BOOST_CHECK_LE(n4,n3);

}

BOOST_AUTO_TEST_SUITE_END()


