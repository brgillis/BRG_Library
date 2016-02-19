/**********************************************************************\
 @file dfa_cache_test.cpp
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
#include "IceBRG_physics/astro_caches.h"

using namespace IceBRG;

BOOST_AUTO_TEST_SUITE (DFA_Cache_Test)

BOOST_AUTO_TEST_CASE( dfa_cache_test )
{
    dfa_cache cache;

    cache.set_range(0,1,0.1);

    auto dfa_1 = cache.get(0.05);
    auto dfa_2 = cache.get(0.15);

    BOOST_CHECK_GT(value_of(dfa_2),value_of(dfa_1));

    // Allow 10% closeness here, due to coarseness of cache in this test
    BOOST_CHECK_CLOSE(value_of(dfa_1),value_of(201.7*unitconv::Mpctom*m),10);
    BOOST_CHECK_CLOSE(value_of(dfa_2),value_of(540.0*unitconv::Mpctom*m),10);

    auto inv_dfa_1 = cache.inverse_get(dfa_1);
    auto inv_dfa_2 = cache.inverse_get(dfa_2);

    BOOST_CHECK_CLOSE(value_of(inv_dfa_1),0.05,0.1);
    BOOST_CHECK_CLOSE(value_of(inv_dfa_2),0.15,0.1);
}

BOOST_AUTO_TEST_SUITE_END()

