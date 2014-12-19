/**********************************************************************\
  @file lensing_tNFW_caches.cpp

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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "brg/global.h"

#include "brg/math/cache/cache_4d.hpp"

#include "brg_lensing/lensing_tNFW_profile.h"

#include "brg_physics/astro.h"
#include "brg_physics/units/unit_obj.h"

#include "lensing_tNFW_caches.h"

// Initialisation for brgastro::tNFW_sig_cache
DEFINE_BRG_CACHE_3D_STATIC_VARS( tNFW_sig_cache,
	std::log(1e7*unitconv::Msuntokg),std::log(1e16*unitconv::Msuntokg),(std::log(1e16)-std::log(1e7))/300,
	0.1,1.5,0.1,
	std::log(0.1*unitconv::kpctom),std::log(20000*unitconv::kpctom),(std::log(20000)-std::log(0.1))/300);

// Initialisation for brgastro::tNFW_offset_sig_cache
DEFINE_BRG_CACHE_4D_STATIC_VARS( tNFW_offset_sig_cache,
	std::log(1e10*unitconv::Msuntokg),std::log(1e16*unitconv::Msuntokg),(std::log(1e16)-std::log(1e10))/100,
	0.1,1.5,0.4,
	std::log(0.1*unitconv::kpctom),std::log(20000*unitconv::kpctom),(std::log(20000)-std::log(0.1))/100,
	std::log(0.01*unitconv::kpctom),std::log(4000*unitconv::kpctom),(std::log(4000)-std::log(0.01))/300);

// Initialisation for brgastro::tNFW_group_sig_cache
DEFINE_BRG_CACHE_4D_STATIC_VARS( tNFW_group_sig_cache,
	std::log(1e10*unitconv::Msuntokg),std::log(1e16*unitconv::Msuntokg),(std::log(1e16)-std::log(1e10))/100,
	0.1,1.5,0.4,
	std::log(0.1*unitconv::kpctom),std::log(20000*unitconv::kpctom),(std::log(20000)-std::log(0.1))/100,
	2.5,10.,0.5);

// Initialisation for brgastro::tNFW_shifted_sig_cache
DEFINE_BRG_CACHE_3D_STATIC_VARS( tNFW_shifted_sig_cache,
	std::log(1e7*unitconv::Msuntokg),std::log(1e16*unitconv::Msuntokg),(std::log(1e16)-std::log(1e7))/100,
	0.1,1.5,0.4,
	std::log(0.1*unitconv::kpctom),std::log(20000*unitconv::kpctom),(std::log(20000)-std::log(0.1))/100);

// brgastro::tNFW_sig_cache class methods
#if (1)
const double brgastro::tNFW_sig_cache::_calculate( const double in_param_1, const double in_param_2,
		const double in_param_3 ) const
{
	const double mass = std::exp(in_param_1);
	const double z = in_param_2;
	const double r = std::exp(in_param_3);
	brgastro::lensing_tNFW_profile profile(mass,z);
	return profile.WLsig( r );
}

#endif // end brgastro::tNFW_sig_cache methods

// brgastro::NFW_offset_sig_cache class methods
#if (1)

const double brgastro::tNFW_offset_sig_cache::_calculate( const double in_param_1, const double in_param_2,
		const double in_param_3, const double in_param_4 ) const
{
	const double mass = std::exp(in_param_1);
	const double z = in_param_2;
	const double r = std::exp(in_param_3);
	const double offset_r = std::exp(in_param_4);
	brgastro::lensing_tNFW_profile profile(mass,z);
	return profile.offset_WLsig( r, offset_r );
}

#endif // end brgastro::NFW_offset_sig_cache functions

// brgastro::tNFW_group_sig_cache class methods
#if (1)

const double brgastro::tNFW_group_sig_cache::_calculate( const double in_param_1, const double in_param_2,
		const double in_param_3, const double in_param_4 ) const
{
	const double mass = std::exp(in_param_1);
	const double z = in_param_2;
	const double r = std::exp(in_param_3);
	const double group_c = in_param_4;
	brgastro::lensing_tNFW_profile profile(mass,z);
	return profile.semiquick_group_WLsig( r, group_c );
}

#endif // end brgastro::tNFW_group_sig_cache functions

// brgastro::tNFW_shifted_sig_cache class methods
#if (1)
const double brgastro::tNFW_shifted_sig_cache::_calculate( const double in_param_1, const double in_param_2,
		const double in_param_3 ) const
{
	const double mass = std::exp(in_param_1);
	const double z = in_param_2;
	const double R = std::exp(in_param_3);
	brgastro::lensing_tNFW_profile profile(mass,z);
	return profile.semiquick_shifted_WLsig( R );
}

#endif // end brgastro::tNFW_sig_cache methods


