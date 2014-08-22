/**       @file tNFW_caches.cpp
 *
 *     Project: brg
 *        Path: /brg/density_profile/tNFW_caches.cpp
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../../brg_global.h"

#include "../../brg_astro.h"
#include "../../brg_units.h"

#include "../../cache/cache_4d.hpp"

#include "lensing_tNFW_profile.h"

#include "lensing_tNFW_caches.h"

// Initialisation for brgastro::tNFW_sig_cache
DEFINE_BRG_CACHE_3D_STATIC_VARS( tNFW_sig_cache,
	std::log(1e7*unitconv::Msuntokg),std::log(1e16*unitconv::Msuntokg),(std::log(1e16)-std::log(1e7))/100,
	0.1,1.5,0.1,
	std::log(10*unitconv::kpctom),std::log(20000*unitconv::kpctom),(std::log(20000)-std::log(0.1))/100);

// Initialisation for brgastro::tNFW_offset_sig_cache
DEFINE_BRG_CACHE_4D_STATIC_VARS( tNFW_offset_sig_cache,
	std::log(1e10*unitconv::Msuntokg),std::log(1e16*unitconv::Msuntokg),(std::log(1e16)-std::log(1e10))/100,
	0.1,1.5,0.2,
	std::log(10*unitconv::kpctom),std::log(20000*unitconv::kpctom),(std::log(20000)-std::log(0.1))/100,
	std::log(0.1*unitconv::kpctom),std::log(4000*unitconv::kpctom),(std::log(4000)-std::log(0.1))/100);

// Initialisation for brgastro::tNFW_group_sig_cache
DEFINE_BRG_CACHE_4D_STATIC_VARS( tNFW_group_sig_cache,
	std::log(1e10*unitconv::Msuntokg),std::log(1e16*unitconv::Msuntokg),(std::log(1e16)-std::log(1e10))/100,
	0.1,1.5,0.2,
	std::log(0.1*unitconv::kpctom),std::log(20000*unitconv::kpctom),(std::log(20000)-std::log(0.1))/100,
	2.5,10.,0.5);

// brgastro::tNFW_sig_cache class methods
#if (1)
const int brgastro::tNFW_sig_cache::_calculate( const double in_param_1, const double in_param_2,
		const double in_param_3, double & out_params ) const
{
	try
	{
		const double mass = std::exp(in_param_1);
		const double z = in_param_2;
		const double r = std::exp(in_param_3);
		brgastro::lensing_tNFW_profile profile(mass,z);
		out_params = profile.WLsig( r );
	}
	catch(const std::out_of_range &e)
	{
		return INVALID_ARGUMENTS_ERROR;
	}
	return 0;
}

#endif // end brgastro::tNFW_sig_cache functions

// brgastro::NFW_offset_sig_cache class methods
#if (1)

const int brgastro::tNFW_offset_sig_cache::_calculate( const double in_param_1, const double in_param_2,
		const double in_param_3, const double in_param_4, double & out_params ) const
{
	try
	{
		const double mass = std::exp(in_param_1);
		const double z = in_param_2;
		const double r = std::exp(in_param_3);
		const double offset_r = std::exp(in_param_4);
		brgastro::lensing_tNFW_profile profile(mass,z);
		out_params = profile.offset_WLsig( r, offset_r );
	}
	catch(const std::out_of_range &e)
	{
		return INVALID_ARGUMENTS_ERROR;
	}
	return 0;
}

#endif // end brgastro::NFW_offset_sig_cache functions

// brgastro::tNFW_group_sig_cache class methods
#if (1)

const int brgastro::tNFW_group_sig_cache::_calculate( const double in_param_1, const double in_param_2,
		const double in_param_3, const double in_param_4, double & out_params ) const
{
	try
	{
		const double mass = std::exp(in_param_1);
		const double z = in_param_2;
		const double r = std::exp(in_param_3);
		const double group_c = in_param_4;
		brgastro::lensing_tNFW_profile profile(mass,z);
		out_params = profile.semiquick_group_WLsig( r, group_c );
	}
	catch(const std::out_of_range &e)
	{
		return INVALID_ARGUMENTS_ERROR;
	}
	return 0;
}

#endif // end brgastro::tNFW_group_sig_cache functions


