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
#include "lensing_tNFW_profile.h"

#include "lensing_tNFW_caches.h"

// Initialisation for brgastro::tNFW_sig_cache
double tNFW_sig_cache_mins[] = {std::log(1e7*unitconv::Msuntokg),
								0.1,
								std::log(0.1*unitconv::kpctom)};
double tNFW_sig_cache_maxes[] = {std::log(1e16*unitconv::Msuntokg),
								 1.5,
								 std::log(20000*unitconv::kpctom)};
double tNFW_sig_cache_steps[] = {(std::log(1e16)-std::log(1e7))/1000,
								 0.1,
								 (std::log(20000)-std::log(0.1))/1000};
DEFINE_BRG_CACHE_ND_STATIC_VARS( tNFW_sig_cache, tNFW_sig_cache_mins,
		tNFW_sig_cache_maxes, tNFW_sig_cache_steps, 3 );

// Initialisation for brgastro::tNFW_offset_sig_cache
double tNFW_offset_sig_cache_mins[] = {std::log(1e10*unitconv::Msuntokg),
									   0.1,
									   std::log(0.1*unitconv::kpctom),
									   std::log(0.1*unitconv::kpctom)};
double tNFW_offset_sig_cache_maxes[] = {std::log(1e16*unitconv::Msuntokg),
										1.5,
										std::log(20000*unitconv::kpctom),
										std::log(4000*unitconv::kpctom)};
double tNFW_offset_sig_cache_steps[] = {(std::log(1e16)-std::log(1e7))/10,
										0.2,
										(std::log(20000)-std::log(0.1))/10,
										(std::log(4000)-std::log(0.1))/10};
DEFINE_BRG_CACHE_ND_STATIC_VARS( tNFW_offset_sig_cache, tNFW_offset_sig_cache_mins, tNFW_offset_sig_cache_maxes,
		tNFW_offset_sig_cache_steps, 4 );

// Initialisation for brgastro::tNFW_group_sig_cache
double tNFW_group_sig_cache_mins[] =  {std::log(1e10*unitconv::Msuntokg),
									   0.1,
									   std::log(0.1*unitconv::kpctom),
									   1};
double tNFW_group_sig_cache_maxes[] = {std::log(1e16*unitconv::Msuntokg),
									   1.5,
									   std::log(20000*unitconv::kpctom),
									   20};
double tNFW_group_sig_cache_steps[] = {(std::log(1e16)-std::log(1e7))/10,
									   0.1,
									   (std::log(20000)-std::log(0.1))/10,
									   0.2};
DEFINE_BRG_CACHE_ND_STATIC_VARS( tNFW_group_sig_cache, tNFW_group_sig_cache_mins, tNFW_group_sig_cache_maxes,
		tNFW_group_sig_cache_steps, 4 );

// brgastro::tNFW_sig_cache class methods
#if (1)
const int brgastro::tNFW_sig_cache::_calculate( const brgastro::vector<double> in_params,
		double & out_params ) const
{
	try
	{
		const double mass = std::exp(in_params.at(0));
		const double z = in_params.at(1);
		const double r = std::exp(in_params.at(2));
		brgastro::lensing_tNFW_profile profile(mass,z);
		out_params = profile.deltasigma( r );
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

const int brgastro::tNFW_offset_sig_cache::_calculate( const brgastro::vector<double> in_params,
		double & out_params ) const
{
	try
	{
		const double mass = std::exp(in_params.at(0));
		const double z = in_params.at(1);
		const double r = std::exp(in_params.at(2));
		const double offset_r = std::exp(in_params.at(3));
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

const int brgastro::tNFW_group_sig_cache::_calculate( const brgastro::vector<double> in_params,
		double & out_params ) const
{
	try
	{
		const double mass = std::exp(in_params.at(0));
		const double z = in_params.at(1);
		const double r = std::exp(in_params.at(2));
		const double group_c = in_params.at(3);
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


