/**       @file tNFW_caches.h
 *
 *     Project: brg
 *        Path: /brg/density_profile/tNFW_caches.h
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#ifndef _BRG_TNFW_CACHES_H_
#define _BRG_TNFW_CACHES_H_

#include <string>
#include <vector>

#include "../../brg_global.h"

namespace brgastro {

class tNFW_sig_cache : public brg_cache_nd<tNFW_sig_cache>
{
	// Weak lensing signal delta-sigma from truncated NFW profile, using default c and tau
private:

	DECLARE_BRG_CACHE_ND_STATIC_VARS();

	friend class brg_cache_nd<tNFW_sig_cache>;

protected:

	const std::string _name_base() const throw()
	{
		char name_base[BRG_CACHE_ND_NAME_SIZE] = "tNFW_sig";
		return name_base;
	}

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0,-2,0,1,0,0,0); // Surface density
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function.
	const int _calculate( const brgastro::vector<double> in_params, double & out_params ) const;

public:

}; // class tNFW_sig_cache

class tNFW_offset_sig_cache : public brg_cache_nd<tNFW_offset_sig_cache>
{
	// Weak lensing signal delta-sigma from truncated NFW profile, using default c and tau,
	// offset from the centre of it by a distance R
private:

	DECLARE_BRG_CACHE_ND_STATIC_VARS();

	friend class brg_cache_nd<tNFW_offset_sig_cache>;

protected:

	const std::string _name_base() const throw()
	{
		char name_base[BRG_CACHE_ND_NAME_SIZE] = "tN_o_sig";
		return name_base;
	}

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0,-2,0,1,0,0,0); // Surface density
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function.
	const int _calculate( const brgastro::vector<double> in_params, double & out_params ) const;

public:

}; // class tNFW_offset_sig_cache

class tNFW_group_sig_cache : public brg_cache_nd<tNFW_group_sig_cache>
{
	// Weak lensing signal delta-sigma from truncated NFW profile, using default c and tau,
	// averaged over the positions of possible satellites of it
private:

	DECLARE_BRG_CACHE_ND_STATIC_VARS();

	friend class brg_cache_nd<tNFW_group_sig_cache>;

protected:

	const std::string _name_base() const throw()
	{
		char name_base[BRG_CACHE_ND_NAME_SIZE] = "tN_g_sig";
		return name_base;
	}

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0,-2,0,1,0,0,0); // Surface density
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function.
	const int _calculate( const brgastro::vector<double> in_params, double & out_params ) const;

public:

}; // class tNFW_group_sig_cache

} // end namespace brgastro

#endif /* _BRG_TNFW_CACHES_H_ */
