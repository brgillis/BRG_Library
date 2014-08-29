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

#include "brg/global.h"

#include "brg/math/cache/cache_3d.hpp"
#include "brg/math/cache/cache_4d.hpp"

namespace brgastro {

class tNFW_sig_cache : public brg_cache_3d<tNFW_sig_cache>
{
	// Weak lensing signal delta-sigma from truncated NFW profile, using default c and tau
private:

	DECLARE_BRG_CACHE_3D_STATIC_VARS();

	friend class brg_cache_3d<tNFW_sig_cache>;

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
	const double _calculate( const double in_param_1, const double in_param_2,
			const double in_param_3 ) const;

public:

}; // class tNFW_sig_cache

class tNFW_offset_sig_cache : public brg_cache_4d<tNFW_offset_sig_cache>
{
	// Weak lensing signal delta-sigma from truncated NFW profile, using default c and tau,
	// offset from the centre of it by a distance R
private:

	DECLARE_BRG_CACHE_4D_STATIC_VARS();

	friend class brg_cache_4d<tNFW_offset_sig_cache>;

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
	const double _calculate( const double in_param_1, const double in_param_2,
			const double in_param_3, const double in_param_4 ) const;

public:

}; // class tNFW_offset_sig_cache

class tNFW_group_sig_cache : public brg_cache_4d<tNFW_group_sig_cache>
{
	// Weak lensing signal delta-sigma from truncated NFW profile, using default c and tau,
	// averaged over the positions of possible satellites of it
private:

	DECLARE_BRG_CACHE_4D_STATIC_VARS();

	friend class brg_cache_4d<tNFW_group_sig_cache>;

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
	const double _calculate( const double in_param_1, const double in_param_2,
			const double in_param_3, const double in_param_4 ) const;

public:

}; // class tNFW_group_sig_cache

class tNFW_shifted_sig_cache : public brg_cache_3d<tNFW_shifted_sig_cache>
{
	// Weak lensing signal delta-sigma from truncated NFW profile, using default c and tau,
	// with correction for shift due to lensing translation
private:

	DECLARE_BRG_CACHE_3D_STATIC_VARS();

	friend class brg_cache_3d<tNFW_shifted_sig_cache>;

protected:

	const std::string _name_base() const throw()
	{
		char name_base[BRG_CACHE_ND_NAME_SIZE] = "tN_s_sig";
		return name_base;
	}

	void _load_cache_dependencies() const
	{
		// This depends on offset_sig, so we'll have to load it through a dummy get
		brgastro::tNFW_offset_sig_cache().get(0,0,0,0,true);
	}

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0,-2,0,1,0,0,0); // Surface density
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function.
	const double _calculate( const double in_param_1, const double in_param_2,
			const double in_param_3 ) const;

public:

}; // class tNFW_shifted_sig_cache

} // end namespace brgastro

#endif /* _BRG_TNFW_CACHES_H_ */
