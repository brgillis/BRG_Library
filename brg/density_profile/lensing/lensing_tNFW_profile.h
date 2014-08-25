/**       @file lensing_tNFW_profile.h
 *
 *     Project: brg
 *        Path: /brg/density_profile/lensing_density_profile/lensing_tNFW_profile.h
 *
 *  Created on: 12 Aug 2014
 *      Author: brg
 */

#ifndef _BRG_LENSING_TNFW_PROFILE_H_
#define _BRG_LENSING_TNFW_PROFILE_H_

#include "../../brg_global.h"

#include "../tNFW_profile.h"
#include "lensing_profile_extension.h"

namespace brgastro {

/*
 *
 */
class lensing_tNFW_profile: public tNFW_profile, public lensing_profile_extension {
public:
	// Constructors and destructor
#if (1)
	lensing_tNFW_profile()
	{
	}
	lensing_tNFW_profile( const BRG_MASS &init_mvir0, const double init_z,
			const double init_c = -1, const double init_tau = -1 )
	: tNFW_profile(init_mvir0, init_z, init_c, init_tau)
	{
	}
	virtual ~lensing_tNFW_profile()
	{
	}
#endif // Constructors and destructor

	// Lensing related methods
#if (1)
	const BRG_UNITS proj_dens( const BRG_DISTANCE &R,
			const bool silent = false ) const;
	const BRG_UNITS proj_enc_dens( const BRG_DISTANCE &R, const bool silent =
			false ) const;
	const BRG_MASS proj_enc_mass( const BRG_DISTANCE &R, const bool silent =
			false ) const;
	const BRG_UNITS quick_WLsig( const BRG_DISTANCE &R, const bool silent =
			false ) const;
	const BRG_UNITS quick_offset_WLsig( const BRG_DISTANCE &R,
			const BRG_DISTANCE &offset_R, const bool silent = false ) const;
	const BRG_UNITS quick_group_WLsig( const BRG_DISTANCE &R,
			const double group_c, const bool silent = false ) const;
	const BRG_UNITS quick_shifted_WLsig( const BRG_DISTANCE &R, const bool silent = false ) const;

#endif // Lensing related methods

	// Implementations of pure virtual parts of lensing_profile_extension
#if (1)
	IMPLEMENT_BRG_LENSING_EXTENSION_METHODS(tNFW_profile);
#endif

	// Implementations of clone functions
#if (1)
	virtual density_profile *density_profile_clone() const
	{
		return new lensing_tNFW_profile( *this );
	}
	virtual tNFW_profile *tNFW_profile_clone() const
	{
		return new lensing_tNFW_profile( *this );
	}
	virtual lensing_profile_extension *lensing_profile_extension_clone() const
	{
		return new lensing_tNFW_profile( *this );
	}
#endif
};

} // end namespace brgastro

#endif /* _BRG_LENSING_TNFW_PROFILE_H_ */
