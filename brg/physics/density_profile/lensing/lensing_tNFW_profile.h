/**********************************************************************\
  @file lensing_tNFW_profile.h

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

// body file: lensing_tNFW_profile.cpp

#ifndef _BRG_LENSING_TNFW_PROFILE_H_
#define _BRG_LENSING_TNFW_PROFILE_H_

#include "brg/global.h"

#include "brg/physics/density_profile/lensing/lensing_profile_extension.hpp"
#include "brg/physics/density_profile/tNFW_profile.h"

namespace brgastro {

/*
 *
 */
class lensing_tNFW_profile: public tNFW_profile, public lensing_profile_extension<lensing_tNFW_profile> {
public:
	// Constructors and destructor
#if (1)
	lensing_tNFW_profile()
	{
	}
	lensing_tNFW_profile( CONST_BRG_MASS_REF init_mvir0, const double init_z,
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
	const BRG_UNITS proj_dens( CONST_BRG_DISTANCE_REF R,
			const bool silent = false ) const;
	const BRG_UNITS proj_enc_dens( CONST_BRG_DISTANCE_REF R, const bool silent =
			false ) const;
	const BRG_MASS proj_enc_mass( CONST_BRG_DISTANCE_REF R, const bool silent =
			false ) const;
	const BRG_UNITS quick_WLsig( CONST_BRG_DISTANCE_REF R, const bool silent =
			false ) const;
	const BRG_UNITS quick_offset_WLsig( CONST_BRG_DISTANCE_REF R,
			CONST_BRG_DISTANCE_REF offset_R, const bool silent = false ) const;
	const BRG_UNITS quick_group_WLsig( CONST_BRG_DISTANCE_REF R,
			const double group_c, const bool silent = false ) const;
	const BRG_UNITS quick_shifted_WLsig( CONST_BRG_DISTANCE_REF R, const bool silent = false ) const;

#endif // Lensing related methods

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
#endif
};

} // end namespace brgastro

#endif /* _BRG_LENSING_TNFW_PROFILE_H_ */
