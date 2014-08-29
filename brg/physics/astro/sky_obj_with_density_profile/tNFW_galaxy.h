/**       @file tNFW_galaxy.h
 *
 *     Project: brg
 *        Path: /brg/sky_obj_with_density_profile/tNFW_galaxy.h
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#ifndef _BRG_TNFW_GALAXY_H_
#define _BRG_TNFW_GALAXY_H_

#include "brg/global.h"

#include "brg/physics/astro/density_profile/density_profile.h"
#include "brg/physics/astro/density_profile/tNFW_profile.h"
#include "brg/physics/astro/sky_obj/sky_obj.h"
#include "brg/physics/astro/sky_obj/galaxy.h"

namespace brgastro {

class tNFW_galaxy: public tNFW_profile, public galaxy
{
	// Simple combination of the two classes
public:
	tNFW_galaxy()
	{
		galaxy();
		tNFW_profile();
	}
	virtual density_profile *density_profile_clone() const
	{
		return new tNFW_galaxy( *this );
	}
	virtual tNFW_profile *tNFW_profile_clone() const
	{
		return new tNFW_galaxy( *this );
	}
	virtual sky_obj *sky_obj_clone() const
	{
		return new tNFW_galaxy( *this );
	}
	virtual galaxy *galaxy_clone() const
	{
		return new tNFW_galaxy( *this );
	}
};

} // end namespace brgastro

#endif /* _BRG_TNFW_GALAXY_H_ */
