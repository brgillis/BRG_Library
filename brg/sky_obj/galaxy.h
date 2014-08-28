/**       @file galaxy.h
 *
 *     Project: brg
 *        Path: /brg/sky_obj/galaxy.h
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#ifndef _BRG_GALAXY_H_
#define _BRG_GALAXY_H_

#include "sky_obj.h"

namespace brgastro {

// Forward declare galaxy_group class so we can point to it
class galaxy_group;

class galaxy: public sky_obj
{
	/**********************************
	 galaxy class
	 -------------

	 A class for galaxies.

	 Parent class: sky_obj

	 **********************************/

private:
	// private member variables (none needed in current implementation. Should make all private for consistency in future though)

public:

	// Public member variables

	BRG_MASS stellar_mass;
	double umag, gmag, rmag, imag, zmag;
	double umag_err, gmag_err, rmag_err, imag_err, zmag_err;

	double z_phot, z_phot_err;
	double odds;
	double phot_template;

	galaxy_group *host_group; // Pointer to group this galaxy resides within
	int host_group_index;

	// Public member functions

	// Constructor
	galaxy();

	// Copy constructor
	//galaxy(const galaxy other_galaxy); // Implicit is fine for us

	// Virtual destructor
	virtual ~galaxy()
	{
	}

	// Clear function
	virtual void clear();

	// Clone functions
	virtual redshift_obj *redshift_obj_clone()
	{
		return new galaxy( *this );
	}
	virtual sky_obj *sky_obj_clone()
	{
		return new galaxy( *this );
	}
	virtual galaxy *galaxy_clone()
	{
		return new galaxy( *this );
	}

};
// class galaxy

} // end namespace brgastro

#endif /* _BRG_GALAXY_H_ */
