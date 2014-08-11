/**       @file galaxy.cpp
 *
 *     Project: brg
 *        Path: /brg/sky_obj/galaxy.cpp
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#include "galaxy.h"

// brgastro::galaxy class methods
#if (1)
brgastro::galaxy::galaxy()
{
	clear();
}

const int brgastro::galaxy::clear()
{

	stellar_mass = 0;
	umag = gmag = rmag = imag = zmag = 0;
	umag_err = gmag_err = rmag_err = imag_err = zmag_err = 0;

	z_phot = z();
	z_phot_err = 0;
	odds = 1;
	phot_template = 0;

	host_group = 0;
	host_group_index = -1;

	return ( sky_obj::clear() );
}

#endif // end brgastro::galaxy class functions
