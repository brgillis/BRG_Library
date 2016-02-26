/**********************************************************************\
  @file galaxy.cpp

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

#include "galaxy.hpp"

#include "IceBRG_main/common.hpp"


// IceBRG::galaxy class methods
#if (1)
IceBRG::galaxy::galaxy( const angle_type & init_ra, const angle_type & init_dec, flt_t init_z,
		const angle_type & init_ra_err, const angle_type & init_dec_err, flt_t init_z_err,
		const mass_type & init_stellar_mass, flt_t init_mag, flt_t init_mag_err )
:	sky_obj(init_ra,init_dec,init_z,init_ra_err,init_dec_err,init_z_err),
 	stellar_mass(init_stellar_mass),
 	umag(init_mag), umag_err(init_mag_err),
	gmag(init_mag), gmag_err(init_mag_err),
 	rmag(init_mag), rmag_err(init_mag_err),
 	imag(init_mag), imag_err(init_mag_err),
 	zmag(init_mag), zmag_err(init_mag_err),
 	z_phot(z()), z_phot_err(0),
 	odds(1), phot_template(0),
 	host_group(NULL), host_group_index(-1)
{

}

void IceBRG::galaxy::clear()
{

	sky_obj::clear();

	stellar_mass = 0;
	umag = gmag = rmag = imag = zmag = 0;
	umag_err = gmag_err = rmag_err = imag_err = zmag_err = 0;

	z_phot = z();
	z_phot_err = 0;
	odds = 1;
	phot_template = 0;

	host_group = 0;
	host_group_index = -1;
}

#endif // end IceBRG::galaxy class functions
