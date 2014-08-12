/**       @file point_mass_profile.h
 *
 *     Project: brg
 *        Path: /brg/density_profile/point_mass_profile.h
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#ifndef _BRG_POINT_MASS_PROFILE_H_
#define _BRG_POINT_MASS_PROFILE_H_

#include <vector>

#include "../brg_global.h"

#include "../brg_units.h"
#include "density_profile.h"

namespace brgastro {

class point_mass_profile: public density_profile
{
	/**********************************
	 point_mass_profile class
	 ------------------------

	 A virtual class for a point mass

	 Defined by two parameters:
	 mass and z

	 Parent class:
	 density profile

	 Derived classes:
	 (none)

	 **********************************/
private:
#if (1) // private member variables specific to this class

	BRG_MASS _mass_;

#endif // end private member variables

public:
#if (1) // public variables and functions

#if (1) // Constructors
	point_mass_profile();

	point_mass_profile( const BRG_MASS init_mass, const double init_z );

#endif // End constructors

	// Destructor
	~point_mass_profile();

#if (1) // Set functions
	virtual const int set_mvir( const BRG_MASS &new_halo_mass, bool silent =
			false );
	virtual const int set_parameters( const unsigned int num_parameters,
			const std::vector< BRG_UNITS > &new_parameters,
			bool silent = false );
#endif // End set functions

#if (1) // Basic get functions
	const BRG_MASS mass() const;

	const BRG_MASS hmvir() const;
	const BRG_MASS mvir() const;
	const BRG_MASS mtot() const;

	const BRG_DISTANCE rvir() const;
	const BRG_DISTANCE rt(const bool silent = false) const;
	const BRG_DISTANCE rs() const;

	const BRG_VELOCITY vvir() const;

#endif // end basic get functions

#if (1) // advanced get functions
	const BRG_UNITS dens( const BRG_DISTANCE &r ) const;
	const BRG_UNITS enc_dens( const BRG_DISTANCE &r,
			const bool silent = false ) const;
	const BRG_MASS enc_mass( const BRG_DISTANCE &r, const bool silent =
				true ) const; // Mass enclosed with sphere of radius r
	const unsigned int num_parameters() const
	{
		return 2; // Mass and redshift
	}
	const int get_parameters( std::vector< BRG_UNITS > & parameters,
			const bool silent = false ) const;

	const int get_parameter_names( std::vector< std::string > & parameter_names,
			const bool silent = false ) const;
#endif // end advanced get functions

#if (1) // Other operations

	const int truncate_to_fraction( const double fraction,
			bool silent = false );
	virtual density_profile *density_profile_clone() const
	{
		return new point_mass_profile( *this );
	}
	virtual point_mass_profile *point_mass_profile_clone() const
	{
		return new point_mass_profile( *this );
	}

#endif

#endif // end public variables and functions
};
// class point_mass_profile

} // end namespace brgastro

#endif /* _BRG_POINT_MASS_PROFILE_H_ */
