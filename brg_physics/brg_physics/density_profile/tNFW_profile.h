/**********************************************************************\
  @file tNFW_profile.h

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

// body file: brg_physics/astro/density_profile/tNFW_profile.cpp

#ifndef _BRG_TNFW_PROFILE_H_
#define _BRG_TNFW_PROFILE_H_

#include <cmath>
#include <cstdlib>
#include <vector>

#include "brg/common.h"

#include "brg/utility.hpp"
#include "brg_physics/units/unit_obj.h"
#include "brg_physics/density_profile/density_profile.h"

namespace brgastro {

class tNFW_profile: public density_profile
{
	/**********************************
	 tNFW_profile class
	 -------------

	 A virtual class for an object whose
	 mass fits a truncated NFW profile.
	 See ... .

	 Defined by four parameters:
	 mvir0, z, c, and tau

	 Parent class:
	 density profile

	 Derived classes:
	 tNFW_galaxy

	 **********************************/
private:
#if (1) // Private member variables specific to this class
	BRG_MASS _mvir0_;
	mutable BRG_DISTANCE _rvir_cache_;
	mutable bool _rvir_cached_;
	flt_type _c_, _tau_;

#endif // end private member variables

protected:

	void _uncache_mass();

	flt_type _taufm( const flt_type mratio, flt_type precision = 0.00001,
			const bool silent = false ) const; // tau from Mtot/Mvir
	flt_type _delta_c() const // Simple function of concentration used as a step in calculating NFW densities
	{
		return ( 200. / 3. ) * cube( _c_ )
				/ ( std::log( 1 + _c_ ) - _c_ / ( 1 + _c_ ) );
	}

	// Functions relating to tNFW profiles
	flt_type _cfm( const BRG_MASS mass, const flt_type z=0 ) const // Concentration from mass relationship, from Neto
	{
		return 4.67 * std::pow( mass * unitconv::kgtottMsun / ( 1e4 ), -0.11 );
	}
	flt_type _cfm() const // Concentration from mass relationship, from Neto
	{
		return _cfm(_mvir0_,z());
	}
	flt_type _mftau( const flt_type tau, const flt_type conc ) const // Mtot/Mvir from tau
	{
		if(tau<=0) return 0;
		flt_type M0oM200 = 1 / ( std::log( 1 + conc ) - conc / ( 1 + conc ) );
		flt_type tautau = tau*tau;
		flt_type result =
				M0oM200 * tautau / square( tautau + 1 )
						* ( ( tautau - 1 ) * std::log( tau ) + tau * pi
								- ( tautau + 1 ) );
		return result;
	}
	flt_type _mftau( const flt_type tau ) const // Mtot/Mvir from tau
	{
		return _mftau(tau,_c_);
	}
	flt_type _mftau() const // Mtot/Mvir from tau
	{
		return _mftau(_tau_,_c_);
	}
public:
#if (1) // public variables and functions

#if (1) // Constructors
	tNFW_profile();

	explicit tNFW_profile( const BRG_MASS & init_mvir0, const flt_type init_z = 0,
			const flt_type init_c = -1, const flt_type init_tau = -1 );

#endif // End constructors

	// Destructor
	~tNFW_profile();

#if (1) // Set functions

	void set_mvir( const BRG_MASS & new_halo_mass, const bool silent =
			false );
	void set_parameters( const std::vector< BRG_UNITS > & new_parameters,
			const bool silent = false );

	void set_z( const flt_type new_z );
	void set_tau( const flt_type new_halo_tau, const bool silent = false );
	void set_c( const flt_type new_halo_c, const bool silent = false );

#endif // End set functions

#if (1) //Basic get functions

	BRG_MASS mvir() const;
	const BRG_MASS & mvir0() const;

	BRG_MASS mtot() const;

	BRG_DISTANCE rvir() const;
	BRG_DISTANCE rvir0() const;
	BRG_DISTANCE rt( const bool silent = false ) const;
	BRG_DISTANCE rs() const;

	BRG_VELOCITY vvir() const;
	BRG_VELOCITY vvir0() const;

	flt_type c() const;
	flt_type tau() const;
#endif // end basic get functions

#if (1) // advanced get functions

	BRG_UNITS dens( const BRG_DISTANCE & r ) const;
	BRG_MASS enc_mass( const BRG_DISTANCE & r,
			const bool silent = false ) const;
	size_t num_parameters() const
	{
		return 4; // Mass, redshift, c, and tau
	}
	std::vector< BRG_UNITS > get_parameters( const bool silent = true ) const;

	std::vector< std::string > get_parameter_names( const bool silent = true ) const;
#endif // end advanced get functions

#if (1) // Other operations

	virtual void truncate_to_fraction( const flt_type fraction,
			const bool silent = false );
	virtual redshift_obj *redshift_obj_clone() const
	{
		return new tNFW_profile( *this );
	}
	virtual density_profile *density_profile_clone() const
	{
		return new tNFW_profile( *this );
	}
	virtual tNFW_profile *tNFW_profile_clone() const
	{
		return new tNFW_profile( *this );
	}

#endif

#endif // end public variables and functions
};

} // end namespace brgastro

#endif /* _BRG_TNFW_PROFILE_H_ */
