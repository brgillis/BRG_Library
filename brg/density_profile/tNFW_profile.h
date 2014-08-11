/**       @file tNFW_profile.h
 *
 *     Project: brg
 *        Path: /brg/density_profile/tNFW_profile.h
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#ifndef _BRG_TNFW_PROFILE_H_
#define _BRG_TNFW_PROFILE_H_

#include <cstdlib>
#include <vector>

#include "../brg_global.h"

#include "../brg_misc_functions.hpp"
#include "../brg_units.h"
#include "density_profile.h"

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
	double _c_, _tau_;

#endif // end private member variables

	const double _taufm( const double mratio, double precision = 0.00001,
			const bool silent = false ) const; // tau from Mtot/Mvir
	const double _delta_c() const // Simple function of concentration used as a step in calculating NFW densities
	{
		return ( 200. / 3. ) * cube( _c_ )
				/ ( log( 1 + _c_ ) - _c_ / ( 1 + _c_ ) );
	}

	// Functions relating to tNFW profiles
	const double _cfm( const BRG_MASS mass, const double z=0 ) const // Concentration from mass relationship, from Neto
	{
		return 4.67 * std::pow( mass * unitconv::kgtottMsun / ( 1e4 ), -0.11 );
	}
	const double _cfm() const // Concentration from mass relationship, from Neto
	{
		return _cfm(_mvir0_,z());
	}
	const double _mftau( const double tau, const double conc ) const // Mtot/Mvir from tau
	{
		if(tau<=0) return 0;
		double M0oM200 = 1 / ( log( 1 + conc ) - conc / ( 1 + conc ) );
		double tautau = tau*tau;
		double result =
				M0oM200 * tautau / square( tautau + 1 )
						* ( ( tautau - 1 ) * log( tau ) + tau * pi
								- ( tautau + 1 ) );
		return result;
	}
	const double _mftau( const double tau ) const // Mtot/Mvir from tau
	{
		return _mftau(tau,_c_);
	}
	const double _mftau() const // Mtot/Mvir from tau
	{
		return _mftau(_tau_,_c_);
	}
public:
#if (1) // public variables and functions

#if (1) // Constructors
	tNFW_profile();

	tNFW_profile( const BRG_MASS &init_mvir0, const double init_z,
			const double init_c = -1, const double init_tau = -1 );

#endif // End constructors

	// Destructor
	~tNFW_profile();

#if (1) // Set functions

	const int set_mvir( const BRG_MASS &new_halo_mass, const bool silent =
			false );
	const int set_parameters( const unsigned int num_parameters,
			const std::vector< BRG_UNITS > & new_parameters,
			const bool silent = false );

	const int set_z( const double new_z );
	const int set_tau( const double new_halo_tau, const bool silent = false );
	const int set_c( const double new_halo_c, const bool silent = false );

#endif // End set functions

#if (1) //Basic get functions

	const BRG_MASS mvir0() const;

	const BRG_MASS hmvir() const;
	const BRG_MASS mvir() const;
	const BRG_MASS mtot() const;

	const BRG_DISTANCE rvir() const;
	const BRG_DISTANCE rt( const bool silent = false ) const;
	const BRG_DISTANCE rs() const;

	const BRG_VELOCITY vvir() const;

	const double c() const;
	const double tau() const;
#endif // end basic get functions

#if (1) // advanced get functions

	const BRG_UNITS dens( const BRG_DISTANCE &r ) const;
	const BRG_UNITS proj_dens( const BRG_DISTANCE &R,
			const bool silent = false ) const;
	const BRG_MASS enc_mass( const BRG_DISTANCE &r,
			const bool silent = false ) const;
	const BRG_UNITS proj_enc_dens( const BRG_DISTANCE &R, const bool silent =
			false ) const;
	const BRG_MASS proj_enc_mass( const BRG_DISTANCE &R, const bool silent =
			false ) const;
	const BRG_UNITS quick_WLsig( const BRG_DISTANCE &R, const bool silent =
			false ) const;
	const BRG_UNITS quick_offset_WLsig( const BRG_DISTANCE &R,
			const BRG_DISTANCE &offset_R, const bool silent = false ) const;
	const BRG_UNITS semiquick_group_WLsig( const BRG_DISTANCE &R,
			const double group_c, const bool silent = false ) const;
	const BRG_UNITS quick_group_WLsig( const BRG_DISTANCE &R,
			const double group_c, const bool silent = false ) const;
	const unsigned int num_parameters() const
	{
		return 4; // Mass, redshift, c, and tau
	}
	const int get_parameters( std::vector< BRG_UNITS > & parameters,
			const bool silent = true ) const;

	const int get_parameter_names( std::vector< std::string > & parameter_names,
			const bool silent = true ) const;
#endif // end advanced get functions

#if (1) // Other operations

	virtual const int truncate_to_fraction( const double fraction,
			const bool silent = false );
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
