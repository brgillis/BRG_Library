/**       @file density_profile.h
 *
 *     Project: brg
 *        Path: /brg/density_profile/density_profile.h
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#ifndef _BRG_DENSITY_PROFILE_H_
#define _BRG_DENSITY_PROFILE_H_

#include <vector>

#include "../brg_global.h"

#include "../brg_astro.h"
#include "../brg_units.h"

namespace brgastro {

class density_profile: public virtual redshift_obj
{
	/**********************************
	 density_profile class
	 -------------

	 An abstract class for anything which has
	 a mass density profile. Meant to be
	 overriden by child classes.

	 The methods for this class are sorted
	 into four groups:
	 -Pure virtual functions (must be implemented
	 by any derived classes)
	 -Virtual functions which must be overwritten if they're going
	 to be used (not pure virtual since they may not be needed)
	 -Virtual functions which should be overwritten if possible and
	 if they're going to be used (those that require integration,
	 but may have an analytic solution for certain profiles)
	 -Virtual functions that don't need to be overwritten in most
	 cases and non-virtual functions non-virtual functions

	 Derived classes:
	 tNFW_profile
	 point_mass profile

	 **********************************/

private:
#if (1)
	int _hm_type_;

	mutable BRG_DISTANCE _rhmvir_cache_, _rhmtot_cache_;

	BRG_VELOCITY _v_from_r( BRG_DISTANCE r ) const
	{
		BRG_UNITS a;

		a = accel( fabs( r ) );
		if ( a >= 0 )
			return 0;
		return std::sqrt( -a * fabs( r ) );
	}
#endif

protected:
	mutable bool hmvir_cached, hmtot_cached;

public:
#if (1)

#if (1) // Constructor
	density_profile()
	{
		_hm_type_ = 1;
		_rhmvir_cache_ = 0;
		hmvir_cached = false;
		_rhmtot_cache_ = 0;
		hmtot_cached = false;
		set_z( 0 );
	}
#endif

#if (1) // Destructor
	virtual ~density_profile()
	{
	}
#endif

#if (1) // Pure virtual functions (must be implemented by any derived classes)

	virtual const BRG_MASS mvir() const =0; // Virial mass (exact definition can be chosen per profile)
	virtual const BRG_UNITS dens( const BRG_DISTANCE &r ) const =0; // Local density at radius r

	virtual density_profile *density_profile_clone() const =0; // Creates a clone of this
#endif

#if (1) // Virtual functions which must be overwritten if they're going to be used

#if (1) // Set functions - will return 1 if profile doesn't support this method of setting
	// All take values in default unit set (m, s, kg, K, rad, C)
	virtual const int set_mvir( const BRG_MASS &new_mvir, bool silent = false )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_mvir(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_vvir( const BRG_VELOCITY &new_vvir, bool silent =
			false )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_vvir(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_rvir( const BRG_DISTANCE &new_rvir, bool silent =
			false )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_rvir(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_rs( const BRG_DISTANCE &new_rs, bool silent = false ) // Scale radius
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_rs(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_rt( const BRG_DISTANCE &new_rt, bool silent = false ) // Tidal/truncation radius
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_rt(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_parameters( const unsigned int num_parameters,
			const std::vector< BRG_UNITS > &parameters, bool silent = false )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_parameters(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_tau( const double new_tau, bool silent = false ) // Truncation parameter
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_tau(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_c( const double new_c, bool silent = false ) // Concentration
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_c(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	const int override_rhmvir( const BRG_DISTANCE &new_rhmvir, bool silent =
			false )
	{
		_rhmvir_cache_ = new_rhmvir;
		hmvir_cached = true;
		return 0;
	}
	const int override_rhmtot( const BRG_DISTANCE &new_rhmtot, bool silent =
			false )
	{
		_rhmtot_cache_ = new_rhmtot;
		hmtot_cached = true;
		return 0;
	}

#endif

#if (1) // Basic get functions
	// All return values in default unit set (m, s, kg, K, rad, C)

	virtual const BRG_MASS mtot() const // Total mass
	{
		return 0;
	}

	virtual const BRG_DISTANCE rt( const bool silent = false ) const // Tidal/truncation radius
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::get_parameters(...) must be overridden to be used.\n";
		return 0;
	}

#endif // end basic get functions

	virtual const unsigned int num_parameters() const
	{
		throw std::runtime_error("ERROR: num_parameters must be overridden if it's used.");
		return 0;
	}

	virtual const int get_parameters( std::vector< BRG_UNITS > & parameters,
			const bool silent = false ) const // Returns a set of parameters for this halo (all the variables needed to define it - must be defined for each child)
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::get_parameters(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}

	virtual const int get_parameter_names( std::vector< std::string > & parameter_names,
			const bool silent = false ) const // Returns a set of names of this halo's parameters
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::get_parameter_names(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}

#endif // end Virtual functions which must be overwritten if they're going to be used

#if (1) // Virtual functions which should be overwritten if at all possible and if they'll be used
	virtual const BRG_MASS enc_mass( const BRG_DISTANCE &r, const bool silent =
			true ) const; // Mass enclosed with sphere of radius r
	virtual const BRG_UNITS proj_dens( const BRG_DISTANCE &R,
			const bool silent = true ) const; // Projected surface density at radius R
	virtual const BRG_MASS proj_enc_mass( const BRG_DISTANCE &R,
			const bool silent = true ) const; // Mass enclosed within a cylinder of radius R
	virtual const BRG_UNITS offset_WLsig( const BRG_DISTANCE &R,
			const BRG_DISTANCE &offset_R, const bool silent = true ) const; // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from position offset by offset_R
	virtual const BRG_UNITS group_WLsig( const BRG_DISTANCE &R,
			const double group_c, const bool silent = true ) const; // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from ensemble of satellites in group with satellite concentration group_c
	virtual const BRG_UNITS quick_WLsig( const BRG_DISTANCE &R,
			const bool silent = true ) const // As deltasigma, but uses cache to speed it up if overwritten
	{
		return deltasigma( R, silent );
	}
	virtual const BRG_UNITS quick_offset_WLsig( const BRG_DISTANCE &R,
			const BRG_DISTANCE &offset_R, const bool silent = true ) const // As offset_WLsig, but uses cache to speed it up if overwritten
	{
		return offset_WLsig( R, offset_R, silent );
	}
	virtual const BRG_UNITS semiquick_group_WLsig( const BRG_DISTANCE &R,
			const double group_c, const bool silent = true ) const // As group_WLsig, but uses offset_WLsig cache to speed it up if overwritten
	{
		return group_WLsig( R, group_c, silent );
	}
	virtual const BRG_UNITS quick_group_WLsig( const BRG_DISTANCE &R,
			const double group_c, const bool silent = true ) const // As deltasigma, but uses group_WLsig cache to speed it up if overwritten
	{
		return group_WLsig( R, group_c, silent );
	}
#endif

#if (1) // Virtual functions which shouldn't be overwritten in most cases and non-virtual functions

	virtual const BRG_DISTANCE rvir() const // Virial radius (exact definition can be chosen per profile if preferred)
	{
		double virial_density_factor = 200;

		return safe_pow(
				2 * mvir() * Gc / ( square( H() ) * virial_density_factor ),
				1. / 3. );
	}
	const BRG_MASS hmvir() const // Half virial mass
	{
		return mvir() / 2;
	}
	const BRG_MASS hmtot() const // Half total mass
	{
		return mtot() / 2;
	}
	virtual const BRG_MASS hm() const // Half mass (which depends on hm_type value)
	{
		return ( _hm_type_ == 0 ? hmvir() : hmtot() );
	}
	virtual const BRG_UNITS enc_dens( const BRG_DISTANCE &r,
			const bool silent = false ) const // Mean density enclosed with sphere of radius r
	{
		BRG_DISTANCE r_to_use = max( std::fabs( r ), SMALL_FACTOR );
		return enc_mass( r_to_use, silent )
				/ ( 4. / 3. * pi * cube( std::fabs( r_to_use ) ) );
	}
	virtual const BRG_UNITS proj_enc_dens( const BRG_DISTANCE &R,
			const bool silent = false ) const // Mean surface density enclosed within a cylinder of radius R
	{
		BRG_DISTANCE R_to_use = max( std::fabs( R ), SMALL_FACTOR );
		return proj_enc_mass( R_to_use, silent )
				/ ( pi * square( std::fabs( R_to_use ) ) );
	}
	virtual const BRG_DISTANCE rhmvir( const bool silent = false ) const; // Half-virial-mass radius
	virtual const BRG_DISTANCE rhmtot( const bool silent = false ) const; // Half-total-mass radius
	virtual const BRG_DISTANCE rhm( const bool silent = false ) const // Half-mass radius
	{
		return ( _hm_type_ == 0 ? rhmvir( silent ) : rhmtot( silent ) );
	}
	virtual const BRG_VELOCITY vvir() const // Orbital velocity at rvir
	{
		return _v_from_r( rvir() );
	}
	virtual const BRG_VELOCITY vhmvir( const bool silent = false ) const // Orbital velocity at rhmvir
	{
		return _v_from_r( rhmvir( silent ) );
	}
	virtual const BRG_VELOCITY vhmtot( const bool silent = false ) const // Orbital velocity at rhmtot
	{
		return _v_from_r( rhmtot( silent ) );
	}
	virtual const BRG_VELOCITY vhm( const bool silent = false ) const // Orbital velocity at rhm
	{
		return ( _hm_type_ == 0 ? vhmvir( silent ) : vhmtot( silent ) );
	}
	virtual const BRG_VELOCITY vt( const bool silent = false ) const // Orbital velocity at rt
	{
		return _v_from_r( rt( silent ) );
	}
	const BRG_TIME otvir() const // Orbital period at rvir
	{
		return 2 * pi * rvir() / vvir();
	}
	const BRG_TIME othmvir( const bool silent = false ) const // Orbital period at rhmvir
	{
		return 2 * pi * rhmvir( silent ) / safe_d(vhmvir( silent ));
	}
	const BRG_TIME othmtot( const bool silent = false ) const // Orbital period at rhmtot
	{
		return 2 * pi * rhmtot( silent ) / safe_d(vhmtot( silent ));
	}
	const BRG_TIME othm( const bool silent = false ) const // Orbital period at rhm
	{
		return 2 * pi * rhm( silent ) / safe_d(vhm( silent ));
	}
	const BRG_TIME ott( const bool silent = false ) const // Orbital period at rt
	{
		return 2 * pi * rt( silent ) / safe_d(vt( silent ));
	}

	const int set_hm_type( int new_hm_type ) // Whether default half-mass is half virial, half total, or something else
	{
		_hm_type_ = new_hm_type;
		return 0;
	}

#if (1) // advanced get functions

	const BRG_UNITS accel( const BRG_DISTANCE &r,
			const bool silent = false ) const // Gravitational acceleration at radius r
	{
		return r == 0 ? 0 : -Gc * enc_mass( r, silent ) / square( r );
	}
	virtual const BRG_UNITS Daccel( const BRG_DISTANCE &r, const bool silent =
			false ) const; // Derivative of acceleration at radius r
	const BRG_UNITS deltasigma( const BRG_DISTANCE &R, const bool silent =
			false ) const // Weak lensing signal in tangential shear Delta-Sigma at radius R
	{
		return proj_enc_dens( R, silent ) - proj_dens( R, silent );
	}

#endif // end advanced get functions

#endif // end Virtual functions which shouldn't be overwritten in most cases and non-virtual functions

#if (1) // Other operations

	virtual const int truncate_to_fraction( const double fraction,
			bool silent = false ) // Adjusts parameters of this class to decrease the mass to fraction of its previous mass. Must be defined for each child
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::truncate_to_fraction(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}

#endif

#endif
}; // end class density_profile

// Function to estimate orbital period from current position and velocity in a density profile
// Note that this is merely an estimate from analogy to calculations in a Keplerian potential
const BRG_TIME period( const density_profile *host, const BRG_DISTANCE &r,
		const BRG_VELOCITY &vr, const BRG_VELOCITY &vt = 0 );

} // end namespace brgastro

#endif /* _BRG_DENSITY_PROFILE_H_ */
