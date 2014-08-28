/**       @file lensing_profile_extension.h
 *
 *     Project: brg
 *        Path: /brg/density_profile/lensing_profile_extension.h
 *
 *  Created on: 12 Aug 2014
 *      Author: brg
 */

#ifndef _BRG_LENSING_PROFILE_EXTENSION_H_
#define _BRG_LENSING_PROFILE_EXTENSION_H_

#include "../../brg_global.h"

#include "../../brg_units.h"
#include "../density_profile.h"

#define IMPLEMENT_VIRTUAL_BRG_LENSING_EXTENSION_METHODS(class_name)   \
	virtual double z() const {return class_name::z();}                \
	virtual BRG_DISTANCE rvir() const {return class_name::rvir();}    \
	virtual BRG_UNITS dens( CONST_BRG_DISTANCE_REF r ) const          \
		{return class_name::dens(r);}                                 \
	virtual BRG_MASS enc_mass( CONST_BRG_DISTANCE_REF r,              \
		const bool silent =	true ) const                              \
			{return class_name::enc_mass(r,silent);}                  \
	virtual double c() const {return class_name::c();}                \
	virtual void set_c( const double new_c, bool silent = false )     \
		{class_name::set_c(new_c,silent);}                            \
	virtual double tau() const {return class_name::tau();}            \
	virtual void set_tau( const double new_tau, bool silent = false ) \
		{class_name::set_tau(new_tau,silent);}

#define IMPLEMENT_BRG_LENSING_EXTENSION_METHODS(class_name)   \
	double z() const {return class_name::z();}                \
	BRG_DISTANCE rvir() const {return class_name::rvir();}    \
	BRG_UNITS dens( CONST_BRG_DISTANCE_REF r ) const          \
		{return class_name::dens(r);}                         \
	BRG_MASS enc_mass( CONST_BRG_DISTANCE_REF r,              \
		const bool silent =	true ) const                      \
			{return class_name::enc_mass(r,silent);}          \
	double c() const {return class_name::c();}                \
	void set_c( const double new_c, bool silent = false )     \
		{class_name::set_c(new_c,silent);}                    \
	double tau() const {return class_name::tau();}            \
	void set_tau( const double new_tau, bool silent = false ) \
		{class_name::set_tau(new_tau,silent);}

namespace brgastro {

/**
 * Extensions to a density profile for weak gravitational lensing. This is a pure virtual class, and
 * only works once inherited by a subclass of density_profile.
 */
class lensing_profile_extension {
private:

	const BRG_DISTANCE _shift_sigma( CONST_BRG_DISTANCE_REF R, const bool silent = true ) const;

public:

	// Constructors and destructors
#if (1)
	lensing_profile_extension()
	{
	}
	virtual ~lensing_profile_extension()
	{
	}
#endif

	// Virtual functions that will need to be implemented by base classes
#if (1)

	// Clone function for this
	virtual lensing_profile_extension * lensing_profile_extension_clone() const = 0;

	// Get functions for values that will be needed for calculations
	virtual double z() const = 0;
	virtual BRG_DISTANCE rvir() const = 0;
	virtual BRG_UNITS dens( CONST_BRG_DISTANCE_REF r ) const = 0;
	virtual BRG_MASS enc_mass( CONST_BRG_DISTANCE_REF r, const bool silent =
			true ) const = 0;
	virtual double c() const = 0;
	virtual void set_c( const double new_c, bool silent = false ) = 0;
	virtual double tau() const = 0;
	virtual void set_tau( const double new_c, bool silent = false ) = 0;

#endif // Virtual functions that will need to be implemented by base classes

	// Projected mass and density functions
#if (1)

	// These ones should be overridden if at all possible, as otherwise they'll have to be
	// integrated
#if (1)
	virtual const BRG_UNITS proj_dens( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const; // Projected surface density at radius R
	virtual const BRG_MASS proj_enc_mass( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const; // Mass enclosed within a cylinder of radius R
#endif

	// These ones typically don't need to be overridden, but they can be if it's convenient
#if (1)
	virtual const BRG_UNITS proj_enc_dens( CONST_BRG_DISTANCE_REF R,
			const bool silent = false ) const // Mean surface density enclosed within a cylinder of radius R
	{
		BRG_DISTANCE R_to_use = max( std::fabs( R ), SMALL_FACTOR );
		return proj_enc_mass( R_to_use, silent )
				/ ( pi * square( R_to_use ) );
	}
#endif

#endif // Projected mass and density functions

	// Weak lensing functions
#if (1)

	// WLsig (delta sigma)
#if (1)
	const BRG_UNITS WLsig( CONST_BRG_DISTANCE_REF R, const bool silent =
			false ) const // Weak lensing signal in tangential shear Delta-Sigma at radius R
	{
		return proj_enc_dens( R, silent ) - proj_dens( R, silent );
	}
	virtual const BRG_UNITS quick_WLsig( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const // As deltasigma, but uses cache to speed it up if overwritten
	{
		return WLsig( R, silent );
	}
#endif

	// Offset WLsig
#if (1)
	virtual const BRG_UNITS offset_WLsig( CONST_BRG_DISTANCE_REF R,
			CONST_BRG_DISTANCE_REF offset_R, const bool silent = true ) const; // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from position offset by offset_R
	virtual const BRG_UNITS quick_offset_WLsig( CONST_BRG_DISTANCE_REF R,
			CONST_BRG_DISTANCE_REF offset_R, const bool silent = true ) const // As offset_WLsig, but uses cache to speed it up if overwritten
	{
		return offset_WLsig( R, offset_R, silent );
	}
#endif

	// Group WLsig
#if (1)
	virtual const BRG_UNITS group_WLsig( CONST_BRG_DISTANCE_REF R,
			const double group_c, const bool silent = true ) const; // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from ensemble of satellites in group with satellite concentration group_c
	virtual const BRG_UNITS semiquick_group_WLsig( CONST_BRG_DISTANCE_REF R,
			const double group_c, const bool silent = true ) const; // As group_WLsig, but uses offset_WLsig cache to speed it up if overwritten
	virtual const BRG_UNITS quick_group_WLsig( CONST_BRG_DISTANCE_REF R,
			const double group_c, const bool silent = true ) const // As deltasigma, but uses group_WLsig cache to speed it up if overwritten
	{
		return semiquick_group_WLsig( R, group_c, silent );
	}
#endif

	// Shifted WLsig
#if (1)
	// This represents the weak-lensing signal after being corrected for errors due to relative
	// shifting of the lens and source due to an intervening mass distribution

	const double shift_factor( CONST_BRG_DISTANCE_REF R, const bool silent = true ) const;

	virtual const BRG_UNITS shifted_WLsig( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const;
	virtual const BRG_UNITS semiquick_shifted_WLsig( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const;
	virtual const BRG_UNITS quick_shifted_WLsig( CONST_BRG_DISTANCE_REF R,
			const bool silent = true ) const
	{
		return semiquick_shifted_WLsig( R, silent );
	}
#endif

#endif

};

} // end namespace brgastro

#endif /* _BRG_LENSING_PROFILE_EXTENSION_H_ */
