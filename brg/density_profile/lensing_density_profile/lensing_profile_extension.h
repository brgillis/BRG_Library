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

#define IMPLEMENT_VIRTUAL_BRG_LENSING_EXTENSION_METHODS(class_name)         \
	virtual const BRG_DISTANCE rvir() const {return class_name::rvir();}    \
	virtual const BRG_UNITS dens( const BRG_DISTANCE &r ) const             \
		{return class_name::dens(r);}                                       \
	virtual const BRG_MASS enc_mass( const BRG_DISTANCE &r,                 \
		const bool silent =	true ) const;                                   \
			{return class_name::enc_mass(r,silent);}                        \
	virtual const int set_c( const double new_c, bool silent = false )      \
		{return class_name::set_c(new_c,silent);}

#define IMPLEMENT_BRG_LENSING_EXTENSION_METHODS(class_name)         \
	const BRG_DISTANCE rvir() const {return class_name::rvir();}    \
	const BRG_UNITS dens( const BRG_DISTANCE &r ) const             \
		{return class_name::dens(r);}                               \
	const BRG_MASS enc_mass( const BRG_DISTANCE &r,                 \
		const bool silent =	true ) const                            \
			{return class_name::enc_mass(r,silent);}                \
	const int set_c( const double new_c, bool silent = false )      \
		{return class_name::set_c(new_c,silent);}

namespace brgastro {

/**
 * Extensions to a density profile for weak gravitational lensing. This is a pure virtual class, and
 * only works once inherited by a subclass of density_profile.
 */
class lensing_profile_extension {
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
	virtual const BRG_DISTANCE rvir() const = 0;
	virtual const BRG_UNITS dens( const BRG_DISTANCE &r ) const = 0;
	virtual const BRG_MASS enc_mass( const BRG_DISTANCE &r, const bool silent =
			true ) const = 0;
	virtual const int set_c( const double new_c, bool silent = false ) = 0;

#endif // Virtual functions that will need to be implemented by base classes

	// Projected mass and density functions
#if (1)

	// These ones should be overridden if at all possible, as otherwise they'll have to be
	// integrated
#if (1)
	virtual const BRG_UNITS proj_dens( const BRG_DISTANCE &R,
			const bool silent = true ) const; // Projected surface density at radius R
	virtual const BRG_MASS proj_enc_mass( const BRG_DISTANCE &R,
			const bool silent = true ) const; // Mass enclosed within a cylinder of radius R
#endif

	// These ones typically don't need to be overridden, but they can be if it's convenient
#if (1)
	virtual const BRG_UNITS proj_enc_dens( const BRG_DISTANCE &R,
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
	const BRG_UNITS deltasigma( const BRG_DISTANCE &R, const bool silent =
			false ) const // Weak lensing signal in tangential shear Delta-Sigma at radius R
	{
		return proj_enc_dens( R, silent ) - proj_dens( R, silent );
	}
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
			const double group_c, const bool silent = true ) const; // As group_WLsig, but uses offset_WLsig cache to speed it up if overwritten
	virtual const BRG_UNITS quick_group_WLsig( const BRG_DISTANCE &R,
			const double group_c, const bool silent = true ) const // As deltasigma, but uses group_WLsig cache to speed it up if overwritten
	{
		return group_WLsig( R, group_c, silent );
	}
#endif

};

} // end namespace brgastro

#endif /* _BRG_LENSING_PROFILE_EXTENSION_H_ */
