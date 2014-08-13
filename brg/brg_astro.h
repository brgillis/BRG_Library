/**********************************************************************\
brg_astro.h
 -----------

 If this header is used, the source file brg_astro.cpp must be included
 and compiled with the project.

 This file includes various classes and functions for astrophysical
 objects and calculations. The file is split into five primary sections:

 -Constant Definitions
 -Class Forward Declarations
 -Static Class Definitions
 -Function Declarations
 -Class Definitions

 These sections are explained in further detail in their respective
 documentation blocks.

 With the exception of the Constant Definitions section, everything
 in this file is declared in the namespace brgastro.

 \**********************************************************************/

#ifndef __BRG_ASTRO_H_INCLUDED__
#define __BRG_ASTRO_H_INCLUDED__

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>

#include "brg_global.h"

#include "brg_units.h"

/** Constant Definitions **/
#if (1)
/**********************************************************************
 This section stores the definitions of various physical and
 astrophysical constants, stored as variables. Key points:

 -H_0 is set to 70 km/s/Mpc. If this is not changed, all results should
 be assumed to be in h_70 units.
 -Cosmological values are based on best-fit WMAP + priors values
 from Hinshaw et al. 2012.
 -If the _BRG_USE_UNITS_ parameter is set (see the brg_global.h file), values
 that have units will be declared as unit_objs.

 \**********************************************************************/

// Physical Constants
#if (1)
#ifndef __PI_DEFINED__
#define __PI_DEFINED__
const double pi = 3.14159265358979323846;
#endif


namespace brgastro
{

#ifndef __PHYS_VARS_DEFINED__
#define __PHYS_VARS_DEFINED__

#ifdef _BRG_USE_UNITS_
const brgastro::unit_obj Gc(6.67384e-11,3,-2,-1,0,0); // In m^3 s^-2 kg^-1
#else
const float Gc = 6.67384e-11; // In m^3 s^-2 kg^-1
#endif
const BRG_VELOCITY c = unitconv::ctomps;

#endif

#endif // end Physical Constants

// Cosmological Parameters
#if (1)
#ifndef __COSMO_VARS_DEFINED__
#define __COSMO_VARS_DEFINED__

#ifdef _BRG_USE_UNITS_
const brgastro::unit_obj H_0(70*unitconv::kmtom/unitconv::stos/unitconv::Mpctom,0,-1,0,0,0,0); // So all results will implicitly be in h_70 units
#else
const double H_0 = 70 * unitconv::kmtom / unitconv::stos / unitconv::Mpctom; // So all results will implicitly be in h_70 units
#endif

const float Omega_m = 0.288; // WMAP9 + priors
const float Omega_r = 0.000086; // WMAP9 + priors
const float Omega_k = 0; // Assuming flat space
const float Omega_l = 1 - Omega_k - Omega_m - Omega_r;
const float Omega_b = 0.0472; // WMAP9 + priors
const float sigma_8 = 0.830; // WMAP9 + priors
const float n_s = 0.971;     // WMAP9 + priors

#endif

#endif // end Cosmological Parameters

// Default parameters
#if (1)
#ifndef __BRG_DEFAULT_ASTRO_VARS_DEFINED__
#define __BRG_DEFAULT_ASTRO_VARS_DEFINED__

const double default_c = 6; // To help prevent crashes. Warning will be output
const double default_tau_factor = 2;

#endif
#endif

#endif // End Constant Definitions

/** Function Declarations **/
#if (1)
/**********************************************************************
 This section defines various functions which I find useful for
 astrophysical calculations. All are declared in the namespace brgastro.

 \**********************************************************************/

const BRG_UNITS H( const double z );

// Functions to get grid integers or grid boundaries from integers
const int get_ra_grid( const BRG_ANGLE &ra );
const int get_dec_grid( const BRG_ANGLE &dec );
const int get_z_grid( const double z );

const BRG_ANGLE get_ra_grid_lower( const int ra_grid );
const BRG_ANGLE get_dec_grid_lower( const int dec_grid );
const double get_z_grid_lower( const int z_grid );

const BRG_ANGLE get_ra_grid_upper( const int ra_grid );
const BRG_ANGLE get_dec_grid_upper( const int dec_grid );
const double get_z_grid_upper( const int z_grid );

const BRG_ANGLE get_ra_grid_mid( const int ra_grid );
const BRG_ANGLE get_dec_grid_mid( const int dec_grid );
const double get_z_grid_mid( const int z_grid );

// Functions to get transverse distance (in m) from angle (in rad) or vice-versa
const BRG_DISTANCE dfa( const BRG_ANGLE &da, const double z );
const BRG_DISTANCE dfa( const BRG_ANGLE &a1, const BRG_ANGLE &a2,
		const double z );
const BRG_DISTANCE dfa( const BRG_ANGLE &a1x, const BRG_ANGLE &a1y,
		const BRG_ANGLE &a2x, const BRG_ANGLE &a2y, const double z );

const BRG_ANGLE afd( const BRG_DISTANCE &dd, const double z );
const BRG_ANGLE afd( const BRG_DISTANCE &d1, const BRG_DISTANCE &d2,
		const double z );
const BRG_ANGLE afd( const BRG_DISTANCE &d1x, const BRG_DISTANCE &d1y,
		const BRG_DISTANCE &d2x, const BRG_DISTANCE &d2y, const double z );

// Functions to work between redshift, scale factor, and time (in s, with zero = present day)
const double zfa( const double a );
const double afz( const double z );

const BRG_TIME tfz( const double z );
const BRG_TIME tfa( const double z );
const double zft( const BRG_TIME &t );
const double aft( const BRG_TIME &t );

// Functions to integrate out distances
const double integrate_add( const double z1, const double z2 );
const double integrate_cmd( const double z1, const double z2 );
const double integrate_Ld( const double z1, const double z2 );
const double integrate_ltd( const double z1, const double z2 );
const double integrate_add( const double z );
const double integrate_cmd( const double z );
const double integrate_Ld( const double z );
const double integrate_ltd( const double z );
const double integrate_distance( const double z1, const double z2,
		const int mode, const int resolution = 10000 );

// Lensing functions
const BRG_DISTANCE ad_distance( double z1, double z2 = 0 );
const BRG_UNITS sigma_crit( const double z_lens, const double z_source );

#endif // end function declarations

/** Class Definitions **/
#if (1)
class redshift_obj
{
	/**********************************
	 redshift_obj class
	 ------------------

	 A base class for anything with a redshift.

	 Derived classes:

	 sky_obj
	 galaxy
	 group
	 density_profile
	 tNFW_profile
	 point_mass_profile

	 **********************************/
private:
	double _z_, _z_err_;
	mutable BRG_UNITS _H_cache_;
	mutable bool _H_cached_;
	mutable int _z_grid_;
	mutable bool _z_grid_cached_;
	mutable int _local_z_grid_change_num_;
public:

	// Constructor
	redshift_obj( const double init_z = 0, const double init_z_err = 0 )
	{
		_z_ = init_z;
		_z_err_ = init_z_err;
		_H_cache_ = 0;
		_H_cached_ = false;
		_z_grid_ = 0;
		_z_grid_cached_ = false;
		_local_z_grid_change_num_ = -1;
	}

	// Copy constructor
	// (implicit is fine for us)

	// Virtual destructor
	virtual ~redshift_obj()
	{
	}

	// Set functions
#if (1)
	virtual const int set_z( const double new_z ) // Sets z
	{
		_z_ = new_z;
		_H_cached_ = false;
		_z_grid_cached_ = false;
		return 0;
	}
	virtual const int set_z_err( const double new_z_err ) // Sets z error
	{
		_z_err_ = new_z_err;
		return 0;
	}
#endif

	// Get functions
#if (1)
	const double z() const
	{
		return _z_;
	}
	const double z_err() const
	{
		return _z_err_;
	}
	const int z_grid() const;
#endif

	// H(z) at the redshift of the object, given in units of m/s^2
	const BRG_UNITS H() const;

	// Clone function - Not needed in current implementation
	// virtual redshift_obj *redshift_obj_clone()=0;
};

#endif // end Class Definitions

} // end namespace brgastro

#endif
