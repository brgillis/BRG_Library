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

#include "brg_global.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include "cache/brg_cache.hpp"
#include "cache/brg_cache_2d.hpp"
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

/** Static Class Definitions **/
#if (1)
/**********************************************************************
 This section defines all "static" classes. That is, classes with
 entirely static variables. These classes are used as caches within the
 program to aid functions which use look-up tables to save time. The use
 of static member variables in this manner ensures that the look-up
 tables only need to be loaded once each.

 \**********************************************************************/

class grid_cache
{
private:
	static int _ra_grid_change_num_, _dec_grid_change_num_,
			_z_grid_change_num_;
	static BRG_ANGLE _ra_grid_min_, _ra_grid_max_, _ra_grid_step_;
	static BRG_ANGLE _dec_grid_min_, dec_grid_max_val, _dec_grid_step_;
	static double _z_grid_min_, _z_grid_max_, _z_grid_step_;
public:
	// Set functions
#if (1)
	const int set_ra_grid( const BRG_ANGLE new_ra_grid_min,
			const BRG_ANGLE new_ra_grid_max, const BRG_ANGLE new_ra_grid_step )
	{
		_ra_grid_min_ = new_ra_grid_min;
		_ra_grid_max_ = new_ra_grid_max;
		_ra_grid_step_ = new_ra_grid_step;
		_ra_grid_change_num_++;
		return 0;
	}

	const int set_dec_grid( const BRG_ANGLE new_dec_grid_min,
			const BRG_ANGLE new_dec_grid_max,
			const BRG_ANGLE new_dec_grid_step )
	{
		_dec_grid_min_ = new_dec_grid_min;
		dec_grid_max_val = new_dec_grid_max;
		_dec_grid_step_ = new_dec_grid_step;
		_dec_grid_change_num_++;
		return 0;
	}

	const int set_z_grid( const double new_z_grid_min,
			const double new_z_grid_max, const double new_z_grid_step )
	{
		_z_grid_min_ = new_z_grid_min;
		_z_grid_max_ = new_z_grid_max;
		_z_grid_step_ = new_z_grid_step;
		_z_grid_change_num_++;
		return 0;
	}
#endif

	// Get functions
#if (1)

	const int ra_grid_change_num()
	{
		return _ra_grid_change_num_;
	}
	const int dec_grid_change_num()
	{
		return _dec_grid_change_num_;
	}
	const int z_grid_change_num()
	{
		return _z_grid_change_num_;
	}
	const BRG_ANGLE ra_grid_min()
	{
		return _ra_grid_min_;
	}
	const BRG_ANGLE dec_grid_min()
	{
		return _dec_grid_min_;
	}
	const double z_grid_min()
	{
		return _z_grid_min_;
	}
	const BRG_ANGLE ra_grid_max()
	{
		return _ra_grid_max_;
	}
	const BRG_ANGLE dec_grid_max()
	{
		return dec_grid_max_val;
	}
	const double z_grid_max()
	{
		return _z_grid_max_;
	}
	const BRG_ANGLE ra_grid_step()
	{
		return _ra_grid_step_;
	}
	const BRG_ANGLE dec_grid_step()
	{
		return _dec_grid_step_;
	}
	const double z_grid_step()
	{
		return _z_grid_step_;
	}
#endif

};
// class grid_cache

class dfa_cache : public brg_cache<dfa_cache>
{
	// "Distance from angle" cache
private:

	DECLARE_BRG_CACHE_STATIC_VARS();

	friend class brg_cache<dfa_cache>;

protected:

	const std::string _name_base() const throw()
	{
		return "dfa";
	}

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0,0,1,0,0,0,0);
	}
	const brgastro::unit_obj _inverse_units() const throw()
	{
		return brgastro::unit_obj(0);
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function.
	const int _calculate( const double in_params, double & out_params ) const;

public:

};
// class dfa_cache

class add_cache : public brg_cache_2d<add_cache>
{
	// Angular diameter distance
private:

	DECLARE_BRG_CACHE_2D_STATIC_VARS();

	friend class brg_cache_2d<add_cache>;

protected:

	const std::string _name_base() const throw()
	{
		char name_base[BRG_CACHE_ND_NAME_SIZE] = "ang_di_d";
		return name_base;
	}

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0,1,0,0,0,0,0); // Distance units
	}
	const brgastro::unit_obj _inverse_units() const throw()
	{
		return brgastro::unit_obj(0); // Unitless (redshift)
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function.
	const int _calculate( const double in_param_1, const double in_param_2, double & out_param ) const;

public:

};
// class add_cache

class tfa_cache : public brg_cache<tfa_cache>
{
	// "Time from a (scale factor)" cache
private:

	DECLARE_BRG_CACHE_STATIC_VARS();

	friend class brg_cache<tfa_cache>;

protected:

	const std::string _name_base() const throw()
	{
		return "tfa";
	}

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0,0,1,0,0,0,0);
	}
	const brgastro::unit_obj _inverse_units() const throw()
	{
		return brgastro::unit_obj(0);
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function.
	const int _calculate( const double in_params, double & out_params ) const;

public:

};
// class tfa_cache



#endif // End static class definitions

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

