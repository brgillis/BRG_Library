/**********************************************************************\
  @file astro.h

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

// body file: brg/physics/astro.cpp

#ifndef _BRG_ASTRO_H_INCLUDED_
#define _BRG_ASTRO_H_INCLUDED_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>

#include "IceBRG_main/common.h"

#include "IceBRG_main/units/unit_conversions.hpp"
#include "IceBRG_main/units/units.hpp"
#include "IceBRG_main/math/misc_math.hpp"


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
#ifndef _BRG_PI_DEFINED_
#define _BRG_PI_DEFINED_
const flt_t & pi = M_PI;
#endif


namespace IceBRG
{

#ifndef _PHYS_VARS_DEFINED_
#define _PHYS_VARS_DEFINED_

#ifdef _BRG_USE_UNITS_
const custom_unit_type<3,-2,-1,0,0> Gc(6.67384e-11 * cube(m) / square(s) / kg); // In m^3 s^-2 kg^-1
#else
constexpr flt_t Gc = 6.67384e-11; // In m^3 s^-2 kg^-1
#endif
const velocity_type c(IceBRG::unitconv::ctomps * mps);

#endif

#endif // end Physical Constants

// Cosmological Parameters
#if (1)
#ifndef _COSMO_VARS_DEFINED_
#define _COSMO_VARS_DEFINED_

const inverse_time_type H_0(
		70*unitconv::kmtom/unitconv::stos/unitconv::Mpctom / IceBRG::second); // So all results will implicitly be in h_70 units

constexpr flt_t Omega_m = 0.288; // WMAP9 + priors
constexpr flt_t Omega_r = 0.000086; // WMAP9 + priors
constexpr flt_t Omega_k = 0; // Assuming flat space
constexpr flt_t Omega_l = 1 - Omega_k - Omega_m - Omega_r;
constexpr flt_t Omega_b = 0.0472; // WMAP9 + priors
constexpr flt_t sigma_8 = 0.830; // WMAP9 + priors
constexpr flt_t n_s = 0.971;     // WMAP9 + priors

#endif

#endif // end Cosmological Parameters

// Default parameters
#if (1)
#ifndef _BRG_DEFAULT_ASTRO_VARS_DEFINED_
#define _BRG_DEFAULT_ASTRO_VARS_DEFINED_

constexpr flt_t default_c = 6.; // To help prevent crashes. Warning will be output
constexpr flt_t default_tau_factor = 2.;

#endif
#endif

#endif // End Constant Definitions

/** Function Declarations **/
#if (1)
/**********************************************************************
 This section defines various functions which I find useful for
 astrophysical calculations. All are declared in the namespace IceBRG.

 \**********************************************************************/

inverse_time_type H( const flt_t & z );

// Functions to get transverse distance_type (in m) from angle_type (in rad) or vice-versa
custom_unit_type<1,0,0,-1,0> dfa( const flt_t & z ); // Just gives conversion factor
distance_type dfa( const angle_type & da, const flt_t & z ); // Performs conversion
distance_type dfa( const angle_type & a1, const angle_type & a2,
		const flt_t & z ); // Performs conversion of dist between two angles
distance_type dfa( const angle_type & a1x, const angle_type & a1y,
		const angle_type & a2x, const angle_type & a2y, const flt_t & z ); // Performs conversion of dist between two positions

custom_unit_type<-1,0,0,1,0> afd( const flt_t & z );
angle_type afd( const distance_type & dd, const flt_t & z );
angle_type afd( const distance_type & d1, const distance_type & d2,
		const flt_t & z );
angle_type afd( const distance_type & d1x, const distance_type & d1y,
		const distance_type & d2x, const distance_type & d2y, const flt_t & z );

// Functions to work between redshift, scale factor, and time_type (in s, with zero = present day)
flt_t zfa( const flt_t & a );
flt_t afz( const flt_t & z );

time_type tfz( const flt_t & z );
time_type tfa( const flt_t & z );
flt_t zft( const time_type & t );
flt_t aft( const time_type & t );

// Functions to integrate out distances
distance_type integrate_add( const flt_t & z1, const flt_t & z2 );
distance_type integrate_cmd( const flt_t & z1, const flt_t & z2 );
distance_type integrate_Ld( const flt_t & z1, const flt_t & z2 );
distance_type integrate_ltd( const flt_t & z1, const flt_t & z2 );
distance_type integrate_add( const flt_t & z );
distance_type integrate_cmd( const flt_t & z );
distance_type integrate_Ld( const flt_t & z );
distance_type integrate_ltd( const flt_t & z );
distance_type integrate_distance( const flt_t & z1, const flt_t & z2,
		const int_t & mode, const int_t & resolution = 10000 );

// Lensing functions
distance_type ad_distance( flt_t z1, flt_t z2 = 0 );
surface_density_type sigma_crit( const flt_t & z_lens, const flt_t & z_source );

// Like dist2d, but using corrections for spherical geometry
template< typename Tr1, typename Td1, typename Tr2, typename Td2 >
inline const Tr1 skydist2d( const Tr1 & ra1, const Td1 & dec1,
		const Tr2 & ra2, const Td2 & dec2 )
{
	return quad_add( ( ra2 - ra1 ) * cos( ( dec2 + dec1 ) / 2. ), dec2 - dec1 );
}

// Luminosity functions

constexpr flt_t bright_abs_mag_max = -19.5;
constexpr flt_t faint_app_mag_max = 27;

distance_type get_lum_distance( flt_t const & z );
flt_t get_abs_mag_from_app_mag( flt_t const & app_mag, flt_t const & z );
custom_unit_type<-3,0,0,0,0> differential_luminosity_function( flt_t const & mag );
custom_unit_type<-3,0,0,0,0> integrated_luminosity_function( flt_t const & mag_lo, flt_t const & mag_hi );
flt_t faint_bright_ratio( flt_t const & z, flt_t const & bright_abs_mag_lim = bright_abs_mag_max,
		flt_t const & faint_app_mag_lim = faint_app_mag_max);

// Mass functions
flt_t delta_c();
density_type rho_bar( flt_t const & z = 0.);
distance_type r_of_m( mass_type const & mass, flt_t const & z );
flt_t sigma_of_r( distance_type const & r);
flt_t sigma_of_m( mass_type const & mass, flt_t const & z = 0. );
flt_t nu_of_m( mass_type const & mass, flt_t const & z = 0. );
custom_unit_type<-3,0,-1,0,0> mass_function( mass_type const & mass, flt_t const & z = 0. );
custom_unit_type<-3,0,0,0,0> log10_mass_function( flt_t const & log10msun_mass, flt_t const & z = 0. );

flt_t cluster_richness( mass_type const & mass, flt_t const & z,
		flt_t const & bright_abs_mag_lim = bright_abs_mag_max,
		flt_t const & faint_app_mag_lim = faint_app_mag_max );
mass_type min_cluster_mass( flt_t const & z,
		flt_t const & bright_abs_mag_lim = bright_abs_mag_max,
		flt_t const & faint_app_mag_lim = faint_app_mag_max );

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
	flt_t _z_, _z_err_;
	mutable inverse_time_type _H_cache_;
	mutable bool _H_cached_;
public:

	// Constructor
	redshift_obj( const flt_t & init_z = 0, const flt_t & init_z_err = 0 )
	: _z_(init_z),
	  _z_err_(init_z_err),
	  _H_cache_(0),
	  _H_cached_(false)
	{
	}

	// Copy constructor
	// (implicit is fine for us)

	// Virtual destructor
	virtual ~redshift_obj()
	{
	}

	// Set functions
#if (1)
	virtual void set_z( const flt_t & new_z ) // Sets z
	{
		_z_ = new_z;
		_H_cached_ = false;
	}
	virtual void set_z_err( const flt_t & new_z_err ) // Sets z error
	{
		_z_err_ = new_z_err;
	}
#endif

	// Get functions
#if (1)
	flt_t z() const
	{
		return _z_;
	}
	flt_t z_err() const
	{
		return _z_err_;
	}
	int_t z_grid() const;
#endif

	// Astro calculations

#if (1)
	// H at the redshift of the object, given in units of m/s^2
	inverse_time_type H() const;

	// Critical density_type at this redshift
	density_type rho_crit() const;

#endif

	// Clone function
	virtual redshift_obj *redshift_obj_clone() const
	{
		return new redshift_obj(*this);
	}
};

#endif // end Class Definitions

} // end namespace IceBRG

#endif
