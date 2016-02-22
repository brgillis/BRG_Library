/**********************************************************************\
  @file astro.cpp

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

#include <cstdlib>
#include <iostream>
#include <string>

#include "IceBRG_main/common.h"

#include "IceBRG_main/math/cache/cache.hpp"
#include "IceBRG_main/math/cache/cache_2d.hpp"
#include "IceBRG_main/math/calculus/differentiate.hpp"

#include "IceBRG_physics/astro_caches.h"
#include "IceBRG_physics/astro.h"

#include "IceBRG_main/units/units.hpp"

#include "IceBRG_physics/sky_obj/position_grid_cache.h"

namespace IceBRG {

// redshift_obj class methods
#if (1)

inverse_time_type redshift_obj::H() const
{
	// If not cached, calculate and cache it
	if(!_H_cached_)
	{
		if(_z_==0)
		{
			_H_cache_ = H_0;
		}
		else
		{
			_H_cache_ = IceBRG::H(_z_);
		}
		_H_cached_ = true;
	}
	return _H_cache_;
}

density_type redshift_obj::rho_crit() const
{
	return 3.*square(H())/(8.*pi*Gc);
}

#endif // redshift_obj class methods

/** Global Function Definitions **/
#if (1)

inverse_time_type H( const flt_t & test_z )
{
	// Friedmann equation, assuming omega = -1
	if(test_z==0) return H_0;
	flt_t zp1 = 1.+test_z;
	return H_0 * sqrt( Omega_r * quart( zp1 )
							+ Omega_m * cube( zp1 )
							+ Omega_k * square( zp1 ) + Omega_l );
}

// grid functions
#if (1)

#endif // end grid functions

// dfa and afd functions
#if (1)

custom_unit_type<1,0,0,-1,0> dfa( const flt_t & z )
{
	return dfa_cache().get( z );
}
distance_type dfa( const angle_type & da, const flt_t & z )
{
	return da * dfa(z);
}
distance_type dfa( const angle_type & a1, const angle_type & a2,
		const flt_t & z )
{
	return dfa( a2 - a1, z );
}
distance_type dfa( const angle_type & a1x, const angle_type & a1y,
		const angle_type & a2x, const angle_type & a2y, const flt_t & z )
{
	return dfa( skydist2d( a1x, a1y, a2x, a2y ), z );
}
custom_unit_type<-1,0,0,1,0> afd( const flt_t & z )
{
	return 1./dfa(z);
}
angle_type afd( const distance_type & dd, const flt_t & z )
{
	return dd / dfa(z);
}
angle_type afd( const distance_type & d1, const distance_type & d2,
		const flt_t & z )
{
	return afd( abs( d2 - d1 ), z );
}
angle_type afd( const distance_type & d1x,
		const distance_type & d1y, const distance_type & d2x,
		const distance_type & d2y, const flt_t & z )
{
	return afd( dist2d( d1x, d1y, d2x, d2y ), z );
}

flt_t zfa( const flt_t & a )
{
	return 1. / safe_d( a ) - 1.;
}
flt_t afz( const flt_t & z )
{
	return 1. / safe_d( 1 + z );
}

time_type tfz( const flt_t & z )
{
	return tfa( afz( z ) );
}
time_type tfa( const flt_t & a )
{
	return tfa_cache().get( a );
}
flt_t zft( const time_type & t )
{
	return zfa( aft( t ) );
}
flt_t aft( const time_type & t )
{
	tfa_cache cache;
	return cache.inverse_get( t );
}

#endif

// Functions to integrate out distances
#if (1)
distance_type integrate_add( const flt_t & z1, const flt_t & z2 )
{
	return integrate_distance( z1, z2, 0, 100000 );
}
distance_type integrate_cmd( const flt_t & z1, const flt_t & z2 )
{
	return integrate_distance( z1, z2, 1, 10000000 );
}
distance_type integrate_Ld( const flt_t & z1, const flt_t & z2 )
{
	return integrate_distance( z1, z2, 2, 10000000 );
}
distance_type integrate_ltd( const flt_t & z1, const flt_t & z2 )
{
	return integrate_distance( z1, z2, 3, 10000000 );
}
distance_type integrate_add( const flt_t & z )
{
	return integrate_distance( 0, z, 0, 100000 );
}
distance_type integrate_cmd( const flt_t & z )
{
	return integrate_distance( 0, z, 1, 10000000 );
}
distance_type integrate_Ld( const flt_t & z )
{
	return integrate_distance( 0, z, 2, 10000000 );
}
distance_type integrate_ltd( const flt_t & z )
{
	return integrate_distance( 0, z, 3, 10000000 );
}
distance_type integrate_distance( const flt_t & z1_init,
		const flt_t & z2_init, const int_t & mode, const int_t & n )
{
	// Function that will help calculate cosmic distances, thanks to Richard Powell - http://www.atlasoftheuniverse.com/
	// NOTE: Will return a negative value if z2 < z1. This is by design.

	flt_t OK0;
	flt_t OM0 = Omega_m, OR0 = Omega_r, OL0 = Omega_l;
	flt_t OM, OR, OL, OK;
	flt_t z1 = z1_init, z2 = z2_init;
	flt_t HD; //Hubble distance_type in billions of lightyears
	flt_t z, a, a1, a2, adot, h1;
	flt_t DC = std::numeric_limits<flt_t>::max(), DCC = 0, DT = std::numeric_limits<flt_t>::max(), DTT = 0, DM;
	int_t i;
	short_int_t sign = 1;

	if(z1==z2) return 0;

	OK0 = 1 - OM0 - OL0 - OR0;

	if ( z2 < z1 )
	{
		std::swap(z1,z2);
		sign = -1;
	}

	if ( z1 == 0 )
	{
		z = z2;
		h1 = value_of(H_0);
		OM = OM0;
		OR = OR0;
		OL = OL0;
		OK = OK0;
	}
	else
	{
		a1 = 1 / ( 1 + z1 );
		a2 = 1 / ( 1 + z2 );
		z = ( a1 / a2 ) - 1;
		h1 = value_of(H_0)
				* std::sqrt( OR0 * inv_quart( a1 ) + OM0 * inv_cube( a1 )
								+ OK0 * inv_square( a1 ) + OL0 );
		OM = OM0 * square( value_of(H_0) / h1 ) * inv_cube( a1 );
		OR = OR0 * square( value_of(H_0) / h1 ) * inv_quart( a1 );
		OL = OL0 * square( value_of(H_0) / h1 );
		OK = 1 - OM - OR - OL;
	}

	HD = value_of(c) / h1;

	for ( i = n; i >= 1; i-- )        // This loop is the numerical integration
	{
		a = ( i - 0.5 ) / n;              // Steadily decrease the scale factor
		// Comoving formula (See section 4 of Hogg, but I've added a radiation term too):
		adot = a * std::sqrt( OM * inv_cube(a) + OK * inv_square(a)
								+ OL + OR * inv_quart(a) ); // Note that "a" is equivalent to 1/(1+z)
		 // Collect DC and DT until the correct scale is reached
		DCC = DCC + 1 / ( a * adot ) / n; // Running total of the comoving distance_type
		DTT = DTT + 1 / adot / n; // Running total of the light travel time_type (see section 10 of Hogg)
		 // Collect DC and DT until the correct scale is reached
		DC = DCC;                 // Comoving distance_type DC
		DT = DTT;                 // Light travel time_type DT
		if ( a <= 1 / ( 1 + z ) ) // Collect DC and DT until the correct scale is reached
		{
			break;
		}
	}

	// Transverse comoving distance_type DM from section 5 of Hogg:
	if ( OK > 0.0001 )
		DM = ( 1 / std::sqrt( OK ) ) * sinh( std::sqrt( OK ) * DC );
	else if ( OK < -0.0001 )
		DM = ( 1 / std::sqrt( fabs( OK ) ) ) * std::sin( std::sqrt( fabs( OK ) ) * DC );
	else
		DM = DC;

	//age = HD*DTT;                 // Age of the universe (billions of years)
	//size = HD*DCC;                // Comoving radius of the observable universe

	flt_t res;

	switch ( mode )
	{
	case 0:
		res = sign * HD * DM / ( 1 + z ); // Angular diameter distance_type (section 6 of Hogg)
		break;
	case 1:
		res = sign * HD * DC;             // Comoving distance_type
		break;
	case 2:
		res = sign * HD * DM * ( 1 + z );  // Luminosity distance_type (section 7 of Hogg)
		break;
	case 3:
		res = sign * HD * DT;             // Light travel distance_type
		break;
	default:
		res = sign * HD * DT;             // Light travel distance_type
	}

	return units_cast<distance_type>(res);
}
#endif

// Lensing functions
#if (1)

distance_type ad_distance( flt_t z1, flt_t z2 )
{
	if ( z2 < z1 )
		std::swap( z1, z2 );
	return add_cache().get( z1, z2 );
}

surface_density_type sigma_crit( const flt_t & z_lens,
		const flt_t & z_source )
{
	return square( c ) / ( 4. * pi * Gc )
			* ad_distance( 0, z_source )
			/ ( ad_distance( 0, z_lens )
					* ad_distance( z_lens, z_source ) );

}
#endif // end lensing functions

// Luminosity-related functions
#if (1)

distance_type get_lum_distance( flt_t const & z )
{
	distance_type ad_distance = dfa( 1*rad, z );
	distance_type lum_distance = ad_distance*square(1+z);

	return lum_distance;
}

flt_t get_abs_mag_from_app_mag( flt_t const & app_mag, flt_t const & z )
{
	distance_type lum_distance_at_z = get_lum_distance(z);

	flt_t res = app_mag - 5*((std::log10(lum_distance_at_z/(1*unitconv::pctom*m)))-1);

	return res;
}

// Using results from Willmer et al. 2005 for the "All" sample with minimal weighting at
// z = 0.7 for the luminosity function

const custom_unit_type<-3,0,0,0,0> phi_star = 24.43 / cube( unitconv::Mpctom * m);
constexpr flt_t mag_star = -21.53;
constexpr flt_t alpha = -1.3;
constexpr flt_t abs_mag_min = -25;

custom_unit_type<-3,0,0,0,0> differential_luminosity_function( flt_t const & mag )
{
	custom_unit_type<-3,0,0,0,0> res = 0.4 * log(10.) * phi_star *
			pow(10.,0.4*(mag_star-mag)*(alpha+1))*exp(-pow(10,0.4*(mag_star-mag)));

	return res;
}
custom_unit_type<-3,0,0,0,0> integrated_luminosity_function( flt_t const & mag_lo, flt_t const & mag_hi )
{
	custom_unit_type<-3,0,0,0,0> res = lum_func_integral_cache().get(mag_lo,mag_hi);
	return res;
}

flt_t faint_bright_ratio( flt_t const & z, flt_t const & bright_abs_mag_lim,
		flt_t const & faint_app_mag_lim)
{
	flt_t faint_abs_mag_lim = get_abs_mag_from_app_mag( faint_app_mag_lim, z );

	flt_t res = integrated_luminosity_function( abs_mag_min, faint_abs_mag_lim ) /
			integrated_luminosity_function( abs_mag_min, bright_abs_mag_lim );

	return res;
}

#endif // end Luminosity-related functions

// Mass-related functions
#if(1)

// Press-Schechter formalism functions
flt_t delta_c()
{
	return 1.686;
}
density_type rho_bar( flt_t const & z)
{
	density_type rho_crit_0 = 3.*square(H_0)/(8.*pi*Gc);
	density_type rho_bar_0 = Omega_m * rho_crit_0;

	density_type res = rho_bar_0*cube(1.+z);

	return res;
}
distance_type r_of_m( mass_type const & mass, flt_t const & z )
{
	distance_type r = ipow<1,3>(3.*mass/(4*pi*rho_bar(z)));

	return r;
}
flt_t sigma_of_r( distance_type const & r)
{
	return 0.8*sigma_r_cache().get(r)/sigma_r_cache().get(8.*unitconv::Mpctom*m);
}
flt_t sigma_of_m( mass_type const & mass, flt_t const & z )
{
	return sigma_of_r(r_of_m(mass,z));
}
flt_t nu_of_m( mass_type const & mass, flt_t const & z )
{
	return delta_c()/sigma_of_m(mass,z);
}
flt_t fps_of_nu(flt_t const & nu)
{
	return sqrt(2./pi)*nu*exp(-square(nu)/2.);
}
flt_t fec_of_nu(flt_t const & nu)
{
	flt_t nu_tilde = 0.84*nu;
	return 0.322*(1+std::pow(nu_tilde,-0.6))*fps_of_nu(nu_tilde);
}
custom_unit_type<-3,0,-1,0,0> mass_function( mass_type const & mass, flt_t const & z )
{
	// Using Sheth, Mo, and Tormen (2001) formalism

	flt_t nu = nu_of_m(mass,z);

	flt_t f_of_nu = fec_of_nu(nu);

	auto ln_nu_of_m = [&z] (mass_type const & mass)
	{
		return log(nu_of_m(mass,z));
	};

	custom_unit_type<0,0,-1,0,0> d_ln_nu_d_m = differentiate(&ln_nu_of_m,mass);

	custom_unit_type<-3,0,-1,0,0> res = rho_bar(z)/mass * f_of_nu * d_ln_nu_d_m;

	return res;
}
inverse_volume_type log10_mass_function( flt_t const & log10msun_mass, flt_t const & z )
{
	mass_type mass = unitconv::Msuntokg*kg*std::pow(10.,log10msun_mass);

	inverse_volume_type res = mass_function(mass,z)*mass*std::log(10.);

	return res;
}

// Cluster richness
const mass_type richness_mstar = 4.07e12*unitconv::Msuntokg*kg;
constexpr flt_t richness_beta = 1.4;

flt_t cluster_richness( mass_type const & mass, flt_t const & z,
		flt_t const & bright_abs_mag_lim, flt_t const & faint_app_mag_lim )
{
	flt_t fb_ratio = faint_bright_ratio(z,bright_abs_mag_lim,faint_app_mag_lim);

	flt_t res = fb_ratio * log(mass/richness_mstar)/log(richness_beta);

	return res;
}
mass_type min_cluster_mass( flt_t const & z, flt_t const & bright_abs_mag_lim,
		flt_t const & faint_app_mag_lim )
{
	flt_t fb_ratio = faint_bright_ratio(z,bright_abs_mag_lim,faint_app_mag_lim);

	mass_type res = richness_mstar * pow(2/fb_ratio,richness_beta);

	return res;
}

custom_unit_type<0,0,0,-2,0> cluster_angular_density_at_z(flt_t const & z)
{
	inverse_volume_type cluster_volume_density = visible_cluster_density_cache().get(z);

	custom_unit_type<3,0,0,-2,0> vol_per_area = square(ad_distance(z)/rad)*c/H(z)/(1.+z);

	custom_unit_type<0,0,0,-2,0> density_per_area = cluster_volume_density*vol_per_area;

	return density_per_area;
};

flt_t visible_clusters( square_angle_type const & area, flt_t const & z1, flt_t const & z2 )
{
	custom_unit_type<0,0,0,-2,0> area_density = visible_clusters_cache().get(z1,z2);

	flt_t res = area*area_density;

	return res;
}

#endif // end Mass-related functions

#endif // end Global function definitions

} // namespace IceBRG
