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

#ifndef _BRG_ASTRO_H_INCLUDED_
#define _BRG_ASTRO_H_INCLUDED_

#include "IceBRG_main/common.h"

#include "IceBRG_main/units/unit_conversions.hpp"
#include "IceBRG_main/units/units.hpp"
#include "IceBRG_main/math/misc_math.hpp"

#include "IceBRG_physics/constants.hpp"
#include "IceBRG_physics/cosmology.hpp"

namespace IceBRG
{

/** Function Declarations **/
#if (1)

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

distance_type ad_distance( flt_t z1, flt_t z2 = 0 );
distance_type comoving_distance( flt_t z );
distance_type luminosity_distance( flt_t z );

volume_type comoving_volume_element( flt_t z );

// Lensing functions
surface_density_type sigma_crit( const flt_t & z_lens, const flt_t & z_source );

// Like dist2d, but using corrections for spherical geometry
template< typename Tr1, typename Td1, typename Tr2, typename Td2 >
inline const Tr1 skydist2d( const Tr1 & ra1, const Td1 & dec1,
		const Tr2 & ra2, const Td2 & dec2 )
{
	return quad_add( ( ra2 - ra1 ) * cos( ( dec2 + dec1 ) / 2. ), dec2 - dec1 );
}

// Luminosity functions
#if(1)

flt_t get_abs_mag_from_app_mag( flt_t const & app_mag, flt_t const & z );
flt_t get_app_mag_from_abs_mag( flt_t const & abs_mag, flt_t const & z );
inverse_volume_type differential_luminosity_function( flt_t const & mag_B );
inverse_volume_type integrated_luminosity_function( flt_t const & mag_B_lo, flt_t const & mag_B_hi );
flt_t faint_bright_ratio( flt_t const & z, flt_t const & bright_abs_mag_i_lim = bright_abs_mag_i_max,
		flt_t const & faint_app_mag_i_lim = faint_app_mag_i_max);

#endif // end luminosity functions

// Mass functions
#if(1)

flt_t delta_c( flt_t const & z = 0. );
density_type rho_bar( flt_t const & z = 0. );
distance_type r_of_m( mass_type const & mass, flt_t const & z = 0. );
flt_t sigma_of_r( distance_type const & r);
flt_t sigma_of_m( mass_type const & mass );
flt_t nu_of_m( mass_type const & mass, flt_t const & z = 0. );
inverse_volume_inverse_mass_type mass_function( mass_type const & mass, flt_t const & z = 0. );
inverse_volume_type log10_mass_function( flt_t const & log10msun_mass, flt_t const & z = 0. );
inverse_volume_type integrated_log10_mass_function( flt_t const & l10_m_lo, flt_t const & l10_m_hi,
		flt_t const & z = 0. );

#endif // end mass functions

// Cluster visibility functions
#if(1)

flt_t cluster_richness( mass_type const & mass, flt_t const & z,
		flt_t const & bright_abs_mag_i_lim = bright_abs_mag_i_max,
		flt_t const & faint_app_mag_i_lim = faint_app_mag_i_max );
mass_type min_cluster_mass( flt_t const & z,
		flt_t const & bright_abs_mag_i_lim = bright_abs_mag_i_max,
		flt_t const & faint_app_mag_i_lim = faint_app_mag_i_max );

/**
 * Get the number density of clusters at a given redshift in units of
 * number per square radian per unit redshift.
 *
 * @param z
 * @return
 */
inverse_square_angle_type cluster_angular_density_at_z(flt_t const & z);
flt_t visible_clusters( square_angle_type const & area, flt_t const & z1 = 0.1, flt_t const & z2 = 1.3 );

flt_t integrate_mean_cluster_richness_at_redshift( flt_t const & z );
flt_t integrate_mean_cluster_richness( flt_t const & z_min, flt_t const & z_max );
flt_t mean_cluster_richness_at_redshift( flt_t const & z );
flt_t mean_cluster_richness( flt_t const & z_min, flt_t const & z_max );

#endif // Cluster visibility functions

// Abundance matching functions
#if(1)

flt_t get_abs_mag_B_from_mass( mass_type const & m, flt_t const & z );
mass_type get_mass_from_abs_mag_B( flt_t const & abs_mag, flt_t const & z );
flt_t get_app_mag_B_from_mass( mass_type const & m, flt_t const & z );
mass_type get_mass_from_app_mag_B( flt_t const & app_mag, flt_t const & z );

#endif // end abundance matching functions

// Galaxy visibility functions
#if(1)

flt_t max_galaxy_abs_mag_B( flt_t const & z,
		flt_t const & faint_app_mag_i_lim = faint_app_mag_i_max );

/**
 * Get the number density of galaxies at a given redshift in units of
 * number per square radian per unit redshift.
 *
 * @param z
 * @return
 */
inverse_square_angle_type galaxy_angular_density_at_z(flt_t const & z);
flt_t visible_galaxies( square_angle_type const & area, flt_t const & z1 = 0.1, flt_t const & z2 = 2.0 );

#endif // end galaxy visibility functions

// Stellar mass v abs mag functions
#if(1)

// Taken from regression of CFHTLenS data. Scatter of 1.2972528920157342
flt_t estimate_abs_mag_g_from_stellar_mass( mass_type const & stellar_mass );

// Taken from regression of CFHTLenS data. Log10 scatter of 0.585461144291201
mass_type estimate_stellar_mass_from_abs_mag_g( flt_t const & abs_mag_g );

// Taken from regression of CFHTLenS data. Scatter of 0.93030163251098796
flt_t estimate_abs_mag_i_from_stellar_mass( mass_type const & stellar_mass );

// Taken from regression of CFHTLenS data. Log10 scatter of 0.41169238482277282
mass_type estimate_stellar_mass_from_abs_mag_i( flt_t const & abs_mag_g );

flt_t estimate_abs_mag_g_from_abs_mag_i( flt_t const & abs_mag_i );
flt_t estimate_abs_mag_i_from_abs_mag_g( flt_t const & abs_mag_g );

#endif // stellar mass functions

#endif // end function declarations

} // end namespace IceBRG

#endif
