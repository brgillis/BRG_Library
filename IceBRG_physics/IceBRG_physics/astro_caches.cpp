/**********************************************************************\
  @file astro_caches.cpp

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
#include <exception>
#include <string>

#include "IceBRG_main/common.h"

#include "IceBRG_main/math/cache/cache.hpp"
#include "IceBRG_main/math/cache/cache_2d.hpp"
#include "IceBRG_main/math/calculus/integrate.hpp"
#include "IceBRG_main/units/units.hpp"

#include "IceBRG_physics/astro.h"
#include "astro_caches.h"

/** Static Class Initialisation **/
#if (1)

namespace IceBRG {

// Initialisation for IceBRG::dfa_cache
DEFINE_BRG_CACHE( dfa_cache, flt_t,
		decltype(custom_unit_type<1, 0, 0, -1, 0>()), 0, 5, 0.001
		,
			return integrate_add( 0, in_param )/radian;
		,

);

// Initialisation for IceBRG::add_cache
DEFINE_BRG_CACHE_2D( add_cache, flt_t, flt_t, distance_type,
		0, 5, 0.01, 0, 5, 0.01,
			return IceBRG::integrate_add(in_param_1,in_param_2);
		,

);

// Initialisation for IceBRG::tfa_cache
DEFINE_BRG_CACHE( tfa_cache, flt_t, time_type, 0.001, 1.02, 0.001
		,
			return -integrate_ltd( 0, zfa( in_param ) ) / c;
		,

);

// Initialisation for IceBRG::lum_func_integral_cache
DEFINE_BRG_CACHE_2D( lum_func_integral_cache, flt_t, flt_t, inverse_volume_type,
		-25, -11, 0.2,
		-25, -11, 0.2
		,
			inverse_volume_type res = integrate_Romberg(differential_luminosity_function,in_param_1,in_param_2);
			return res;
		,

);

// Initialisation for IceBRG::tfa_cache
DEFINE_BRG_CACHE( sigma_r_cache, distance_type, flt_t,
		0.1*unitconv::Mpctom*m, 100.*unitconv::Mpctom*m, 0.1*unitconv::Mpctom*m
		,
			constexpr flt_t nsp = 0.958;
			const distance_type Mpc = unitconv::Mpctom * m;
			const flt_t h = H_0 / (100. * (unitconv::kmtom*m)/s/Mpc);
			const distance_type Mpc_o_h = Mpc/h;
			const flt_t xi = Omega_m * h;

			distance_type aR = in_param/h;

			auto power_spectrum = [&] (inverse_distance_type const & k)
			{
				using std::log;
				using std::sin;
				using std::cos;

				flt_t res = pow(k*Mpc_o_h,2. + nsp) * square(log(1. + 2.34*Mpc_o_h * k / xi) /
						(2.34*Mpc_o_h * k / xi) /
						pow(1. + 3.89*Mpc_o_h * k / xi + square(16.2*Mpc_o_h * k / xi) + cube(5.47*Mpc_o_h * k / xi) +
								quart(6.71*Mpc_o_h * k / xi),0.25)) *
								square(3. * (sin(k * aR) - k * aR * cos(k * aR)) / cube(k * aR));

				return res;
			};

			flt_t res = Mpc_o_h*integrate_Romberg(power_spectrum,0.01/Mpc_o_h,100./Mpc_o_h)/(2.*square(pi));
			return res;
		,

);

// Initialisation for IceBRG::visible_clusters_cache
DEFINE_BRG_CACHE_2D( l10_mass_function_cache, flt_t, flt_t, inverse_volume_type,
		8, 16, 0.01,
		0.1, 1.3, 0.01
		,
			return log10_mass_function( in_param_1, in_param_2 );
		,
			sigma_r_cache().load();
);

// Initialisation for IceBRG::visible_cluster_density_cache
DEFINE_BRG_CACHE( visible_cluster_density_cache, flt_t, inverse_volume_type,
		0.1, 1.3, 0.01
		,
			mass_type min_mass = min_cluster_mass(in_param);
			flt_t l10_min_mass = std::log10(min_mass/(unitconv::Msuntokg*kg));

			l10_mass_function_cache l10_mf_cache;

			auto l10_mass_function_at_z = [&] (flt_t const & l10_m) {return l10_mf_cache.get(l10_m,in_param);};

			inverse_volume_type res = integrate_Romberg(l10_mass_function_at_z,l10_min_mass,16.);

			return res;
		,
			sigma_r_cache().load();
			l10_mass_function_cache().load();
);

// Initialisation for IceBRG::visible_clusters_cache
DEFINE_BRG_CACHE_2D( visible_clusters_cache, flt_t, flt_t, decltype(custom_unit_type<0,0,0,-2,0>()),
		0.1, 1.3, 0.01,
		0.1, 1.3, 0.01
		,

			decltype(custom_unit_type<0,0,0,-2,0>()) res = integrate_Romberg(cluster_angular_density_at_z,
					in_param_1,in_param_2);

			return res;
		,
			sigma_r_cache().load();
			l10_mass_function_cache().load();
			visible_cluster_density_cache().load();
);

} // namespace IceBRG

#endif // end Static Class Initialisation
