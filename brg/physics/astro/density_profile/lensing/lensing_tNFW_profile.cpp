/**********************************************************************\
  @file lensing_tNFW_profile.cpp

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

#include <vector>

#include "brg/global.h"

#include "brg/physics/astro/density_profile/lensing/lensing_tNFW_caches.h"

#include "lensing_tNFW_profile.h"

const double min_x = 0.000001;

const BRG_UNITS brgastro::lensing_tNFW_profile::quick_WLsig( CONST_BRG_DISTANCE_REF r,
		const bool silent ) const
{
	BRG_UNITS result = brgastro::tNFW_sig_cache().get( std::log(mvir0()), z(), std::log(r), silent );
	return result;
}
const BRG_UNITS brgastro::lensing_tNFW_profile::quick_offset_WLsig(
		CONST_BRG_DISTANCE_REF r, CONST_BRG_DISTANCE_REF offset_r,
		const bool silent ) const
{
	if(offset_r<=0) return WLsig(r,silent);
	BRG_UNITS result = brgastro::tNFW_offset_sig_cache().get( std::log(mvir0()), z(), std::log(r),
			std::log(offset_r), silent );
	return result;
}
const BRG_UNITS brgastro::lensing_tNFW_profile::quick_group_WLsig(
		CONST_BRG_DISTANCE_REF r, const double group_c, const bool silent ) const
{
	BRG_UNITS result = brgastro::tNFW_group_sig_cache().get( std::log(mvir0()), z(), std::log(r),
			group_c, silent );
	return result;
}
const BRG_UNITS brgastro::lensing_tNFW_profile::quick_shifted_WLsig(
		CONST_BRG_DISTANCE_REF R, const bool silent ) const
{
	BRG_UNITS result = brgastro::tNFW_shifted_sig_cache().get( std::log(mvir0()), z(), std::log(R),
			silent );
	return result;
}
const BRG_UNITS brgastro::lensing_tNFW_profile::proj_dens( CONST_BRG_DISTANCE_REF r,
		const bool silent ) const
{
	BRG_UNITS result, rho_c_t_4pi;
	double d_c, x, tau_use, fx, lx;

	if ( tau() <= 0 )
		tau_use = default_tau_factor * c();
	else
		tau_use = tau();

	d_c = _delta_c();
	rho_c_t_4pi = 3. * square(H()) / ( 2. * Gc );
	x = max(r / rs(),min_x);
	double xx = x*x;
	double tautau = tau_use*tau_use;
	double tautaup1 = tautau + 1.;
	double sqrt_tautaupxx = std::sqrt(tautau + xx);

	if ( x == 1. )
		fx = 1.;
	else if ( x > 1. )
		fx = acos( 1. / x ) / std::sqrt( xx - 1. );
	else
		fx = -log( 1. / x - std::sqrt( 1. / ( xx ) - 1. ) ) / std::sqrt( 1. - xx );
	lx = log( x / ( sqrt_tautaupxx + tau_use ) );
	if ( x == 1 )
		result =
				( rs() * d_c * rho_c_t_4pi ) * tautau
						/ ( 2. * pi * square( tautaup1 ) )
						* ( tautaup1/3. + 2. * fx
								- pi / sqrt_tautaupxx
								+ ( tautau - 1. ) * lx
										/ ( tau_use
												* sqrt_tautaupxx ) );
	else
		result =
				( rs() * d_c * rho_c_t_4pi ) * tautau
						/ ( 2. * pi * square( tautaup1 ) )
						* ( tautaup1 / ( xx - 1. )
								* ( 1. - fx ) + 2. * fx
								- pi / sqrt_tautaupxx
								+ ( tautau - 1. ) * lx
										/ ( tau_use
												* sqrt_tautaupxx ) );
	return result;
}
const BRG_UNITS brgastro::lensing_tNFW_profile::proj_enc_dens( CONST_BRG_DISTANCE_REF r,
		const bool silent ) const
{
	BRG_UNITS result, rho_c_t_4pi;
	//Takes M in kg, r in kpc
	long double d_c, x, tau_use, fx, lx;
	if ( tau() <= 0 )
		tau_use = default_tau_factor * c();
	else
		tau_use = tau();

	d_c = _delta_c();
	rho_c_t_4pi = 3. * square(H()) / ( 2. * Gc );
	x = max(r / rs(),min_x);
	long double xx = x*x;
	long double tautau = tau_use*tau_use;
	long double tautaup1 = tautau + 1.;
	long double sqrt_tautaupxx = std::sqrt(tautau + xx);

	if ( x == 1. )
		fx = 1.;
	else if ( x > 1. )
		fx = acos( 1. / x ) / std::sqrt( xx - 1. );
	else
		fx = -log( 1 / x - std::sqrt( 1 / ( xx ) - 1 ) ) / std::sqrt( 1 - xx );
	lx = log( x / ( sqrt_tautaupxx + tau_use ) );
	double log_tau = log(tau_use);

	result = ( rs() * d_c * rho_c_t_4pi ) * tautau
					/ ( pi * xx * tautaup1 * tautaup1 )
					* ( ( tautaup1 + 2 * ( xx - 1 ) ) * fx
							+ pi * tau_use
							+ ( tautau - 1 ) * log_tau
							+ sqrt_tautaupxx
									* ( lx * ( tautau - 1 )
											/ tau_use - pi ) );
	return result;
}
const BRG_MASS brgastro::lensing_tNFW_profile::proj_enc_mass( CONST_BRG_DISTANCE_REF r,
		const bool silent ) const
{
	return proj_enc_dens( r, silent ) * pi * r * r;
}
