/**       @file tNFW_profile.cpp
 *
 *     Project: brg
 *        Path: /brg/density_profile/tNFW_profile.cpp
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#include <iostream>
#include <cstdlib>
#include <vector>

#include "../brg_misc_functions.hpp"
#include "../brg_units.h"
#include "tNFW_profile.h"
#include "tNFW_caches.h"

// brgastro::tNFW_profile class methods
#if (1)

const double brgastro::tNFW_profile::_taufm( const double m_ratio,
		double precision, const bool silent ) const
{
	double m_target = m_ratio * _mftau();
	double taustepsize = _tau_ / 2;
	double tautest[3];
	double mtest, mbest;
	int i, ibest;

	tautest[0] = _tau_ / 2;
	mbest = 1e99;

	while ( taustepsize > precision * _c_ )
	{
		taustepsize /= 2;
		tautest[1] = tautest[0] - taustepsize;
		tautest[2] = tautest[0] + taustepsize;
		ibest = 0;
		for ( i = 0; i < 3; i++ )
		{
			mtest = _mftau( tautest[i] );
			if ( fabs( mtest - m_target ) <= fabs( mbest - m_target ) )
			{
				ibest = i;
				mbest = mtest;
			}
		}
		tautest[0] = tautest[ibest];
	}

	return tautest[0];

}

#if (1) // Constructors
brgastro::tNFW_profile::tNFW_profile()
{
	_mvir0_ = 0;
	_c_ = 0;
	_tau_ = 0;
}

brgastro::tNFW_profile::tNFW_profile( const BRG_MASS &init_mvir0,
		const double init_z, const double init_c, const double init_tau ) :
		brgastro::redshift_obj( init_z )
{
	_mvir0_ = init_mvir0;
	if ( init_c <= 0 )
	{
		_c_ = _cfm();
	}
	else
	{
		_c_ = init_c;
	}
	if ( init_tau < 0 )
	{
		_tau_ = default_tau_factor * _c_;
	}
	else
	{
		_tau_ = init_tau;
	}
}

#endif // End constructors

// Destructor
brgastro::tNFW_profile::~tNFW_profile()
{
}

#if (1) // Set functions

const int brgastro::tNFW_profile::set_mvir( const BRG_MASS &new_halo_mass,
		const bool silent )
{
	_mvir0_ = new_halo_mass;
	hmvir_cached = false;
	hmtot_cached = false;
	return 0;
}
const int brgastro::tNFW_profile::set_tau( const double new_halo_tau,
		const bool silent )
{
	_tau_ = new_halo_tau;
	hmvir_cached = false;
	hmtot_cached = false;
	return 0;
}
const int brgastro::tNFW_profile::set_c( const double new_halo_c,
		const bool silent )
{
	_c_ = new_halo_c;
	hmvir_cached = false;
	hmtot_cached = false;
	return 0;
}
const int brgastro::tNFW_profile::set_z( const double new_z )
{
	redshift_obj::set_z( new_z );
	hmvir_cached = false;
	hmtot_cached = false;
	return 0;
}
const int brgastro::tNFW_profile::set_parameters(
		const unsigned int num_parameters,
		const std::vector< BRG_UNITS > &parameters, const bool silent )
{
	if ( ( num_parameters != 4 ) || ( num_parameters != parameters.size() ) )
	{
		if ( !silent )
			if ( !silent )
				std::cerr
						<< "ERROR: Invalid number of parameters passed to tNFW_profile::set_parameters.\n"
						<< "Four are required for both num_parameters and parameters.size().\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	if ( set_mvir( parameters.at( 0 ) ) )
		return errorNOS( silent );
	if ( set_z( parameters.at( 1 ) ) )
		return errorNOS( silent );
	if ( parameters.at( 2 ) <= 0 )
	{
		if ( set_c( _cfm( parameters.at( 0 ), parameters.at( 1 ) ) ) )
			return errorNOS( silent );
	}
	else
	{
		if ( set_c( parameters.at( 2 ) ) )
			return errorNOS( silent );
	}
	if ( parameters.at( 3 ) <= 0 )
	{
		if ( set_tau( default_tau_factor * _c_ ) )
			return errorNOS( silent );
	}
	else
	{
		if ( set_tau( parameters.at( 3 ) ) )
			return errorNOS( silent );
	}
	return 0;
}

#endif // end set functions

const BRG_MASS brgastro::tNFW_profile::mvir() const
{
	return enc_mass(rvir()); // Not technically correct, but close enough for our purposes
}
const BRG_MASS brgastro::tNFW_profile::mvir0() const
{
	return _mvir0_;
}

const double brgastro::tNFW_profile::tau() const
{
	return _tau_;
}
const double brgastro::tNFW_profile::c() const
{
	return _c_;
}

const BRG_MASS brgastro::tNFW_profile::mtot() const
{
	return _mvir0_ * _mftau();
}

const BRG_VELOCITY brgastro::tNFW_profile::vvir() const
{
	return std::pow( 10 * Gc * H() * _mvir0_, 1. / 3. );
}
const BRG_DISTANCE brgastro::tNFW_profile::rvir() const
{
	return vvir() / H() / 10;
}
const BRG_DISTANCE brgastro::tNFW_profile::rs() const
{
	return rvir() / _c_;
}
const BRG_DISTANCE brgastro::tNFW_profile::rt( const bool silent ) const
{
	return rvir() / _tau_;
}

const BRG_MASS brgastro::tNFW_profile::hmvir() const
{
	return enc_mass( rvir() ) / 2;
}

const BRG_UNITS brgastro::tNFW_profile::quick_WLsig( const BRG_DISTANCE &r,
		const bool silent ) const
{
	BRG_UNITS result;

	result = brgastro::tNFW_sig_cache().get( z(), _mvir0_, r, silent );
	return result;
}
const BRG_UNITS brgastro::tNFW_profile::quick_offset_WLsig(
		const BRG_DISTANCE &r, const BRG_DISTANCE &offset_r,
		const bool silent ) const
{
	BRG_UNITS result;
	result = brgastro::tNFW_offset_sig_cache().get( z(), _mvir0_, r, offset_r,
			silent );
	return result;
}
const BRG_UNITS brgastro::tNFW_profile::semiquick_group_WLsig(
		const BRG_DISTANCE &r, const double group_c, const bool silent ) const
{
	BRG_UNITS result;
	if ( !silent )
		std::cerr << "ERROR: Placeholder function being used.\n";
	result = 0; // Placeholder!!!
	return result;
}
const BRG_UNITS brgastro::tNFW_profile::quick_group_WLsig(
		const BRG_DISTANCE &r, const double group_c, const bool silent ) const
{
	BRG_UNITS result;
	result = brgastro::tNFW_group_sig_cache().get( z(), _mvir0_, r, group_c,
			silent );
	return result;
}

const BRG_UNITS brgastro::tNFW_profile::dens( const BRG_DISTANCE &r ) const
{
	BRG_UNITS result, rho_c;

	double d_c, x, tau_use;
	if ( _tau_ <= 0 )
		tau_use = default_tau_factor * _c_;
	else
		tau_use = _tau_;
	d_c = _delta_c();
	rho_c = 3 * square(H()) / ( 8 * pi * Gc );
	x = r / rs();

	result = ( d_c * rho_c ) / ( x * square( 1 + x ) )
			* square( tau_use )
			/ ( square( tau_use ) + square( x ) );

	return result;
}
const BRG_UNITS brgastro::tNFW_profile::proj_dens( const BRG_DISTANCE &r,
		const bool silent ) const
{
	BRG_UNITS result, rho_c;
	double d_c, x, tau_use, fx, lx;

	if ( _tau_ <= 0 )
		tau_use = default_tau_factor * _c_;
	else
		tau_use = _tau_;

	d_c = _delta_c();
	rho_c = 3 * square(H()) / ( 8 * pi * Gc );
	x = r / rs();
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
				( 4 * pi * rs() * d_c * rho_c ) * tautau
						/ ( 2 * pi * square( tautaup1 ) )
						* ( tautaup1/3 + 2 * fx
								- pi / sqrt_tautaupxx
								+ ( tautau - 1 ) * lx
										/ ( tau_use
												* sqrt_tautaupxx ) );
	else
		result =
				( 4 * pi * rs() * d_c * rho_c ) * tautau
						/ ( 2 * pi * square( tautaup1 ) )
						* ( tautaup1 / ( xx - 1 )
								* ( 1 - fx ) + 2 * fx
								- pi / sqrt_tautaupxx
								+ ( tautau - 1 ) * lx
										/ ( tau_use
												* sqrt_tautaupxx ) );
	return result;
}
const BRG_MASS brgastro::tNFW_profile::enc_mass( const BRG_DISTANCE &r,
		const bool silent ) const
{
	using brgastro::square;
	using brgastro::cube;
	using std::log;
	using std::atan;

	BRG_UNITS rho_c;
	BRG_MASS m0, mx;
	double d_c, x, tau_use;
	if ( _tau_ < 0 )
		tau_use = default_tau_factor * _c_;
	else if (_tau_ == 0)
		return 0;
	else
		tau_use = _tau_;

	double tau_sq = square(tau_use);

	d_c = _delta_c();
	rho_c = 3 * square(H()) / ( 8 * pi * Gc );
	x = r / rs();

	// Result here integrated with Wolfram Alpha
	m0 = (2 * (1 + tau_sq) - (-1 + tau_sq) * 2 * log(tau_use));
	mx = ((2 * (1 + tau_sq)) / (1 + x)
							+ 4 * tau_use * atan(x / tau_use)
							+ 2 * (-1 + tau_sq) * log(1 + x)
							- (-1 + tau_sq) * log(tau_sq + square(x))) - m0;
	return cube(rs()) * d_c * rho_c * 2 * pi * tau_sq * mx
			/ square(1 + tau_sq);
}
const BRG_UNITS brgastro::tNFW_profile::proj_enc_dens( const BRG_DISTANCE &r,
		const bool silent ) const
{
	BRG_UNITS result, rho_c;
	//Takes M in kg, r in kpc
	double d_c, x, tau_use, fx, lx;
	if ( _tau_ <= 0 )
		tau_use = default_tau_factor * _c_;
	else
		tau_use = _tau_;

	d_c = _delta_c();
	rho_c = 3 * square(H()) / ( 8 * pi * Gc );
	x = r / rs();
	double xx = x*x;
	double tautau = tau_use*tau_use;
	double tautaup1 = tautau + 1.;
	double sqrt_tautaupxx = std::sqrt(tautau + xx);

	if ( x == 1. )
		fx = 1.;
	else if ( x > 1. )
		fx = acos( 1. / x ) / std::sqrt( xx - 1. );
	else
		fx = -log( 1 / x - std::sqrt( 1 / ( xx ) - 1 ) ) / std::sqrt( 1 - xx );
	lx = log( x / ( sqrt_tautaupxx + tau_use ) );

	result =
			( 4 * pi * rs() * d_c * rho_c ) * tautau
					/ ( pi * xx * square( tautaup1 ) )
					* ( ( tautaup1 + 2 * ( xx - 1 ) ) * fx
							+ pi * tau_use
							+ ( tautau - 1 ) * log( tau_use )
							+ sqrt_tautaupxx
									* ( lx * ( tautau - 1 )
											/ tau_use - pi ) );
	return result;
}
const BRG_MASS brgastro::tNFW_profile::proj_enc_mass( const BRG_DISTANCE &r,
		const bool silent ) const
{
	return proj_enc_dens( r, silent ) * pi * r * r;
}

const int brgastro::tNFW_profile::get_parameters( std::vector< BRG_UNITS > & parameters,
		const bool silent ) const
{
	parameters.resize( num_parameters() );

	try
	{
		parameters.at( 0 ) = _mvir0_;
		parameters.at( 1 ) = z();
		parameters.at( 2 ) = _c_;
		parameters.at( 3 ) = _tau_;
	}
	catch ( const std::exception & )
	{
		return errorNOS( silent );
	}
	return 0;
}

const int brgastro::tNFW_profile::get_parameter_names( std::vector< std::string > & parameter_names,
		const bool silent ) const
{
	parameter_names.resize( num_parameters() );

	try
	{
		parameter_names.at( 0 ) = "mvir0";
		parameter_names.at( 1 ) = "z";
		parameter_names.at( 2 ) = "c";
		parameter_names.at( 3 ) = "tau";
	}
	catch ( const std::exception & )
	{
		return errorNOS( silent );
	}
	return 0;
}

const int brgastro::tNFW_profile::truncate_to_fraction( const double f,
		const bool silent )
{
	if ( f <= 0 )
	{
		if ( f < 0 )
			if ( !silent )
				std::cerr
						<< "WARNING: Cannot truncate to negative fraction. Truncating to zero instead.\n";
		_tau_ = 0;
		override_rhmvir( 0 );
		override_rhmtot( 0 );
	}
	else
	{
		double new_tau_val = _taufm( f, bound(SMALL_FACTOR, std::fabs(0.1*(1-f)),0.00001) );
		if ( new_tau_val < 0 )
		{
			if ( !silent )
				std::cerr
						<< "WARNING: Could not solve for new tau value in truncate_to_fraction.\n"
						<< "tau will remain unchanged.\n";
		}
		else
		{
			_tau_ = new_tau_val;
		}
		hmvir_cached = false;
		hmtot_cached = false;
	}
	return 0;

}

#endif // end tNFW_profile functions
