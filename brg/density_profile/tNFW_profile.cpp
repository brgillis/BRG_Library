/**       @file tNFW_profile.cpp
 *
 *     Project: brg
 *        Path: /brg/density_profile/tNFW_profile.cpp
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../brg_misc_functions.hpp"
#include "../brg_solvers.hpp"
#include "../brg_units.h"
#include "tNFW_profile_functors.hpp"

#include "tNFW_profile.h"

const double min_x = 0.0001;

// brgastro::tNFW_profile class methods
#if (1)

void brgastro::tNFW_profile::_uncache_mass()
{
	_rvir_cached_ = false;
	hmvir_cached = false;
	hmtot_cached = false;
}

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
	_rvir_cache_ = 0;
	_rvir_cached_ = false;
	_c_ = 0;
	_tau_ = 0;
}

brgastro::tNFW_profile::tNFW_profile( CONST_BRG_MASS_REF init_mvir0,
		const double init_z, const double init_c, const double init_tau )
{
	_mvir0_ = init_mvir0;
	_rvir_cache_ = 0;
	_rvir_cached_ = false;
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
	set_z(init_z);
}

#endif // End constructors

// Destructor
brgastro::tNFW_profile::~tNFW_profile()
{
}

#if (1) // Set functions

const int brgastro::tNFW_profile::set_mvir( CONST_BRG_MASS_REF new_halo_mass,
		const bool silent )
{
	_mvir0_ = new_halo_mass;
	_uncache_mass();
	return 0;
}
const int brgastro::tNFW_profile::set_tau( const double new_halo_tau,
		const bool silent )
{
	_tau_ = new_halo_tau;
	_uncache_mass();
	return 0;
}
const int brgastro::tNFW_profile::set_c( const double new_halo_c,
		const bool silent )
{
	_c_ = new_halo_c;
	_uncache_mass();
	return 0;
}
void brgastro::tNFW_profile::set_z( const double new_z )
{
	redshift_obj::set_z( new_z );
	_uncache_mass();
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
	set_mvir( parameters.at( 0 ) );
	set_z( parameters.at( 1 ) );
	if ( parameters.at( 2 ) <= 0 )
	{
		set_c( _cfm( parameters.at( 0 ), parameters.at( 1 ) ) );
	}
	else
	{
		set_c( parameters.at( 2 ) );
	}
	if ( parameters.at( 3 ) <= 0 )
	{
		set_tau( default_tau_factor * _c_ );
	}
	else
	{
		set_tau( parameters.at( 3 ) );
	}
	return 0;
}

#endif // end set functions

const BRG_MASS brgastro::tNFW_profile::mvir() const
{
	return enc_mass(rvir());
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
	return std::pow( 10 * Gc * H() * mvir(), 1. / 3. );
}
const BRG_VELOCITY brgastro::tNFW_profile::vvir0() const
{
	return std::pow( 10 * Gc * H() * _mvir0_, 1. / 3. );
}
const BRG_DISTANCE brgastro::tNFW_profile::rvir() const
{
	if(!_rvir_cached_)
	{

		tNFW_solve_rvir_iterative_functor it_solver(this);

		_rvir_cache_ = 0;

		// First, we try solving iteratively
		_rvir_cache_ = solve_iterate( &it_solver, rvir0(), 1, 0.0001, 1000, true );
		if ( ( _rvir_cache_ == 0 ) || ( isbad( _rvir_cache_ ) ) )
		{
			// Iteratively didn't work, so we go to the grid option
			tNFW_solve_rvir_minimize_functor min_solver(this);

			BRG_DISTANCE max_rvir = rvir0();
			try
			{
				_rvir_cache_ =  solve_grid( &min_solver, 0., max_rvir, 100, 0.);
			}
			catch(const std::exception &e)
			{
				std::cerr << "ERROR: Cannot solve virial radius for tNFW profile.\n";
				throw e;
			}
		}
		_rvir_cached_ = true;
	}
	return _rvir_cache_;
}
const BRG_DISTANCE brgastro::tNFW_profile::rvir0() const
{
	return vvir0() / H() / 10;
}
const BRG_DISTANCE brgastro::tNFW_profile::rs() const
{
	return rvir0() / _c_;
}
const BRG_DISTANCE brgastro::tNFW_profile::rt( const bool silent ) const
{
	return rvir0() / _tau_;
}

const BRG_UNITS brgastro::tNFW_profile::dens( CONST_BRG_DISTANCE_REF r ) const
{
	BRG_UNITS result, rho_c;

	double d_c, x, tau_use;
	if ( _tau_ <= 0 )
		tau_use = default_tau_factor * _c_;
	else
		tau_use = _tau_;
	d_c = _delta_c();
	rho_c = 3 * square(H()) / ( 8 * pi * Gc );
	x = max(r / rs(),min_x);

	result = ( d_c * rho_c ) / ( x * square( 1 + x ) )
			* square( tau_use )
			/ ( square( tau_use ) + square( x ) );

	return result;
}
const BRG_MASS brgastro::tNFW_profile::enc_mass( CONST_BRG_DISTANCE_REF r,
		const bool silent ) const
{
	using brgastro::square;
	using brgastro::cube;
	using std::log;
	using std::atan;

	BRG_UNITS rho_c;
	BRG_MASS m0, mx;
	double d_c, x, tau_use;
	if(r<=0) return 0;
	if ( _tau_ < 0 )
		tau_use = default_tau_factor * _c_;
	else if (_tau_ == 0)
		return 0;
	else
		tau_use = _tau_;

	double tau_sq = square(tau_use);

	d_c = _delta_c();
	rho_c = 3 * square(H()) / ( 8 * pi * Gc );
	x = max(r / rs(),min_x);

	// Result here integrated with Wolfram Alpha
	m0 = (2 * (1 + tau_sq) - (-1 + tau_sq) * 2 * log(tau_use));
	mx = ((2 * (1 + tau_sq)) / (1 + x)
							+ 4 * tau_use * atan(x / tau_use)
							+ 2 * (-1 + tau_sq) * log(1 + x)
							- (-1 + tau_sq) * log(tau_sq + square(x))) - m0;
	return cube(rs()) * d_c * rho_c * 2 * pi * tau_sq * mx
			/ square(1 + tau_sq);
}

const int brgastro::tNFW_profile::get_parameters( std::vector< BRG_UNITS > & parameters,
		const bool silent ) const
{
	parameters.resize( num_parameters() );
	parameters.at( 0 ) = _mvir0_;
	parameters.at( 1 ) = z();
	parameters.at( 2 ) = _c_;
	parameters.at( 3 ) = _tau_;
	return 0;
}

const int brgastro::tNFW_profile::get_parameter_names( std::vector< std::string > & parameter_names,
		const bool silent ) const
{
	parameter_names.resize( num_parameters() );

	parameter_names.at( 0 ) = "mvir0";
	parameter_names.at( 1 ) = "z";
	parameter_names.at( 2 ) = "c";
	parameter_names.at( 3 ) = "tau";
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
		_rvir_cache_ = 0;
		_rvir_cached_ = true;
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
		_uncache_mass();
	}
	return 0;

}

#endif // end tNFW_profile functions
