/**       @file point_mass_profile.cpp
 *
 *     Project: brg
 *        Path: /brg/density_profile/point_mass_profile.cpp
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#include <cstdlib>
#include <iostream>
#include <vector>

#include "../brg_global.h"

#include "../brg_units.h"
#include "point_mass_profile.h"

// brgastro::point_mass_profile class methods
#if (1)

#if (1) // Constructors
brgastro::point_mass_profile::point_mass_profile()
{
	_mass_ = 0;
}

brgastro::point_mass_profile::point_mass_profile( const BRG_MASS init_mass,
		const double init_z )
{
	set_mvir( init_mass );
	set_z( init_z );
}

#endif // End constructors

// Destructor
brgastro::point_mass_profile::~point_mass_profile()
{
}

#if (1) // Set functions

const int brgastro::point_mass_profile::set_mvir(
		const BRG_MASS &new_halo_mass, bool silent )
{
	_mass_ = new_halo_mass;
	return 0;
}
const int brgastro::point_mass_profile::set_parameters(
		const unsigned int num_parameters,
		const std::vector< BRG_UNITS > &parameters, bool silent )
{

	if ( ( num_parameters != 2 ) || ( num_parameters != parameters.size() ) )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Invalid number of parameters passed to tNFW_profile::set_parameters.\n"
					<< "Four are required for both num_parameters and parameters.size().\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	if ( set_mvir( parameters.at( 0 ) ) )
		return errorNOS();
	if ( set_z( parameters.at( 1 ) ) )
		return errorNOS();
	return 0;
}

#endif // end set functions

const BRG_MASS brgastro::point_mass_profile::mvir() const
{
	return _mass_;
}
const BRG_MASS brgastro::point_mass_profile::mass() const
{
	return _mass_;
}

const BRG_MASS brgastro::point_mass_profile::mtot() const
{
	return _mass_;
}

const BRG_VELOCITY brgastro::point_mass_profile::vvir() const
{
	return std::pow( 10 * Gc * H() * mvir(), 1. / 3. );
}
const BRG_DISTANCE brgastro::point_mass_profile::rvir() const
{
	return vvir() / H() / 10;
}
const BRG_DISTANCE brgastro::point_mass_profile::rs() const
{
	return 0;
}
const BRG_DISTANCE brgastro::point_mass_profile::rt(
		const bool silent) const
{
	return 0;
}

const BRG_MASS brgastro::point_mass_profile::hmvir() const
{
	return 0;
}

const BRG_UNITS brgastro::point_mass_profile::dens(
		const BRG_DISTANCE &r ) const
{
#ifdef _BRG_USE_UNITS_
	BRG_UNITS result(0,-3,0,1,0,0,0);
#else
	double result = 0;
#endif

	result = ( r == 0 ? DBL_MAX : 0 );

	return result;
}
const BRG_UNITS brgastro::point_mass_profile::enc_dens(
		const BRG_DISTANCE &r,
		const bool silent ) const
{
	return enc_mass( r ) / ( 4. / 3. * pi * cube( r ) );
}
const BRG_MASS brgastro::point_mass_profile::enc_mass(
		const BRG_DISTANCE &r,
		const bool silent ) const
{
	return _mass_;
}

const int brgastro::point_mass_profile::get_parameters( std::vector< BRG_UNITS > & parameters,
		const bool silent ) const
{
	parameters.resize( num_parameters() );

	try
	{
		parameters.at( 0 ) = _mass_;
		parameters.at( 1 ) = z();
	}
	catch ( std::exception & )
	{
		return errorNOS();
	}
	return 0;
}

const int brgastro::point_mass_profile::get_parameter_names(std::vector< std::string > & parameter_names,
		const bool silent ) const
{
	parameter_names.resize( num_parameters() );

	try
	{
		parameter_names.at( 0 ) = "mass";
		parameter_names.at( 1 ) = "z";
	}
	catch ( const std::exception & )
	{
		return errorNOS();
	}
	return 0;
}

const int brgastro::point_mass_profile::truncate_to_fraction( const double f,
		const bool silent )
{
	if ( ( f <= 0 ) || ( isbad( f ) ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Bad, negative, or zero f passed to truncate_to_fraction.\n"
					<< "Will truncate to zero instead.\n";
		_mass_ = 0;
	}
	else
	{
		_mass_ *= f;
	}
	return 0;

}

#endif // end point_mass profile functions
