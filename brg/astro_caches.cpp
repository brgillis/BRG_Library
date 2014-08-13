#include <cstdlib>
#include <iostream>
#include <exception>
#include <string>

#include "brg_global.h"

#include "brg_units.h"
#include "brg_astro.h"
#include "cache/cache.hpp"
#include "cache/cache_2d.hpp"
#include "astro_caches.h"

/** Static Class Initialisation **/
#if (1)

// Initialisation for brgastro::dfa_cache
DEFINE_BRG_CACHE_STATIC_VARS( dfa_cache, 0, 5, 0.01 );

// Initialisation for brgastro::add_cache
DEFINE_BRG_CACHE_2D_STATIC_VARS( add_cache, 0, 5, 0.01,
		                                    0, 5, 0.01);

// Initialisation for brgastro::tfa_cache
DEFINE_BRG_CACHE_STATIC_VARS( tfa_cache, 0.001, 1.02, 0.001 );

#endif // end Static Class Initialisation

/** Class Method Definitions **/
#if (1)

// brgastro::dfa_cache class methods
#if (1)
const int brgastro::dfa_cache::_calculate( const double in_params, double & out_params ) const
{
	try
	{
		out_params = brgastro::integrate_add( 0, in_params );
	}
	catch( const std::exception &e)
	{
		return LOWER_LEVEL_ERROR;
	}

	return 0;
}

#endif // end brgastro::dfa_cache functions

// brgastro::add_cache class methods
const int brgastro::add_cache::_calculate( const double in_param_1, const double in_param_2,
		double  & out_param ) const
{
	out_param = brgastro::integrate_add(in_param_1,in_param_2);
	return 0;
}

const int brgastro::tfa_cache::_calculate( const double in_params, double & out_params ) const
{
	try
	{
		out_params = -brgastro::integrate_ltd( 0, brgastro::zfa( in_params ) ) / c;
	}
	catch(const std::exception &e)
	{
		std::cerr << "ERROR: Could not calculate cache for " << _name_base() << "\n"
				<< "Exception: " << e.what() << "\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

#endif // end class methods
