#include <cstdlib>
#include <iostream>
#include <exception>
#include <string>

#include "brg/global.h"

#include "brg/math/cache/cache.hpp"
#include "brg/math/cache/cache_2d.hpp"

#include "brg/physics/astro/astro.h"
#include "brg/physics/units/unit_obj.h"

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
const double brgastro::dfa_cache::_calculate( const double in_params ) const
{
	return brgastro::integrate_add( 0, in_params );
}

#endif // end brgastro::dfa_cache functions

// brgastro::add_cache class methods
const double brgastro::add_cache::_calculate( const double in_param_1, const double in_param_2) const
{
	return brgastro::integrate_add(in_param_1,in_param_2);
}

const double brgastro::tfa_cache::_calculate( const double in_params ) const
{
	return -brgastro::integrate_ltd( 0, brgastro::zfa( in_params ) ) / c;
}

#endif // end class methods
