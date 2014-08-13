/**********************************************************************\
brg_astro_caches.h
 -----------

 If this header is used, the source file astro_caches.cpp must be included
 and compiled with the project.

 \**********************************************************************/

#ifndef __BRG_ASTRO_CACHES_H_INCLUDED__
#define __BRG_ASTRO_CACHES_H_INCLUDED__

#include <string>

#include "brg_global.h"

#include "cache/brg_cache.hpp"
#include "cache/brg_cache_2d.hpp"
#include "brg_units.h"

namespace brgastro
{

class dfa_cache : public brg_cache<dfa_cache>
{
	// "Distance from angle" cache
private:

	DECLARE_BRG_CACHE_STATIC_VARS();

	friend class brg_cache<dfa_cache>;

protected:

	const std::string _name_base() const throw()
	{
		return "dfa";
	}

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0,0,1,0,0,0,0);
	}
	const brgastro::unit_obj _inverse_units() const throw()
	{
		return brgastro::unit_obj(0);
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function.
	const int _calculate( const double in_params, double & out_params ) const;

public:

};
// class dfa_cache

class add_cache : public brg_cache_2d<add_cache>
{
	// Angular diameter distance
private:

	DECLARE_BRG_CACHE_2D_STATIC_VARS();

	friend class brg_cache_2d<add_cache>;

protected:

	const std::string _name_base() const throw()
	{
		char name_base[BRG_CACHE_ND_NAME_SIZE] = "ang_di_d";
		return name_base;
	}

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0,1,0,0,0,0,0); // Distance units
	}
	const brgastro::unit_obj _inverse_units() const throw()
	{
		return brgastro::unit_obj(0); // Unitless (redshift)
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function.
	const int _calculate( const double in_param_1, const double in_param_2, double & out_param ) const;

public:

};
// class add_cache

class tfa_cache : public brg_cache<tfa_cache>
{
	// "Time from a (scale factor)" cache
private:

	DECLARE_BRG_CACHE_STATIC_VARS();

	friend class brg_cache<tfa_cache>;

protected:

	const std::string _name_base() const throw()
	{
		return "tfa";
	}

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0,0,1,0,0,0,0);
	}
	const brgastro::unit_obj _inverse_units() const throw()
	{
		return brgastro::unit_obj(0);
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function.
	const int _calculate( const double in_params, double & out_params ) const;

public:

};
// class tfa_cache

} // end namespace brgastro

#endif // __BRG_ASTRO_CACHES_H_INCLUDED__

