/**       @file density_profile_functors.h
 *
 *     Project: brg
 *        Path: /brg/density_profile/density_profile_functors.h
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#ifndef _BRG_DENSITY_PROFILE_FUNCTORS_H_
#define _BRG_DENSITY_PROFILE_FUNCTORS_H_

#include "brg/brg_global.h"

#include "brg/brg_units.h"
#include "brg/density_profile/density_profile.h"

namespace brgastro {

class accel_functor
{
	/**********************************
	 accel_functor class
	 -----------------------------

	 Function class for acceleration within a density profile

	 Parent class: function_class (from brg_functors)

	 **********************************/
private:
	const density_profile *_host_ptr_;
public:

	void set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param, const bool silent = false ) const;

	accel_functor();
	accel_functor( const density_profile *init_host_ptr );
	virtual ~accel_functor()
	{
	}
};
// class accel_functor

class solve_rhm_functor
{
	/**********************************
	 solve_rhm_functor class
	 -----------------------------

	 Function class for solving the half-mass
	 radius of a halo.

	 Parent class: function_class (from brg_functors)

	 **********************************/
private:
	const density_profile *_host_ptr_;BRG_MASS _target_mass_;

public:

	void set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}
	void set_target_mass( const BRG_MASS &new_target_mass );
	const BRG_MASS & target_mass()
	{
		return _target_mass_;
	}

	BRG_UNITS operator ()( CONST_BRG_UNITS_REF in_param, const bool silent = false ) const;

	solve_rhm_functor();
	solve_rhm_functor( const density_profile *init_host,
			const BRG_MASS &init_target_mass );

};
// end class unitless_solve_rhm_functor

class spherical_density_functor
{
	/**********************************
	 spherical_density_functor class
	 -----------------------------

	 Function class integrating density in a sphere

	 Parent class: function_class (from brg_functors)

	 **********************************/
private:
	const density_profile *_host_ptr_;

public:

	void set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF in_param, const bool silent = false ) const;

	spherical_density_functor();
	spherical_density_functor( const density_profile *init_host );
	virtual ~spherical_density_functor()
	{
	}
};

} // namespace brgastro

#endif /* _BRG_DENSITY_PROFILE_FUNCTORS_H_ */
