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

#include "../brg_global.h"

#include "../brg_functor.hpp"
#include "../brg_units.h"
#include "density_profile.h"

namespace brgastro {

class accel_functor: public functor< BRG_UNITS >
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

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	accel_functor();
	accel_functor( const density_profile *init_host_ptr );
	virtual ~accel_functor()
	{
	}
};
// class accel_functor

class solve_rhm_functor: public functor< BRG_UNITS >
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

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}
	const int set_target_mass( const BRG_MASS &new_target_mass );
	const BRG_MASS & target_mass()
	{
		return _target_mass_;
	}

	const int operator ()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	solve_rhm_functor();
	solve_rhm_functor( const density_profile *init_host,
			const BRG_MASS &init_target_mass );

};
// end class unitless_solve_rhm_functor

class spherical_density_functor: public functor< BRG_UNITS >
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

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	spherical_density_functor();
	spherical_density_functor( const density_profile *init_host );
	virtual ~spherical_density_functor()
	{
	}
};

class projected_density_functor: public functor< BRG_UNITS >
{
	/**********************************
	 projected_density_functor class
	 -----------------------------

	 Function class integrating density along a projected line

	 Parent class: function_class (from brg_functors)

	 **********************************/
private:
	const density_profile *_host_ptr_;BRG_UNITS _offset_R_;

public:

	const int set_host_ptr( const density_profile *new_host );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	const int set_offset_R( const BRG_DISTANCE &new_offset_R );
	const BRG_DISTANCE offset_R()
	{
		return _offset_R_;
	}

	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	projected_density_functor();
	projected_density_functor( const density_profile *init_host,
			const BRG_DISTANCE &init_offset_R );
	virtual ~projected_density_functor()
	{
	}
};

class cylindrical_density_functor: public functor< BRG_UNITS >
{
	/**********************************
	 cylindrical_density_functor class
	 -----------------------------

	 Function class integrating density in a cylinder

	 Parent class: function_class (from brg_functors)

	 **********************************/
private:
	const density_profile *_host_ptr_;

public:

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	cylindrical_density_functor();
	cylindrical_density_functor( const density_profile *init_host );
	virtual ~cylindrical_density_functor()
	{
	}
};

class offset_ring_dens_functor: public functor< BRG_UNITS >
{

	const density_profile *_host_ptr_;BRG_DISTANCE _R0_, _R_;
public:

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	const int set_R0( const BRG_DISTANCE &new_R_0 );
	const BRG_DISTANCE & R0()
	{
		return _R0_;
	}

	const int set_R( const BRG_DISTANCE &new_R );
	const BRG_DISTANCE & R()
	{
		return _R_;
	}

	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	offset_ring_dens_functor();
	offset_ring_dens_functor( const density_profile *new_host,
			const BRG_DISTANCE &new_R_0 = 0, const BRG_DISTANCE &new_R = 0 );

};

class offset_circ_dens_functor: public functor< std::vector< BRG_UNITS > >
{
private:
	const density_profile *_host_ptr_;BRG_DISTANCE _R0_;

public:

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	const int set_R0( const BRG_DISTANCE &new_R0 );
	const BRG_DISTANCE & R0()
	{
		return _R0_;
	}

	const int operator()( const std::vector< BRG_UNITS > & in_params,
			std::vector< BRG_UNITS > & out_params,
			const bool silent = false ) const;

	offset_circ_dens_functor();
	offset_circ_dens_functor( const density_profile *new_host,
			const BRG_DISTANCE &new_R_0 = 0 );
};

class offset_WLsig_functor: public functor< BRG_UNITS >
{

private:

	const density_profile *_host_ptr_;BRG_DISTANCE _R_;

public:

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	const int set_R( const BRG_DISTANCE &new_R );
	const BRG_DISTANCE & R()
	{
		return _R_;
	}

	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	offset_WLsig_functor();
	offset_WLsig_functor( const density_profile *init_host,
			const BRG_DISTANCE &new_R = 0 );

};

class offset_WLsig_weight_functor: public functor< BRG_UNITS >
{

private:

	const density_profile *_host_ptr_;
	double _c_;

public:

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	const int set_c( const double new_c );
	const double c()
	{
		return _c_;
	}

	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	offset_WLsig_weight_functor();
	offset_WLsig_weight_functor( const density_profile *new_host,
			const double init_c = -1 );

};

} // namespace brgastro

#endif /* _BRG_DENSITY_PROFILE_FUNCTORS_H_ */
