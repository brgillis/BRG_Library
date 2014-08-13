/**       @file lensing_profile_extension_functors.h
 *
 *     Project: brg
 *        Path: /brg/density_profile/lensing_profile_extension_functors.h
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#ifndef _BRG_lensing_profile_extension_FUNCTORS_H_
#define _BRG_lensing_profile_extension_FUNCTORS_H_

#include "../../brg_global.h"

#include "../../brg_functor.hpp"
#include "../../brg_units.h"
#include "lensing_profile_extension.h"

namespace brgastro {

class projected_density_functor: public functor< BRG_UNITS >
{
	/**********************************
	 projected_density_functor class
	 -----------------------------

	 Function class integrating density along a projected line

	 Parent class: function_class (from brg_functors)

	 **********************************/
private:
	const lensing_profile_extension *_host_ptr_;
	BRG_UNITS _offset_R_;

public:

	const int set_host_ptr( const lensing_profile_extension *new_host );
	const lensing_profile_extension * host_ptr()
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
	projected_density_functor( const lensing_profile_extension *init_host,
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

	 Function class for integrating density in a cylinder

	 Parent class: function_class (from brg_functors)

	 **********************************/
private:
	const lensing_profile_extension *_host_ptr_;

public:

	const int set_host_ptr( const lensing_profile_extension *new_host_ptr );
	const lensing_profile_extension * host_ptr()
	{
		return _host_ptr_;
	}

	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	cylindrical_density_functor();
	cylindrical_density_functor( const lensing_profile_extension *init_host );
	virtual ~cylindrical_density_functor()
	{
	}
};

class offset_ring_dens_functor: public functor< BRG_UNITS >
{

	const lensing_profile_extension *_host_ptr_;
	BRG_DISTANCE _R0_, _R_;

public:

	const int set_host_ptr( const lensing_profile_extension *new_host_ptr );
	const lensing_profile_extension * host_ptr()
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
	offset_ring_dens_functor( const lensing_profile_extension *new_host,
			const BRG_DISTANCE &new_R_0 = 0, const BRG_DISTANCE &new_R = 0 );

};

class offset_circ_dens_functor: public functor< BRG_UNITS >
{
private:
	const lensing_profile_extension *_host_ptr_;
	BRG_DISTANCE _R0_, _R_;

	CONST_BRG_DISTANCE_REF _arc_length_in_circle( CONST_BRG_DISTANCE_REF R2 ) const;

public:

	const int set_host_ptr( const lensing_profile_extension *new_host_ptr );
	const lensing_profile_extension * host_ptr()
	{
		return _host_ptr_;
	}

	const int set_R0( const BRG_DISTANCE &new_R0 );
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
			BRG_UNITS & out_param,
			const bool silent = false ) const;

	offset_circ_dens_functor();
	offset_circ_dens_functor( const lensing_profile_extension *new_host,
			const BRG_DISTANCE & new_R_0 = 0, const BRG_DISTANCE &new_R = 0  );
};

class offset_WLsig_functor: public functor< BRG_UNITS >
{

private:

	const lensing_profile_extension *_host_ptr_;BRG_DISTANCE _R_;

public:

	const int set_host_ptr( const lensing_profile_extension *new_host_ptr );
	const lensing_profile_extension * host_ptr()
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
	offset_WLsig_functor( const lensing_profile_extension *init_host,
			const BRG_DISTANCE &new_R = 0 );

};

class quick_offset_WLsig_functor: public functor< BRG_UNITS >
{

private:

	const lensing_profile_extension *_host_ptr_;
	BRG_DISTANCE _R_;

public:

	const int set_host_ptr( const lensing_profile_extension *new_host_ptr );
	const lensing_profile_extension * host_ptr()
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

	quick_offset_WLsig_functor();
	quick_offset_WLsig_functor( const lensing_profile_extension *init_host,
			const BRG_DISTANCE &new_R = 0 );

};

class offset_WLsig_weight_functor: public functor< BRG_UNITS >
{

private:

	const lensing_profile_extension *_host_ptr_;
	double _c_;

public:

	const int set_host_ptr( const lensing_profile_extension *new_host_ptr );
	const lensing_profile_extension * host_ptr()
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
	offset_WLsig_weight_functor( const lensing_profile_extension *new_host,
			const double init_c = -1 );

};

} // namespace brgastro

#endif /* _BRG_lensing_profile_extension_FUNCTORS_H_ */
