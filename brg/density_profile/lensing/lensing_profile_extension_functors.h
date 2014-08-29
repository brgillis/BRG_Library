/**       @file lensing_profile_extension_functors.h
 *
 *     Project: brg
 *        Path: /brg/density_profile/lensing_profile_extension_functors.h
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#ifndef _BRG_LENSING_PROFILE_EXTENSION_FUNCTORS_H_
#define _BRG_LENSING_PROFILE_EXTENSION_FUNCTORS_H_

#include "brg/brg_global.h"

#include "brg/brg_units.h"
#include "brg/density_profile/lensing/lensing_profile_extension.h"

namespace brgastro {

class projected_density_functor
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

	void set_host_ptr( const lensing_profile_extension *new_host );
	const lensing_profile_extension * host_ptr()
	{
		return _host_ptr_;
	}

	void set_offset_R( CONST_BRG_DISTANCE_REF new_offset_R );
	const BRG_DISTANCE offset_R()
	{
		return _offset_R_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF in_param,
			const bool silent = false ) const;

	projected_density_functor();
	projected_density_functor( const lensing_profile_extension *init_host,
			CONST_BRG_DISTANCE_REF init_offset_R );
	virtual ~projected_density_functor()
	{
	}
};

class cylindrical_density_functor
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

	void set_host_ptr( const lensing_profile_extension *new_host_ptr );
	const lensing_profile_extension * host_ptr()
	{
		return _host_ptr_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF in_param, const bool silent = false ) const;

	cylindrical_density_functor();
	cylindrical_density_functor( const lensing_profile_extension *init_host );
	virtual ~cylindrical_density_functor()
	{
	}
};

class offset_ring_dens_functor
{

	const lensing_profile_extension *_host_ptr_;
	BRG_DISTANCE _R0_, _R_;

public:

	void set_host_ptr( const lensing_profile_extension *new_host_ptr );
	const lensing_profile_extension * host_ptr()
	{
		return _host_ptr_;
	}

	void set_R0( CONST_BRG_DISTANCE_REF new_R_0 );
	CONST_BRG_DISTANCE_REF  R0()
	{
		return _R0_;
	}

	void set_R( CONST_BRG_DISTANCE_REF new_R );
	CONST_BRG_DISTANCE_REF  R()
	{
		return _R_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param,
			const bool silent = false ) const;

	offset_ring_dens_functor();
	offset_ring_dens_functor( const lensing_profile_extension *new_host,
			CONST_BRG_DISTANCE_REF new_R_0 = 0, CONST_BRG_DISTANCE_REF new_R = 0 );

};

class offset_circ_dens_functor
{
private:
	const lensing_profile_extension *_host_ptr_;
	BRG_DISTANCE _R0_, _R_;

	CONST_BRG_DISTANCE_REF _arc_length_in_circle( CONST_BRG_DISTANCE_REF R2 ) const;

public:

	void set_host_ptr( const lensing_profile_extension *new_host_ptr );
	const lensing_profile_extension * host_ptr()
	{
		return _host_ptr_;
	}

	void set_R0( CONST_BRG_DISTANCE_REF new_R0 );
	CONST_BRG_DISTANCE_REF R0()
	{
		return _R0_;
	}

	void set_R( CONST_BRG_DISTANCE_REF new_R );
	CONST_BRG_DISTANCE_REF R()
	{
		return _R_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF in_param,
			const bool silent = false ) const;

	offset_circ_dens_functor();
	offset_circ_dens_functor( const lensing_profile_extension *new_host,
			CONST_BRG_DISTANCE_REF new_R_0 = 0, CONST_BRG_DISTANCE_REF new_R = 0  );
};

class offset_WLsig_functor
{

private:

	const lensing_profile_extension *_host_ptr_;
	BRG_DISTANCE _R_;

public:

	void set_host_ptr( const lensing_profile_extension *new_host_ptr );
	const lensing_profile_extension * host_ptr()
	{
		return _host_ptr_;
	}

	void set_R( CONST_BRG_DISTANCE_REF new_R );
	CONST_BRG_DISTANCE_REF  R()
	{
		return _R_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param,
			const bool silent = false ) const;

	offset_WLsig_functor();
	offset_WLsig_functor( const lensing_profile_extension *init_host,
			CONST_BRG_DISTANCE_REF new_R = 0 );

};

class quick_offset_WLsig_functor
{

private:

	const lensing_profile_extension *_host_ptr_;
	BRG_DISTANCE _R_;

public:

	void set_host_ptr( const lensing_profile_extension *new_host_ptr );
	const lensing_profile_extension * host_ptr()
	{
		return _host_ptr_;
	}

	void set_R( CONST_BRG_DISTANCE_REF new_R );
	CONST_BRG_DISTANCE_REF  R()
	{
		return _R_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param,
			const bool silent = false ) const;

	quick_offset_WLsig_functor();
	quick_offset_WLsig_functor( const lensing_profile_extension *init_host,
			CONST_BRG_DISTANCE_REF new_R = 0 );

};

class group_WLsig_weight_functor
{

private:

	const lensing_profile_extension *_host_ptr_;
	double _c_;

public:

	void set_host_ptr( const lensing_profile_extension *new_host_ptr );
	const lensing_profile_extension * host_ptr()
	{
		return _host_ptr_;
	}

	void set_c( const double new_c );
	const double c()
	{
		return _c_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param,
			const bool silent = false ) const;

	group_WLsig_weight_functor();
	group_WLsig_weight_functor( const lensing_profile_extension *new_host,
			const double init_c = -1 );

};

class shifted_WLsig_weight_functor
{

private:

	BRG_DISTANCE _sigma_;

public:

	void set_sigma( CONST_BRG_DISTANCE_REF new_sigma );
	CONST_BRG_DISTANCE_REF sigma()
	{
		return _sigma_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param,
			const bool silent = false ) const;

	shifted_WLsig_weight_functor();
	shifted_WLsig_weight_functor( CONST_BRG_DISTANCE_REF new_sigma );

};

class shifted_WLsig_circ_functor
{

private:

	const lensing_profile_extension *_host_ptr_;
	BRG_DISTANCE _R_, _R_shift_;

public:

	void set_host_ptr( const lensing_profile_extension *new_host_ptr );
	const lensing_profile_extension * host_ptr()
	{
		return _host_ptr_;
	}

	void set_R_shift( CONST_BRG_DISTANCE_REF new_R_shift );
	CONST_BRG_DISTANCE_REF R_shift()
	{
		return _R_shift_;
	}

	void set_R( CONST_BRG_DISTANCE_REF new_R );
	CONST_BRG_DISTANCE_REF R()
	{
		return _R_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param, const bool silent = false ) const;

	shifted_WLsig_circ_functor( const lensing_profile_extension *new_host,
			CONST_BRG_DISTANCE_REF new_R_shift, CONST_BRG_DISTANCE_REF new_R );

};

class shifted_WLsig_functor
{

private:

	const lensing_profile_extension *_host_ptr_;
	BRG_DISTANCE _R_;

public:

	void set_host_ptr( const lensing_profile_extension *new_host_ptr );
	const lensing_profile_extension * host_ptr()
	{
		return _host_ptr_;
	}

	void set_R( CONST_BRG_DISTANCE_REF new_R );
	CONST_BRG_DISTANCE_REF R()
	{
		return _R_;
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF in_param, const bool silent = false ) const;

	shifted_WLsig_functor( const lensing_profile_extension *new_host,
			CONST_BRG_DISTANCE_REF new_R );

};

} // namespace brgastro

#endif /* _BRG_lensing_profile_extension_FUNCTORS_H_ */
