/**       @file tNFW_profile_functors.h
 *
 *     Project: brg
 *        Path: /brg/density_profile/tNFW_profile_functors.h
 *
 *  Created on: 19 Aug 2014
 *      Author: brg
 */

#ifndef _BRG_TNFW_PROFILE_FUNCTORS_HPP_
#define _BRG_TNFW_PROFILE_FUNCTORS_HPP_

#include <cstdlib>

#include "brg/brg_global.h"

#include "density_profile.h"
#include "tNFW_profile.h"

namespace brgastro
{

class tNFW_solve_rvir_iterative_functor
{
private:
	const tNFW_profile *_halo_;

public:

	tNFW_solve_rvir_iterative_functor(const tNFW_profile *halo)
	: _halo_(halo)
	{
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF in_param,	const bool silent = false ) const
	{
		if(in_param<=0)
		{
			return _halo_->rvir0();
		}
		else
		{
			return in_param * _halo_->enc_dens(in_param)/(virial_density_factor*_halo_->rho_crit());
		}
	}

};

class tNFW_solve_rvir_minimize_functor
{
private:
	const tNFW_profile *_halo_;

public:

	tNFW_solve_rvir_minimize_functor(const tNFW_profile *halo)
	: _halo_(halo)
	{
	}

	BRG_UNITS operator()( CONST_BRG_UNITS_REF  in_param, const bool silent = false ) const
	{
		return std::fabs(virial_density_factor*_halo_->rho_crit()-_halo_->enc_dens(in_param));
	}

};

} // end namespace brgastro

#endif /* _BRG_TNFW_PROFILE_FUNCTORS_HPP_ */
