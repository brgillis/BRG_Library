/**********************************************************************\
 @file pair_bin.cpp
 ------------------

 Source file for the pair_bin class.

 **********************************************************************

 Copyright (C) 2014  Bryan R. Gillis

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 \**********************************************************************/

#include <limits>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include "brg/math/statistics/effective_count.hpp"
#include "brg/math/statistics/standard_error_of_weighted_mean.hpp"

#include "brg/global.h"

#include "brg/math/calculus/integrate.hpp"
#include "brg/physics/lensing/lens_source_pair.h"
#include "brg/physics/lensing/pair_bin_summary.h"
#include "brg/physics/lensing/magnification/magnification_alpha.h"
#include "brg/physics/lensing/magnification/magnification_functors.h"
#include "brg/physics/lensing/magnification/mag_signal_integral_cache.h"
#include "brg/physics/lensing/magnification/mag_weight_integral_cache.h"
#include "brg/physics/units/unit_obj.h"
#include "brg/utility.hpp"

#include "pair_bin.h"

namespace brgastro {

void pair_bin::_uncache_values()
{
	_mu_hat_cached_value_ = std::numeric_limits<double>::infinity();
	_mu_W_cached_value_ = std::numeric_limits<double>::infinity();
}

pair_bin::pair_bin( CONST_BRG_DISTANCE_REF init_R_min, CONST_BRG_DISTANCE_REF init_R_max,
		CONST_BRG_MASS_REF init_m_min, CONST_BRG_MASS_REF init_m_max,
		double init_z_min, double init_z_max,
		double init_mag_min, double init_mag_max )
:	pair_bin_summary(init_R_min,init_R_max,init_m_min,init_m_max,
		init_z_min,init_z_max,init_mag_min,init_mag_max),
	_mu_hat_cached_value_(std::numeric_limits<double>::infinity()),
	_mu_W_cached_value_(std::numeric_limits<double>::infinity())
{
}

// Adding and clearing data
#if(1)

void pair_bin::add_pair( const lens_source_pair & new_pair)
{
	double shear_weight = new_pair.shear_weight();
	double mag_weight = new_pair.mag_weight();
	_R_values_(new_pair.R_proj(), boost::accumulators::weight = mag_weight);
	_m_values_(new_pair.m_lens(), boost::accumulators::weight = mag_weight);
	_z_values_(new_pair.z_lens(), boost::accumulators::weight = mag_weight);
	_source_z_values_(new_pair.z_source(), boost::accumulators::weight = mag_weight);
	_mag_lens_values_(new_pair.mag_lens(), boost::accumulators::weight = mag_weight);
	_mu_obs_values_(brgastro::magnification_alpha(new_pair.mag_source(),new_pair.z_lens())-1,
			boost::accumulators::weight = mag_weight);
	_delta_Sigma_t_values_(new_pair.delta_Sigma_t(), boost::accumulators::weight = shear_weight);
	_delta_Sigma_x_values_(new_pair.delta_Sigma_x(), boost::accumulators::weight = shear_weight);
	_distinct_lens_ids_.insert(new_pair.lens()->index());
	_uncache_values();
}
void pair_bin::clear()
{
	set_zero(_R_values_);
	set_zero(_m_values_);
	set_zero(_z_values_);
	set_zero(_source_z_values_);
	set_zero(_mag_lens_values_);
	set_zero(_mu_obs_values_);
	set_zero(_delta_Sigma_t_values_);
	set_zero(_delta_Sigma_x_values_);
	set_zero(_distinct_lens_ids_);
	_uncache_values();
	pair_bin_summary::clear();
}

#endif

// Calculations on stored values
#if (1)

BRG_UNITS pair_bin::delta_Sigma_t_mean() const
{
	return boost::accumulators::mean(_delta_Sigma_t_values_);
}
BRG_UNITS pair_bin::delta_Sigma_x_mean() const
{
	return boost::accumulators::mean(_delta_Sigma_x_values_);
}

BRG_UNITS pair_bin::delta_Sigma_t_mean_square() const
{
	return boost::accumulators::moment<2>(_delta_Sigma_t_values_);
}
BRG_UNITS pair_bin::delta_Sigma_x_mean_square() const
{
	return boost::accumulators::moment<2>(_delta_Sigma_x_values_);
}

BRG_UNITS pair_bin::delta_Sigma_t_std() const
{
	return std::sqrt(boost::accumulators::variance(_delta_Sigma_t_values_));
}
BRG_UNITS pair_bin::delta_Sigma_x_std() const
{
	return std::sqrt(boost::accumulators::variance(_delta_Sigma_x_values_));
}

BRG_UNITS pair_bin::delta_Sigma_t_stderr() const
{
	return boost::accumulators::standard_error_of_weighted_mean(_delta_Sigma_t_values_);
}
BRG_UNITS pair_bin::delta_Sigma_x_stderr() const
{
	return boost::accumulators::standard_error_of_weighted_mean(_delta_Sigma_x_values_);
}

#endif

// Magnification calculations
#if(1)

double pair_bin::mu_W() const
{
	if(_mu_W_cached_value_==std::numeric_limits<double>::infinity())
	{
		_mu_W_cached_value_ = area()*mag_weight_integral_cache().get(z_mean());
	}
	return _mu_W_cached_value_;
}

double pair_bin::mu_hat() const
{
	if(_mu_hat_cached_value_==std::numeric_limits<double>::infinity())
	{
		// Not cached, so calculate and cache it

		auto mu_observed = boost::accumulators::weighted_mean(_mu_obs_values_)/mu_W();

		const double mu_base = area()*mag_signal_integral_cache().get(z_mean())/mu_W();

		_mu_hat_cached_value_ = mu_observed+mu_base;
	}
	return _mu_hat_cached_value_;

}

#endif

} // end namespace brgastro
