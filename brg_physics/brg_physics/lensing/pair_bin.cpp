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

#include "brg/global.h"

#include "brg/math/calculus/integrate.hpp"
#include "brg/math/statistics/statistic_extractors.hpp"
#include "brg/physics/lensing/lens_source_pair.h"
#include "brg/physics/lensing/pair_bin_summary.h"
#include "brg/physics/lensing/magnification/mag_global_values.h"
#include "brg/physics/lensing/magnification/magnification_alpha.h"
#include "brg/physics/lensing/magnification/magnification_functors.h"
#include "brg/physics/lensing/magnification/mag_signal_integral_cache.h"
#include "brg/physics/lensing/magnification/mag_weight_integral_cache.h"
#include "brg/math/statistics/effective_count.hpp"
#include "brg/physics/units/unit_obj.h"
#include "brg/utility.hpp"

#include "pair_bin.h"

#include "../../math/statistics/error_of_weighted_mean.hpp"

namespace brgastro {

constexpr double pair_bin::_z_buffer_default_value_;

void pair_bin::_uncache_values()
{
	_mu_hat_cached_value_ = std::numeric_limits<double>::infinity();
	_mu_W_cached_value_ = std::numeric_limits<double>::infinity();
}

pair_bin::pair_bin( CONST_BRG_DISTANCE_REF init_R_min, CONST_BRG_DISTANCE_REF init_R_max,
		CONST_BRG_MASS_REF init_m_min, CONST_BRG_MASS_REF init_m_max,
		const double & init_z_min, const double & init_z_max,
		const double & init_mag_min, const double & init_mag_max,
		const double & init_z_buffer)
:	pair_bin_summary(init_R_min,init_R_max,init_m_min,init_m_max,
		init_z_min,init_z_max,init_mag_min,init_mag_max),
	_mu_hat_cached_value_(std::numeric_limits<double>::infinity()),
	_mu_W_cached_value_(std::numeric_limits<double>::infinity()),
	_z_buffer_(init_z_buffer)
{
}

// Adding and clearing data
#if(1)

void pair_bin::add_pair( const lens_source_pair & new_pair)
{
	const double & shear_weight = new_pair.shear_weight();
	const double & mag_weight = new_pair.mag_weight();

	// Shear first
	if(shear_weight>0)
	{
		// General info
		_shear_R_values_(new_pair.R_proj(), boost::accumulators::weight = shear_weight);

		_shear_lens_m_values_(new_pair.m_lens(), boost::accumulators::weight = shear_weight);
		_shear_lens_mag_values_(new_pair.mag_lens(), boost::accumulators::weight = shear_weight);
		_shear_lens_z_values_(new_pair.z_lens(), boost::accumulators::weight = shear_weight);
		_shear_source_z_values_(new_pair.z_source(), boost::accumulators::weight = shear_weight);

		// Shear info
		_delta_Sigma_t_values_(new_pair.delta_Sigma_t(), boost::accumulators::weight = shear_weight);
		_delta_Sigma_x_values_(new_pair.delta_Sigma_x(), boost::accumulators::weight = shear_weight);
	}

	// Now magnification
	const auto & mag_source = new_pair.mag_source();
	const auto & z_source = new_pair.z_source();
	if((mag_weight>0) && (mag_source>=mag_m_min) && (mag_source<=mag_m_max) &&
		(z_source>=mag_z_min) && (z_source<=mag_z_max) )
	{
		// General info
		_magf_R_values_(new_pair.R_proj(), boost::accumulators::weight = mag_weight);
		_magf_source_z_values_(new_pair.z_source(), boost::accumulators::weight = mag_weight);

		// Magnification info
		_mu_obs_values_(brgastro::magnification_alpha(new_pair.mag_source(),new_pair.z_lens()+_z_buffer_)-1,
				boost::accumulators::weight = mag_weight);
		_uncache_values();
	}
}
void pair_bin::add_lens( const lens_id & lens )
{
	if(_distinct_lens_ids_.find(lens.id)==_distinct_lens_ids_.end())
	{
		_distinct_lens_ids_.insert(lens.id);
		_magf_lens_m_values_(lens.m, boost::accumulators::weight = lens.weight);
		_magf_lens_mag_values_(lens.mag, boost::accumulators::weight = lens.weight);
		_magf_lens_z_values_(lens.z, boost::accumulators::weight = lens.weight);
		_magf_unmasked_fracs_(lens.unmasked_frac(R_amid()), boost::accumulators::weight = lens.weight);
		_uncache_values();
	}
}
void pair_bin::clear()
{
	set_zero(_magf_R_values_);
	set_zero(_shear_R_values_);
	set_zero(_magf_lens_m_values_);
	set_zero(_shear_lens_m_values_);
	set_zero(_magf_lens_z_values_);
	set_zero(_shear_lens_z_values_);
	set_zero(_magf_lens_mag_values_);
	set_zero(_shear_lens_mag_values_);

	set_zero(_magf_source_z_values_);
	set_zero(_shear_source_z_values_);

	set_zero(_magf_unmasked_fracs_);
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

BRG_UNITS pair_bin::area() const
{
	double mean_frac = unmasked_frac();
	BRG_UNITS mean_area = pi*(square(afd(R_max(),magf_lens_z_mean()))-square(afd(R_min(),magf_lens_z_mean())));

	BRG_UNITS result = mean_frac*magf_num_lenses()*mean_area;
	brgastro::fixbad(result);
	return result;
}

BRG_UNITS pair_bin::delta_Sigma_t_mean() const
{
	return safe_extract_mean(_delta_Sigma_t_values_);
}
BRG_UNITS pair_bin::delta_Sigma_x_mean() const
{
	return safe_extract_mean(_delta_Sigma_x_values_);
}

BRG_UNITS pair_bin::delta_Sigma_t_mean_square() const
{
	return safe_extract_moment<2>(_delta_Sigma_t_values_);
}
BRG_UNITS pair_bin::delta_Sigma_x_mean_square() const
{
	return safe_extract_moment<2>(_delta_Sigma_x_values_);
}

BRG_UNITS pair_bin::delta_Sigma_t_std() const
{
	return std::sqrt(safe_extract_variance(_delta_Sigma_t_values_));
}
BRG_UNITS pair_bin::delta_Sigma_x_std() const
{
	return std::sqrt(safe_extract_variance(_delta_Sigma_x_values_));
}

BRG_UNITS pair_bin::delta_Sigma_t_stderr() const
{
	return safe_extract_error_of_weighted_mean(_delta_Sigma_t_values_);
}
BRG_UNITS pair_bin::delta_Sigma_x_stderr() const
{
	return safe_extract_error_of_weighted_mean(_delta_Sigma_x_values_);
}

#endif

// Magnification calculations
#if(1)

double pair_bin::mu_W() const
{
	if(_mu_W_cached_value_==std::numeric_limits<double>::infinity())
	{
		_mu_W_cached_value_ = area()*mag_weight_integral_cache().get(magf_lens_z_mean()+_z_buffer_);
	}
	return _mu_W_cached_value_;
}

double pair_bin::mu_hat() const
{
	if(_mu_hat_cached_value_==std::numeric_limits<double>::infinity())
	{
		// Not cached, so calculate and cache it

		const double mu_observed = extract_weighted_sum(_mu_obs_values_)/mu_W();

		const double mu_base = area()*mag_signal_integral_cache().get(magf_lens_z_mean()+_z_buffer_)/mu_W();

		_mu_hat_cached_value_ = mu_observed+mu_base;

		brgastro::fixbad(_mu_hat_cached_value_);
	}
	return _mu_hat_cached_value_;

}

#endif

} // end namespace brgastro
