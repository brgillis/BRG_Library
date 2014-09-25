/**********************************************************************\
 @file pair_bin_summary.cpp
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

#include <cassert>

#include "brg/global.h"

#include "brg/math/calculus/integrate.hpp"
#include "brg/physics/lensing/magnification/magnification_alpha.h"
#include "brg/physics/lensing/magnification/magnification_functors.h"
#include "brg/physics/lensing/pair_bin.h"
#include "brg/physics/units/unit_obj.h"
#include "brg/vector/summary_functions.hpp"

#include "pair_bin_summary.h"

namespace brgastro {


pair_bin_summary::pair_bin_summary( CONST_BRG_DISTANCE_REF init_R_min, CONST_BRG_DISTANCE_REF init_R_max,
		CONST_BRG_MASS_REF init_m_min, CONST_BRG_MASS_REF init_m_max,
		double init_z_min, double init_z_max,
		double init_mag_min, double init_mag_max )
:	_R_min_(init_R_min),
	_R_max_(init_R_max),
	_m_min_(init_m_min),
	_m_max_(init_m_max),
	_z_min_(init_z_min),
	_z_max_(init_z_max),
	_mag_min_(init_mag_min),
	_mag_max_(init_mag_max),
	_count_(0),
	_sum_of_weights_(0),
	_sum_of_square_weights_(0),
	_R_mean_(0),
	_m_mean_(0),
	_z_mean_(0),
	_mag_mean_(0),
	_delta_Sigma_t_mean_(0),
	_delta_Sigma_x_mean_(0),
	_delta_Sigma_t_mean_square_(0),
	_delta_Sigma_x_mean_square_(0),
	_mu_hat_cached_value_(std::numeric_limits<double>::infinity()),
	_mu_W_cached_value_(std::numeric_limits<double>::infinity()),
	_num_lenses_(0)
{
}

pair_bin_summary::pair_bin_summary( const pair_bin & bin)
: 	_R_min_(bin.R_min()),
	_R_max_(bin.R_max()),
	_m_min_(bin.m_min()),
	_m_max_(bin.m_max()),
	_z_min_(bin.z_min()),
	_z_max_(bin.z_max()),
	_mag_min_(bin.mag_min()),
	_mag_max_(bin.mag_max()),
	_count_(bin.count()),
	_sum_of_weights_(bin.sum_of_weights()),
	_sum_of_square_weights_(bin.sum_of_square_weights()),
	_R_mean_(bin.R_mean()),
	_m_mean_(bin.m_mean()),
	_z_mean_(bin.z_mean()),
	_mag_mean_(bin.mag_mean()),
	_delta_Sigma_t_mean_(bin.delta_Sigma_t_mean()),
	_delta_Sigma_x_mean_(bin.delta_Sigma_x_mean()),
	_delta_Sigma_t_mean_square_(bin.delta_Sigma_t_mean_square()),
	_delta_Sigma_x_mean_square_(bin.delta_Sigma_x_mean_square()),
	_mu_hat_cached_value_(std::numeric_limits<double>::infinity()),
	_mu_W_cached_value_(std::numeric_limits<double>::infinity()),
	_num_lenses_(bin.num_lenses())
{
}

// Adding and clearing data
#if(1)

void pair_bin_summary::clear()
{
	_count_ = 0;
	_sum_of_weights_ = 0;
	_sum_of_square_weights_ = 0;
	_R_mean_ = 0;
	_m_mean_ = 0;
	_z_mean_ = 0;
	_mag_mean_ = 0;
	_delta_Sigma_t_mean_ = 0;
	_delta_Sigma_x_mean_ = 0;
	_delta_Sigma_t_mean_square_ = 0;
	_delta_Sigma_x_mean_square_ = 0;
	_num_lenses_ = 0;
	_uncache_values();
}

#endif // Adding and clearing data

// Magnification calculations
#if(1)

double pair_bin_summary::mu_W() const
{
	if(_mu_W_cached_value_==std::numeric_limits<double>::infinity())
	{
		// Not cached, so calculate and cache it

		mu_weight_integration_functor func(area(),z_mean());

		_mu_W_cached_value_ = integrate_Romberg(&func,15,25);
	}
	return _mu_W_cached_value_;

}

double pair_bin_summary::mu_hat() const
{
	if(_mu_hat_cached_value_==std::numeric_limits<double>::infinity())
	{
		// Not cached, so calculate and cache it

		double mu_observed = 0;
		for(size_t source_i=0;source_i<_source_magnitudes_.size();++source_i)
		{
			mu_observed += magnification_alpha(_source_magnitudes_[source_i])-1;
		}

		mu_signal_integration_functor func(area(),z_mean());
		const double mu_base = integrate_Romberg(&func,15,25);

		_mu_hat_cached_value_ = (mu_observed-mu_base)/mu_W();
	}
	return _mu_hat_cached_value_;

}

#endif

// Combining summaries together
#if(1)

pair_bin_summary & pair_bin_summary::operator+=( const pair_bin_summary & other )
{
	// Check bin limits are all the same
	assert(_R_min_==other.R_min());
	assert(_R_max_==other.R_max());
	assert(_m_min_==other.m_min());
	assert(_m_max_==other.m_max());
	assert(_z_min_==other.z_min());
	assert(_z_max_==other.z_max());
	assert(_mag_min_==other.mag_min());
	assert(_mag_max_==other.mag_max());

	// Check for zero count in this or the other
	if(other.count()==0) return *this;

	_uncache_values(); // If we get this far, something is changing with this bin

	if(_count_==0)
	{
		// Simply copy the other bin summary
		*this = other;
		return *this;
	}

	// Combine the source magnitude lists together
	_source_magnitudes_.insert( _source_magnitudes_.end(),
			other.source_magnitudes().begin(), other.source_magnitudes().end() );

	// Add the count and number of distinct lenses
	_count_ += other.count();
	_num_lenses_ += other.num_lenses();

	// Check for zero total weight in this or the other
	if(other.sum_of_weights()==0)
	{
		return *this;
	}
	if(_sum_of_weights_==0)
	{
		// Simply copy the other bin summary, except for the list of source magnitudes
		auto tmp_source_magnitudes = std::move(_source_magnitudes_);
		auto tmp_count = std::move(_count_);
		auto tmp_num_lenses = std::move(_num_lenses_);
		*this = other;
		_source_magnitudes_ = std::move(tmp_source_magnitudes);
		_count_ = std::move(tmp_count);
		_num_lenses_ = std::move(tmp_num_lenses);
		return *this;
	}

	double new_sum_of_weights = sum_of_weights() + other.sum_of_weights();

	_R_mean_ = ( _R_mean_*sum_of_weights() + other.R_mean()*other.sum_of_weights())/new_sum_of_weights;
	_m_mean_ = ( _m_mean_*sum_of_weights() + other.m_mean()*other.sum_of_weights())/new_sum_of_weights;
	_z_mean_ = ( _z_mean_*sum_of_weights() + other.z_mean()*other.sum_of_weights())/new_sum_of_weights;
	_mag_mean_ = ( _mag_mean_*sum_of_weights() + other.mag_mean()*other.sum_of_weights())/new_sum_of_weights;

	_delta_Sigma_t_mean_ = ( _delta_Sigma_t_mean_*sum_of_weights() + other.delta_Sigma_t_mean()*other.sum_of_weights())
			/new_sum_of_weights;
	_delta_Sigma_t_mean_square_ = ( _delta_Sigma_t_mean_square_*sum_of_weights() +
			other.delta_Sigma_t_mean_square()*other.sum_of_weights())
			/new_sum_of_weights;
	_delta_Sigma_x_mean_ = ( _delta_Sigma_x_mean_*sum_of_weights() + other.delta_Sigma_x_mean()*other.sum_of_weights())
					/new_sum_of_weights;
	_delta_Sigma_x_mean_square_ = ( _delta_Sigma_x_mean_square_*sum_of_weights() +
			other.delta_Sigma_x_mean_square()*other.sum_of_weights())
			/new_sum_of_weights;
	_sum_of_weights_ = new_sum_of_weights;
	_sum_of_square_weights_ += other.sum_of_square_weights();

	return *this;
}

#endif // Combining summaries together

} // end namespace brgastro
