/**********************************************************************\
 @file pair_bin.h
 ------------------

 A class representing a bin of lens-source pairs, which can be used for
 calculating statistics on the contained pairs.

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

// body file: pair_bin.cpp
#ifndef _BRG_PAIR_BIN_H_INCLUDED_
#define _BRG_PAIR_BIN_H_INCLUDED_

#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/accumulators/statistics/weighted_moment.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/container/flat_set.hpp>

#include "brg/math/statistics/effective_count.hpp"
#include "brg/math/statistics/statistic_extractors.hpp"

#include "brg/global.h"

#include "brg/physics/lensing/lens_source_pair.h"
#include "brg/physics/lensing/pair_bin_summary.h"
#include "brg/physics/units/unit_obj.h"
#include "brg/vector/limit_vector.hpp"
#include "brg/vector/summary_functions.hpp"
#include "../../math/statistics/error_of_weighted_mean.hpp"

namespace brgastro {

// lens_id struct
#if(1)
struct lens_id
{
	size_t id;
	BRG_MASS m;
	double z;
	double mag;
	double weight;

	// Data on the unmasked fraction of annuli
#if(1)

	brgastro::limit_vector<BRG_DISTANCE> unmasked_frac_bin_limits;
	std::vector<double> unmasked_fracs;

	double unmasked_frac(const BRG_DISTANCE & R_proj) const;

#endif

	lens_id(const size_t & id, const BRG_MASS & m, const double & z, const double & mag,
			const std::vector<BRG_DISTANCE> & unmasked_frac_bin_limits,
			const std::vector<double> & unmasked_fracs,
			const double & weight=1)
	: id(id),
	  m(m),
	  z(z),
	  mag(mag),
	  weight(weight),
	  unmasked_frac_bin_limits(unmasked_frac_bin_limits),
	  unmasked_fracs(unmasked_fracs)
	{
	}
};
#endif

/**
 *
 */
class pair_bin: public pair_bin_summary {
private:

	template <typename T>
	using bin_stat_vec_t = boost::accumulators::accumulator_set<T,
			boost::accumulators::stats<
				boost::accumulators::tag::count,
				boost::accumulators::tag::weighted_mean >,
				double >;
	template <typename T>
	using stat_vec_t = boost::accumulators::accumulator_set<T,
			boost::accumulators::stats<
				boost::accumulators::tag::count,
				boost::accumulators::tag::effective_count,
				boost::accumulators::tag::weighted_mean,
				boost::accumulators::tag::weighted_moment<2>,
				boost::accumulators::tag::weighted_variance,
				boost::accumulators::tag::error_of_weighted_mean >,
				double >;

	// Pair data
#if(1)
	bin_stat_vec_t<BRG_DISTANCE> _magf_R_values_;
	bin_stat_vec_t<BRG_DISTANCE> _shear_R_values_;

	bin_stat_vec_t<BRG_MASS> _magf_lens_m_values_;
	bin_stat_vec_t<BRG_MASS> _shear_lens_m_values_;
	bin_stat_vec_t<double> _magf_lens_z_values_;
	bin_stat_vec_t<double> _shear_lens_z_values_;
	bin_stat_vec_t<double> _magf_lens_mag_values_;
	bin_stat_vec_t<double> _shear_lens_mag_values_;

	bin_stat_vec_t<double> _magf_source_z_values_;
	bin_stat_vec_t<double> _shear_source_z_values_;

	bin_stat_vec_t<double> _magf_unmasked_fracs_;
	bin_stat_vec_t<double> _mu_obs_values_;

	stat_vec_t<BRG_UNITS> _delta_Sigma_t_values_;
	stat_vec_t<BRG_UNITS> _delta_Sigma_x_values_;

	boost::container::flat_set<size_t> _distinct_lens_ids_;

	mutable double _mu_hat_cached_value_;
	mutable double _mu_W_cached_value_;

	double _z_buffer_;
	static constexpr double _z_buffer_default_value_ = 0.1;

#endif // Pair data

protected:
	void _uncache_values();

public:

	// Constructors and destructor
#if(1)
	pair_bin( CONST_BRG_DISTANCE_REF init_R_min=0, CONST_BRG_DISTANCE_REF init_R_max=0,
			CONST_BRG_MASS_REF init_m_min=0, CONST_BRG_MASS_REF init_m_max=0,
			const double & init_z_min=0, const double & init_z_max=0,
			const double & init_mag_min=0, const double & init_mag_max=0,
			const double & init_z_buffer=_z_buffer_default_value_);
	virtual ~pair_bin()
	{
	}
#endif

	// Setting and accessing z_buffer
#if(1)
	void set_z_buffer(const double & new_z_buffer)
	{
		_z_buffer_ = new_z_buffer;
	}
	double z_buffer() const
	{
		return _z_buffer_;
	}
#endif

	// Adding and clearing data
#if(1)

	void add_pair( const lens_source_pair & new_pair);
	void add_lens( const lens_id & lens );
	void clear();

#endif

	// Count
#if(1)
	size_t count() const
	{
		return shear_count();
	}
	size_t shear_count() const
	{
		return extract_count(_shear_R_values_);
	}
	size_t magf_count() const
	{
		return extract_count(_magf_R_values_);
	}
	double effective_count() const
	{
		return shear_effective_count();
	}
	double shear_effective_count() const
	{
		return safe_extract_effective_count(_delta_Sigma_t_values_);
	}
	double sum_of_weights() const
	{
		return shear_sum_of_weights();
	}
	double shear_sum_of_weights() const
	{
		return safe_extract_sum_of_weights(_delta_Sigma_t_values_);
	}
	double sum_of_square_weights() const
	{
		return shear_sum_of_square_weights();
	}
	double shear_sum_of_square_weights() const
	{
		return safe_extract_sum_of_square_weights(_delta_Sigma_t_values_);
	}
	size_t num_lenses() const
	{
		return magf_num_lenses();
	}
	size_t magf_num_lenses() const
	{
		return _distinct_lens_ids_.size();
	}
#endif

	// Limits and means accessors
#if(1)

	BRG_DISTANCE R_mean() const
	{
		return shear_R_mean();
	}
	BRG_DISTANCE shear_R_mean() const
	{
		return safe_extract_weighted_mean(_shear_R_values_);
	}
	BRG_DISTANCE magf_R_mean() const
	{
		return safe_extract_weighted_mean(_magf_R_values_);
	}

	BRG_MASS m_lens_mean() const
	{
		return shear_lens_m_mean();
	}
	BRG_MASS shear_lens_m_mean() const
	{
		return safe_extract_weighted_mean(_shear_lens_m_values_);
	}
	BRG_MASS magf_lens_m_mean() const
	{
		return safe_extract_weighted_mean(_magf_lens_m_values_);
	}

	double lens_z_mean() const
	{
		return shear_lens_z_mean();
	}
	double shear_lens_z_mean() const
	{
		return safe_extract_weighted_mean(_shear_lens_z_values_);
	}
	double magf_lens_z_mean() const
	{
		return safe_extract_weighted_mean(_magf_lens_z_values_);
	}

	double lens_mag_mean() const
	{
		return shear_lens_mag_mean();
	}
	double shear_lens_mag_mean() const
	{
		return safe_extract_weighted_mean(_shear_lens_mag_values_);
	}
	double magf_lens_mag_mean() const
	{
		return safe_extract_weighted_mean(_magf_lens_mag_values_);
	}

	double source_z_mean() const
	{
		return shear_source_z_mean();
	}
	double shear_source_z_mean() const
	{
		return safe_extract_weighted_mean(_shear_source_z_values_);
	}
	double magf_source_z_mean() const
	{
		return safe_extract_weighted_mean(_magf_source_z_values_);
	}

	double unmasked_frac() const
	{
		return safe_extract_weighted_mean(_magf_unmasked_fracs_);
	}

#endif

	// Calculations on stored values
#if (1)

	BRG_UNITS area() const;

	BRG_UNITS delta_Sigma_t_mean() const;
	BRG_UNITS delta_Sigma_x_mean() const;

	BRG_UNITS delta_Sigma_t_mean_square() const;
	BRG_UNITS delta_Sigma_x_mean_square() const;

	BRG_UNITS delta_Sigma_t_std() const;
	BRG_UNITS delta_Sigma_x_std() const;

	BRG_UNITS delta_Sigma_t_stderr() const;
	BRG_UNITS delta_Sigma_x_stderr() const;

	double mu_hat() const;
	double mu_W() const;

#endif
};

} // end namespace brgastro

#endif // _BRG_PAIR_BIN_H_INCLUDED_
