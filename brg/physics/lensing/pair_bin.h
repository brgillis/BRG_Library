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
#include "brg/math/statistics/standard_error_of_weighted_mean.hpp"

#include "brg/global.h"

#include "brg/physics/lensing/lens_source_pair.h"
#include "brg/physics/lensing/pair_bin_summary.h"
#include "brg/physics/units/unit_obj.h"
#include "brg/vector/summary_functions.hpp"

namespace b_acc = boost::accumulators;

namespace brgastro {

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
				boost::accumulators::tag::standard_error_of_weighted_mean >,
				double >;

	// Pair data
#if(1)
	bin_stat_vec_t<BRG_DISTANCE> _R_values_;
	bin_stat_vec_t<BRG_MASS> _m_values_;
	bin_stat_vec_t<double> _z_values_;
	bin_stat_vec_t<double> _lens_z_values_;
	bin_stat_vec_t<double> _lens_unmasked_fracs_;
	bin_stat_vec_t<double> _source_z_values_;
	bin_stat_vec_t<double> _mag_lens_values_;
	stat_vec_t<BRG_UNITS> _delta_Sigma_t_values_;
	stat_vec_t<BRG_UNITS> _delta_Sigma_x_values_;

	stat_vec_t<double> _mu_obs_values_;
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
			double init_z_min=0, double init_z_max=0,
			double init_mag_min=0, double init_mag_max=0,
			double init_z_buffer=_z_buffer_default_value_);
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
	void add_lens( const size_t & lens_id, const double & lens_z, const double & unmasked_frac=1);
	void clear();

#endif

	// Count
	size_t count() const
	{
		return boost::accumulators::count(_R_values_);
	}
	double effective_count() const
	{
		return boost::accumulators::effective_count(_delta_Sigma_t_values_);
	}
	double sum_of_weights() const
	{
		return boost::accumulators::sum_of_weights(_delta_Sigma_t_values_);
	}
	double sum_of_square_weights() const
	{
		return boost::accumulators::sum_of_square_weights(_delta_Sigma_t_values_);
	}
	size_t num_lenses() const
	{
		return _distinct_lens_ids_.size();
	}

	// Limits and means accessors
#if(1)

	BRG_DISTANCE R_mean() const
	{
		return boost::accumulators::weighted_mean(_R_values_);
	}

	BRG_MASS m_mean() const
	{
		return boost::accumulators::weighted_mean(_m_values_);
	}

	double z_mean() const
	{
		return boost::accumulators::weighted_mean(_z_values_);
	}

	double lens_z_mean() const
	{
		return boost::accumulators::weighted_mean(_lens_z_values_);
	}

	double source_z_mean() const
	{
		return boost::accumulators::weighted_mean(_source_z_values_);
	}

	double mag_mean() const
	{
		return boost::accumulators::weighted_mean(_mag_lens_values_);
	}

	double unmasked_frac() const
	{
		double result = boost::accumulators::weighted_mean(_lens_unmasked_fracs_);
		if(isgood(result)) return result;
		return 0;
	}

	BRG_UNITS area() const
	{
		return unmasked_frac()*num_lenses()*pi*(square(afd(R_max(),lens_z_mean()))-square(afd(R_min(),lens_z_mean())));
	}

#endif

	// Accessors to pair values
#if (1)

	bin_stat_vec_t<BRG_DISTANCE> R_values()
	{
		return _R_values_;
	}
	bin_stat_vec_t<BRG_MASS> m_values()
	{
		return _m_values_;
	}
	bin_stat_vec_t<double> z_values()
	{
		return _z_values_;
	}
	bin_stat_vec_t<double> mag_lens_values()
	{
		return _mag_lens_values_;
	}
	stat_vec_t<BRG_UNITS> delta_Sigma_t_values()
	{
		return _delta_Sigma_t_values_;
	}
	stat_vec_t<BRG_UNITS> delta_Sigma_x_values()
	{
		return _delta_Sigma_x_values_;
	}
	stat_vec_t<double> mu_obs_values()
	{
		return _mu_obs_values_;
	}

#endif

	// Calculations on stored values
#if (1)

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
