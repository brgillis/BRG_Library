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

#include <set>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/error_of_mean.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/stats.hpp>

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
				boost::accumulators::tag::mean > >;
	template <typename T>
	using stat_vec_t = boost::accumulators::accumulator_set<T,
			boost::accumulators::stats<
				boost::accumulators::tag::count,
				boost::accumulators::tag::mean,
				boost::accumulators::tag::moment<2>,
				boost::accumulators::tag::variance,
				boost::accumulators::tag::error_of<boost::accumulators::tag::mean> > >;

	// Pair data
#if(1)
	bin_stat_vec_t<BRG_DISTANCE> _R_values_;
	bin_stat_vec_t<BRG_MASS> _m_values_;
	bin_stat_vec_t<double> _z_values_;
	bin_stat_vec_t<double> _mag_lens_values_;
	stat_vec_t<double> _mag_source_values_;
	stat_vec_t<BRG_UNITS> _delta_Sigma_t_values_;
	stat_vec_t<BRG_UNITS> _delta_Sigma_x_values_;

	std::set<size_t> _distinct_lens_ids_;

#endif // Pair data

public:

	// Constructors and destructor
#if(1)
	pair_bin( CONST_BRG_DISTANCE_REF init_R_min=0, CONST_BRG_DISTANCE_REF init_R_max=0,
			CONST_BRG_MASS_REF init_m_min=0, CONST_BRG_MASS_REF init_m_max=0,
			double init_z_min=0, double init_z_max=0,
			double init_mag_min=0, double init_mag_max=0 );
	virtual ~pair_bin()
	{
	}
#endif

	// Adding and clearing data
#if(1)

	void add_pair( const lens_source_pair & new_pair);
	void clear();

#endif

	// Count
	size_t count() const
	{
		return boost::accumulators::count(_R_values_);
	}
	size_t num_lenses() const
	{
		return _distinct_lens_ids_.size();
	}

	// Limits and means accessors
#if(1)

	BRG_DISTANCE R_mean() const
	{
		return boost::accumulators::mean(_R_values_);
	}

	BRG_MASS m_mean() const
	{
		return boost::accumulators::mean(_m_values_);
	}

	double z_mean() const
	{
		return boost::accumulators::mean(_z_values_);
	}

	double mag_mean() const
	{
		return boost::accumulators::mean(_mag_lens_values_);
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
	stat_vec_t<double> mag_source_values()
	{
		return _mag_source_values_;
	}
	stat_vec_t<BRG_UNITS> delta_Sigma_t_values()
	{
		return _delta_Sigma_t_values_;
	}
	stat_vec_t<BRG_UNITS> delta_Sigma_x_values()
	{
		return _delta_Sigma_x_values_;
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

#endif
};

} // end namespace brgastro

#endif // _BRG_PAIR_BIN_H_INCLUDED_
