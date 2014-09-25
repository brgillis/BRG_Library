/**********************************************************************\
 @file pair_bin_summary.h
 ------------------

 A class representing a bin of lens-source pairs, which can be used for
 calculating statistics on the contained pairs. This base class only
 includes the summary data.

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
#ifndef _BRG_PAIR_BIN_SUMMARY_H_INCLUDED_
#define _BRG_PAIR_BIN_SUMMARY_H_INCLUDED_

#include <utility>
#include <vector>

#include "brg/global.h"

#include "brg/physics/units/unit_obj.h"
#include "brg/vector/summary_functions.hpp"

// TODO: Fix effective count use in summation (should add sum of weights and sum of square weights)

namespace brgastro {

// Forward-declare pair_bin
class pair_bin;

/**
 *
 */
class pair_bin_summary {
private:

	// Bin limits
#if(1)

	BRG_DISTANCE _R_min_, _R_max_;
	BRG_MASS _m_min_, _m_max_;
	double _z_min_, _z_max_;
	double _mag_min_, _mag_max_;

#endif // Bin limits

	// Summary data
#if (1)

	size_t _count_;
	double _sum_of_weights_;
	double _sum_of_square_weights_;
	BRG_DISTANCE _R_mean_;
	BRG_MASS _m_mean_;
	double _z_mean_;
	double _mag_mean_;

	BRG_UNITS _delta_Sigma_t_mean_, _delta_Sigma_x_mean_;
	BRG_UNITS _delta_Sigma_t_mean_square_, _delta_Sigma_x_mean_square_;

	mutable double _mu_hat_cached_value_;
	mutable double _mu_W_cached_value_;
	std::vector<double> _source_magnitudes_;
	size_t _num_lenses_;

#endif

protected:

	void _uncache_values()
	{
		_mu_hat_cached_value_ = std::numeric_limits<double>::infinity();
		_mu_W_cached_value_ = std::numeric_limits<double>::infinity();
	}

public:

	// Constructors and destructor
#if(1)
	pair_bin_summary( CONST_BRG_DISTANCE_REF init_R_min=0, CONST_BRG_DISTANCE_REF init_R_max=0,
			CONST_BRG_MASS_REF init_m_min=0, CONST_BRG_MASS_REF init_m_max=0,
			double init_z_min=0, double init_z_max=0,
			double init_mag_min=0, double init_mag_max=0 );
	pair_bin_summary( const pair_bin & bin);
	pair_bin_summary & operator=(const pair_bin & bin)
	{
		*this = pair_bin_summary(bin);
		return *this;
	}
	virtual ~pair_bin_summary()
	{
	}
#endif

	// Setting bin limits
#if(1)

	template <typename T1, typename T2>
	void set_R_limits(T1&& R_min, T2&& R_max)
	{
		_R_min_ = std::forward<T1>(R_min);
		_R_max_ = std::forward<T2>(R_max);
	}
	template <typename T1, typename T2>
	void set_m_limits(T1&& m_min, T2&& m_max)
	{
		_m_min_ = std::forward<T1>(m_min);
		_m_max_ = std::forward<T2>(m_max);
	}
	template <typename T1, typename T2>
	void set_z_limits(T1&& z_min, T2&& z_max)
	{
		_z_min_ = std::forward<T1>(z_min);
		_z_max_ = std::forward<T2>(z_max);
	}
	template <typename T1, typename T2>
	void set_mag_limits(T1&& mag_min, T2&& mag_max)
	{
		_mag_min_ = std::forward<T1>(mag_min);
		_mag_max_ = std::forward<T2>(mag_max);
	}

	template <typename T1, typename T2,
				typename T3, typename T4,
				typename T5, typename T6,
				typename T7, typename T8>
	void set_limits(T1&& R_min, T2&& R_max,
			T3&& m_min, T4&& m_max,
			T5&& z_min, T6&& z_max,
			T7&& mag_min, T8&& mag_max)
	{
		set_R_limits(std::forward<T1>(R_min),std::forward<T2>(R_max));
		set_m_limits(std::forward<T3>(m_min),std::forward<T4>(m_max));
		set_z_limits(std::forward<T5>(z_min),std::forward<T6>(z_max));
		set_mag_limits(std::forward<T7>(mag_min),std::forward<T8>(mag_max));
	}


#endif

	// Adding and clearing data
#if(1)

	virtual void clear();

#endif

	// General statistics accessors
#if(1)

	// Count
	virtual size_t count() const
	{
		return _count_;
	}
	// Effective count
	virtual double effective_count() const
	{
		return square(sum_of_weights())/sum_of_square_weights();
	}
	// Sum of weights
	virtual double sum_of_weights() const
	{
		return _sum_of_weights_;
	}
	// Sum of square weights
	virtual double sum_of_square_weights() const
	{
		return _sum_of_square_weights_;
	}

#endif

	// Limits and means accessors
#if(1)

	CONST_BRG_DISTANCE_REF R_min() const
	{
		return _R_min_;
	}
	CONST_BRG_DISTANCE_REF R_max() const
	{
		return _R_max_;
	}
	virtual BRG_DISTANCE R_mean() const
	{
		return _R_mean_;
	}

	CONST_BRG_MASS_REF m_min() const
	{
		return _m_min_;
	}
	CONST_BRG_MASS_REF m_max() const
	{
		return _m_max_;
	}
	virtual BRG_MASS m_mean() const
	{
		return _m_mean_;
	}

	double z_min() const
	{
		return _z_min_;
	}
	double z_max() const
	{
		return _z_max_;
	}
	virtual double z_mean() const
	{
		return _z_mean_;
	}

	double mag_min() const
	{
		return _mag_min_;
	}
	double mag_max() const
	{
		return _mag_max_;
	}
	virtual double mag_mean() const
	{
		return _mag_mean_;
	}

#endif

	// Summary values
#if (1)

	// Shear
#if (1)
	virtual BRG_UNITS delta_Sigma_t_mean() const
	{
		return _delta_Sigma_t_mean_;
	}
	virtual BRG_UNITS delta_Sigma_x_mean() const
	{
		return _delta_Sigma_x_mean_;
	}

	virtual BRG_UNITS delta_Sigma_t_mean_square() const
	{
		return _delta_Sigma_t_mean_square_;
	}
	virtual BRG_UNITS delta_Sigma_x_mean_square() const
	{
		return _delta_Sigma_x_mean_square_;
	}

	virtual BRG_UNITS delta_Sigma_t_std() const
	{
		return std::sqrt( _delta_Sigma_t_mean_square_ - square(_delta_Sigma_t_mean_) );
	}
	virtual BRG_UNITS delta_Sigma_x_std() const
	{
		return std::sqrt( _delta_Sigma_x_mean_square_ - square(_delta_Sigma_x_mean_) );
	}

	virtual BRG_UNITS delta_Sigma_t_stderr() const
	{
		if(_count_<2) return std::numeric_limits<double>::max();
		return delta_Sigma_t_std()*std::sqrt(_count_/(effective_count()*(_count_-1)) );
	}
	virtual BRG_UNITS delta_Sigma_x_stderr() const
	{
		if(_count_<2) return std::numeric_limits<double>::max();
		return delta_Sigma_x_std()*std::sqrt(_count_/(effective_count()*(_count_-1)) );
	}

#endif // Shear

	// Magnification
#if (1)

	BRG_UNITS area_per_lens() const
	{
		return pi*(square(_R_max_)-square(_R_min_));
	}
	BRG_UNITS area() const
	{
		return area_per_lens()*num_lenses();
	}
	virtual size_t num_lenses() const
	{
		return _num_lenses_;
	}
	double mu_hat() const;
	double mu_W() const;
	const std::vector<double> & source_magnitudes() const
	{
		return _source_magnitudes_;
	}

#endif // Magnification

#endif

	// Combining summaries together
#if(1)

	pair_bin_summary & operator+=( const pair_bin_summary & other );
	pair_bin_summary operator+( pair_bin_summary other ) const
	{
		return other += *this;
	}

#endif // Combining summaries together
};

} // end namespace brgastro

#endif // _BRG_PAIR_BIN_H_INCLUDED_
