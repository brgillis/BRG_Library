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

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "brg/global.h"

#include "brg/file_access/table_typedefs.hpp"
#include "brg/physics/units/unit_obj.h"
#include "brg/vector/summary_functions.hpp"

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
	double _source_z_mean_;
	double _mag_mean_;

	BRG_UNITS _delta_Sigma_t_mean_, _delta_Sigma_x_mean_;
	BRG_UNITS _delta_Sigma_t_mean_square_, _delta_Sigma_x_mean_square_;

	double _mu_hat_;
	double _mu_W_;
	BRG_UNITS _total_area_;
	size_t _num_lenses_;

#endif

	// Serialization (to save it to a file)
#if(1)
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & _R_min_ & _R_max_;
		ar & _m_min_ & _m_max_;
		ar & _z_min_ & _z_max_;
		ar & _mag_min_ & _mag_max_;

		ar & _count_;
		ar & _sum_of_weights_;
		ar & _sum_of_square_weights_;
		ar & _R_mean_;
		ar & _m_mean_;
		ar & _z_mean_;
		ar & _source_z_mean_;
		ar & _mag_mean_;

		ar & _delta_Sigma_t_mean_ & _delta_Sigma_x_mean_;
		ar & _delta_Sigma_t_mean_square_ & _delta_Sigma_x_mean_square_;

		ar & _mu_hat_;
		ar & _mu_W_;
		ar & _total_area_;
		ar & _num_lenses_;
	}
#endif

protected:

	virtual void _uncache_values()
	{
	}

public:

	// Constructors and destructor
#if(1)
	pair_bin_summary( CONST_BRG_DISTANCE_REF init_R_min=0, CONST_BRG_DISTANCE_REF init_R_max=0,
			CONST_BRG_MASS_REF init_m_min=0, CONST_BRG_MASS_REF init_m_max=0,
			double init_z_min=0, double init_z_max=0,
			double init_mag_min=0, double init_mag_max=0 );
	pair_bin_summary( const pair_bin & bin);
	pair_bin_summary( std::istream & in );
	pair_bin_summary( std::string & filename );
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

	// Fix bad values
	void fixbad();

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
	virtual double lens_z_mean() const
	{
		return _z_mean_;
	}
	virtual double source_z_mean() const
	{
		return _source_z_mean_;
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

	double sigma_crit() const
	{
		return brgastro::sigma_crit(z_mean(),source_z_mean());
	}

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

	double gamma_t_mean() const
	{
		return delta_Sigma_t_mean()/sigma_crit();
	}
	double gamma_x_mean() const
	{
		return delta_Sigma_x_mean()/sigma_crit();
	}
	double gamma_mean() const
	{
		return std::sqrt(gamma_mean_square());
	}
	double gamma_mean_square() const
	{
		return square(gamma_t_mean())+square(gamma_x_mean());
	}

	double gamma_t_stderr() const
	{
		return delta_Sigma_t_stderr()/sigma_crit();
	}
	double gamma_x_stderr() const
	{
		return delta_Sigma_x_stderr()/sigma_crit();
	}
	double gamma_stderr() const
	{
		return quad_add_err(gamma_t_mean(),gamma_t_stderr(),gamma_x_mean(),gamma_x_stderr());
	}
	double gamma_square_stderr() const
	{
		return 2*gamma_mean()*gamma_stderr();
	}

#endif // Shear

	// Magnification
#if (1)

	BRG_UNITS area_per_lens() const
	{
		return area()/num_lenses();
		return pi*(square(afd(_R_max_,_z_mean_))-square(afd(_R_max_,_z_mean_)));
	}
	virtual BRG_UNITS area() const
	{
		return _total_area_;
	}
	virtual size_t num_lenses() const
	{
		return _num_lenses_;
	}
	virtual double mu_hat() const
	{
		return _mu_hat_;
	}
	virtual double mu_W() const
	{
		return _mu_W_;
	}
	double mu_stderr() const
	{
		// TODO Correct this for weighted lenses and pairs used for magnification if needed
		return 1./std::sqrt(mu_W());
	}
	double kappa() const
	{
		return 1.-std::sqrt(1/mu_hat()+gamma_mean_square());
	}
	double kappa_stderr() const
	{
		return sqrt_err(1/mu_hat()+gamma_mean_square(),quad_add(mu_stderr(),gamma_square_stderr()));
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

	// Saving/loading data
#if(1)

	void save(std::ostream &) const;
	void save(const std::string &) const;
	void save(std::ostream &);
	void save(const std::string &);
	void load(std::istream &);
	void load(const std::string &);

#endif
};

} // end namespace brgastro

#endif // _BRG_PAIR_BIN_H_INCLUDED_
