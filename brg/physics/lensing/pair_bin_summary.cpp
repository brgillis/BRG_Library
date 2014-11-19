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
#include <fstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "brg/global.h"

#include "brg/file_access/open_file.hpp"
#include "brg/math/calculus/integrate.hpp"
#include "brg/math/misc_math.hpp"
#include "brg/physics/lensing/magnification/magnification_alpha.h"
#include "brg/physics/lensing/magnification/magnification_functors.h"
#include "brg/physics/lensing/pair_bin.h"
#include "brg/physics/units/unit_obj.h"
#include "brg/vector/summary_functions.hpp"

#include "pair_bin_summary.h"

namespace brgastro {

pair_bin_summary::pair_bin_summary( CONST_BRG_DISTANCE_REF init_R_min, CONST_BRG_DISTANCE_REF init_R_max,
		CONST_BRG_MASS_REF init_m_min, CONST_BRG_MASS_REF init_m_max,
		const double & init_z_min, const double & init_z_max,
		const double & init_mag_min, const double & init_mag_max )
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
	_source_z_mean_(0),
	_mag_mean_(0),
	_delta_Sigma_t_mean_(0),
	_delta_Sigma_x_mean_(0),
	_delta_Sigma_t_mean_square_(0),
	_delta_Sigma_x_mean_square_(0),
	_mu_hat_(0),
	_mu_W_(0),
	_total_area_(0),
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
	_source_z_mean_(bin.source_z_mean()),
	_mag_mean_(bin.mag_mean()),
	_delta_Sigma_t_mean_(bin.delta_Sigma_t_mean()),
	_delta_Sigma_x_mean_(bin.delta_Sigma_x_mean()),
	_delta_Sigma_t_mean_square_(bin.delta_Sigma_t_mean_square()),
	_delta_Sigma_x_mean_square_(bin.delta_Sigma_x_mean_square()),
	_mu_hat_(bin.mu_hat()),
	_mu_W_(bin.mu_W()),
	_total_area_(bin.area()),
	_num_lenses_(bin.num_lenses())
{
}

pair_bin_summary::pair_bin_summary( std::istream & in )
{
	load(in);
}
pair_bin_summary::pair_bin_summary( std::string & filename )
{
	load(filename);
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
	_source_z_mean_ = 0;
	_mag_mean_ = 0;
	_delta_Sigma_t_mean_ = 0;
	_delta_Sigma_x_mean_ = 0;
	_delta_Sigma_t_mean_square_ = 0;
	_delta_Sigma_x_mean_square_ = 0;
	_num_lenses_ = 0;
	_uncache_values();
}

#endif // Adding and clearing data

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

	// Check for zero area in this or the other
	if(other.area()==0) return *this;

	_uncache_values(); // If we get this far, something is changing with this bin

	if(area()==0)
	{
		// Simply copy the other bin summary to this one
		*this = other;
		return *this;
	}

	// Add the count and magnification data
	_count_ += other.count();
	_num_lenses_ += other.num_lenses();

	auto new_W = _mu_W_ + other.mu_W();
	_mu_hat_ = (_mu_hat_*_mu_W_ + other.mu_hat()*other.mu_W())/new_W;
	_mu_W_ = new_W;

	// Check for zero total weight in this or the other
	if(other.sum_of_weights()==0)
	{
		return *this;
	}
	if(_sum_of_weights_==0)
	{
		// Simply copy the other bin summary, except for the magnification data
		auto tmp_count = std::move(_count_);
		auto tmp_num_lenses = std::move(_num_lenses_);
		auto tmp_mu_hat = std::move(_mu_hat_);
		auto tmp_mu_W = std::move(_mu_W_);
		*this = other;
		_count_ = std::move(tmp_count);
		_num_lenses_ = std::move(tmp_num_lenses);
		_mu_hat_ = std::move(tmp_mu_hat);
		_mu_W_ = std::move(tmp_mu_W);
		return *this;
	}

	auto new_sum_of_weights = sum_of_weights() + other.sum_of_weights();

	_R_mean_ = ( _R_mean_*sum_of_weights() + other.R_mean()*other.sum_of_weights())/new_sum_of_weights;
	_m_mean_ = ( _m_mean_*sum_of_weights() + other.m_mean()*other.sum_of_weights())/new_sum_of_weights;
	_z_mean_ = ( _z_mean_*sum_of_weights() + other.z_mean()*other.sum_of_weights())/new_sum_of_weights;
	_source_z_mean_ = ( _source_z_mean_*sum_of_weights() + other.source_z_mean()*other.sum_of_weights())/new_sum_of_weights;
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

// Fix bad values
void pair_bin_summary::fixbad()
{
	brgastro::fixbad(_R_min_);
	brgastro::fixbad(_R_max_);
	brgastro::fixbad(_m_min_);
	brgastro::fixbad(_m_max_);
	brgastro::fixbad(_z_min_);
	brgastro::fixbad(_z_max_);
	brgastro::fixbad(_mag_min_);
	brgastro::fixbad(_mag_max_);
	brgastro::fixbad(_count_);
	brgastro::fixbad(_sum_of_weights_);
	brgastro::fixbad(_sum_of_square_weights_);
	brgastro::fixbad(_R_mean_);
	brgastro::fixbad(_m_mean_);
	brgastro::fixbad(_z_mean_);
	brgastro::fixbad(_source_z_mean_);
	brgastro::fixbad(_mag_mean_);
	brgastro::fixbad(_delta_Sigma_t_mean_);
	brgastro::fixbad(_delta_Sigma_x_mean_);
	brgastro::fixbad(_delta_Sigma_t_mean_square_);
	brgastro::fixbad(_delta_Sigma_x_mean_square_);
	brgastro::fixbad(_mu_hat_);
	brgastro::fixbad(_mu_W_);
	brgastro::fixbad(_total_area_);
	brgastro::fixbad(_num_lenses_);
}

// Saving/loading data
#if(1)

void pair_bin_summary::save(std::ostream & out) const
{
	boost::archive::text_oarchive ar(out);
	ar << *this;
}
void pair_bin_summary::save(const std::string & filename) const
{
	std::ofstream out;
	open_file_output(out,filename);
	save(out);
}
void pair_bin_summary::save(std::ostream & out)
{
	fixbad();
	boost::archive::text_oarchive ar(out);
	ar << *this;
}
void pair_bin_summary::save(const std::string & filename)
{
	std::ofstream out;
	open_file_output(out,filename);
	save(out);
}
void pair_bin_summary::load(std::istream & in)
{
	boost::archive::text_iarchive ar(in);
	ar >> *this;
}
void pair_bin_summary::load(const std::string & filename)
{
	std::ifstream in;
	open_file_input(in,filename);
	load(in);
}

#endif



// Summary values
#if (1)

BRG_UNITS pair_bin_summary::sigma_crit() const
{
	#ifdef _BRG_USE_UNITS_
	BRG_UNITS result(brgastro::sigma_crit(z_mean(),source_z_mean()),-2,0,1,0,0,0);
	#else
	BRG_UNITS result = brgastro::sigma_crit(z_mean(),source_z_mean());
	#endif
	return result;
}

// Shear
#if (1)

BRG_UNITS pair_bin_summary::delta_Sigma_t_std() const
{
	return std::sqrt( _delta_Sigma_t_mean_square_ - square(_delta_Sigma_t_mean_) );
}
BRG_UNITS pair_bin_summary::delta_Sigma_x_std() const
{
	return std::sqrt( _delta_Sigma_x_mean_square_ - square(_delta_Sigma_x_mean_) );
}

BRG_UNITS pair_bin_summary::delta_Sigma_t_stderr() const
{
	if(_count_<2) return std::numeric_limits<double>::max();
	BRG_UNITS result = delta_Sigma_t_std()*std::sqrt(_count_/(effective_count()*(_count_-1)) );
	if(isgood(result)) return result;
	return 0;
}
BRG_UNITS pair_bin_summary::delta_Sigma_x_stderr() const
{
	if(_count_<2) return std::numeric_limits<double>::max();
	BRG_UNITS result = delta_Sigma_x_std()*std::sqrt(_count_/(effective_count()*(_count_-1)) );
	if(isgood(result)) return result;
	return 0;
}

double pair_bin_summary::gamma_t_mean() const
{
	return delta_Sigma_t_mean()/sigma_crit();
}
double pair_bin_summary::gamma_x_mean() const
{
	return delta_Sigma_x_mean()/sigma_crit();
}
double pair_bin_summary::gamma_mean() const
{
	return std::sqrt(gamma_mean_square());
}
double pair_bin_summary::gamma_mean_square() const
{
	return square(gamma_t_mean())+square(gamma_x_mean());

}

double pair_bin_summary::gamma_t_stderr() const
{
	return delta_Sigma_t_stderr()/sigma_crit();
}
double pair_bin_summary::gamma_x_stderr() const
{
	return delta_Sigma_x_stderr()/sigma_crit();
}
double pair_bin_summary::gamma_stderr() const
{
	return quad_add_err(gamma_t_mean(),gamma_t_stderr(),gamma_x_mean(),gamma_x_stderr());
}
double pair_bin_summary::gamma_square_stderr() const
{
	return 2*gamma_mean()*gamma_stderr();
}

#endif // Shear

// Magnification
#if (1)

BRG_UNITS pair_bin_summary::area_per_lens() const
{
	return area()/num_lenses();
}
double pair_bin_summary::mu_stderr() const
{
	// TODO Correct this for weighted lenses and pairs used for magnification if needed
	return 1./std::sqrt(mu_W());
}
double pair_bin_summary::kappa() const
{
	double gms = gamma_mean_square();
	brgastro::fixbad(gms);
	return 1.-std::sqrt(1/mu_hat()+gms);
}
double pair_bin_summary::kappa_stderr() const
{
	double gms = gamma_mean_square();
	double gserr = gamma_square_stderr();
	brgastro::fixbad(gms);
	brgastro::fixbad(gserr);
	return sqrt_err(1/mu_hat()+gamma_mean_square(),quad_add(mu_stderr(),gamma_square_stderr()));
}

#endif // Magnification

#endif

} // end namespace brgastro
