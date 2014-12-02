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

// Private methods
#if(1)

double pair_bin_summary::_magf_gamma_t_mean() const
{
	return delta_Sigma_t_mean()/magf_sigma_crit();
}
double pair_bin_summary::_magf_gamma_x_mean() const
{
	return delta_Sigma_x_mean()/magf_sigma_crit();
}
double pair_bin_summary::_magf_gamma_mean() const
{
	return std::sqrt(gamma_mean_square());
}
double pair_bin_summary::_magf_gamma_mean_square() const
{
	return square(gamma_t_mean())+square(gamma_x_mean());

}

double pair_bin_summary::_magf_gamma_t_stderr() const
{
	return delta_Sigma_t_stderr()/magf_sigma_crit();
}
double pair_bin_summary::_magf_gamma_x_stderr() const
{
	return delta_Sigma_x_stderr()/magf_sigma_crit();
}
double pair_bin_summary::_magf_gamma_stderr() const
{
	return quad_add_err(gamma_t_mean(),gamma_t_stderr(),gamma_x_mean(),gamma_x_stderr());
}
double pair_bin_summary::_magf_gamma_square_stderr() const
{
	return 2*gamma_mean()*gamma_stderr();
}

lensing_tNFW_profile pair_bin_summary::_shear_model_profile(const double & MLratio) const
{
	return lensing_tNFW_profile(MLratio*shear_m_mean(),shear_z_mean());
}

lensing_tNFW_profile pair_bin_summary::_magf_model_profile(const double & MLratio) const
{
	return lensing_tNFW_profile(MLratio*magf_m_mean(),magf_z_mean());
}

#endif

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

	_shear_pair_count_(0),
	_magf_pair_count_(0),

	_shear_R_mean_(0),
	_magf_R_mean_(0),

	_shear_lens_m_mean_(0),
	_magf_lens_m_mean_(0),
	_shear_lens_z_mean_(0),
	_magf_lens_z_mean_(0),
	_shear_lens_mag_mean_(0),
	_magf_lens_mag_mean_(0),

	_shear_source_z_mean_(0),
	_magf_source_z_mean_(0),

	_shear_sum_of_weights_(0),
	_shear_sum_of_square_weights_(0),
	_delta_Sigma_t_mean_(0),
	_delta_Sigma_x_mean_(0),
	_delta_Sigma_t_mean_square_(0),
	_delta_Sigma_x_mean_square_(0),

	_mu_hat_(0),
	_mu_W_(0),
	_area_(0),
	_magf_lens_count_(0)
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

	_shear_pair_count_(bin.shear_count()),
	_magf_pair_count_(bin.magf_count()),

	_shear_R_mean_(bin.shear_R_mean()),
	_magf_R_mean_(bin.magf_R_mean()),

	_shear_lens_m_mean_(bin.shear_lens_m_mean()),
	_magf_lens_m_mean_(bin.magf_lens_m_mean()),
	_shear_lens_z_mean_(bin.shear_lens_z_mean()),
	_magf_lens_z_mean_(bin.magf_lens_z_mean()),
	_shear_lens_mag_mean_(bin.shear_lens_mag_mean()),
	_magf_lens_mag_mean_(bin.magf_lens_mag_mean()),

	_shear_source_z_mean_(bin.shear_source_z_mean()),
	_magf_source_z_mean_(bin.magf_source_z_mean()),

	_shear_sum_of_weights_(bin.shear_sum_of_weights()),
	_shear_sum_of_square_weights_(bin.shear_sum_of_square_weights()),
	_delta_Sigma_t_mean_(bin.delta_Sigma_t_mean()),
	_delta_Sigma_x_mean_(bin.delta_Sigma_x_mean()),
	_delta_Sigma_t_mean_square_(bin.delta_Sigma_t_mean_square()),
	_delta_Sigma_x_mean_square_(bin.delta_Sigma_x_mean_square()),

	_mu_hat_(bin.mu_hat()),
	_mu_W_(bin.mu_W()),
	_area_(bin.area()),
	_magf_lens_count_(bin.magf_num_lenses())
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
	_shear_pair_count_ = 0;
	_magf_pair_count_ = 0;

	_shear_R_mean_ = 0;
	_magf_R_mean_ = 0;

	_shear_lens_m_mean_ = 0;
	_magf_lens_m_mean_ = 0;
	_shear_lens_z_mean_ = 0;
	_magf_lens_z_mean_ = 0;
	_shear_lens_mag_mean_ = 0;
	_magf_lens_mag_mean_ = 0;

	_shear_source_z_mean_ = 0;
	_magf_source_z_mean_ = 0;

	_shear_sum_of_weights_ = 0;
	_shear_sum_of_square_weights_ = 0;

	_delta_Sigma_t_mean_ = 0;
	_delta_Sigma_x_mean_ = 0;
	_delta_Sigma_t_mean_square_ = 0;
	_delta_Sigma_x_mean_square_ = 0;
	_magf_lens_count_ = 0;
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

	// Check whether this or the other is empty in shear or magnification data
	const bool this_magf_zero = (area()==0);
	const bool other_magf_zero = (other.area()==0);
	const bool this_shear_zero = (shear_sum_of_weights()==0);
	const bool other_shear_zero = (other.shear_sum_of_weights()==0);

	// Check if the other is completely empty
	if(other_magf_zero && other_shear_zero) return *this;

	// Check if this is completely empty
	if(this_magf_zero && this_shear_zero)
	{
		// Simply copy the other bin summary to this one
		*this = other;
		return *this;
	}

	_uncache_values(); // If we get this far, something is changing with this bin

	// If the other's magnification data isn't empty
	if(!other_magf_zero)
	{

		// If this is empty, simply copy everything
		if(this_magf_zero)
		{
			// Add the magnification data
			_magf_pair_count_ = other.magf_pair_count();
			_magf_lens_count_ = other.magf_num_lenses();
			_area_ = other.area();

			_magf_R_mean_ = other.magf_R_mean();

			_magf_lens_m_mean_ = other.magf_m_mean();
			_magf_lens_z_mean_ = other.magf_z_mean();
			_magf_lens_mag_mean_ = other.magf_mag_mean();

			_magf_source_z_mean_ = other.magf_source_z_mean();

			_mu_hat_ = other.mu_hat();

			_mu_W_ = other.mu_W();
		}
		else
		{
			// Add the magnification data
			_magf_pair_count_ += other.magf_pair_count();
			_magf_lens_count_ += other.magf_num_lenses();
			_area_ += other.area();

			brgastro::fixbad(_mu_W_);
			const auto o_mu_W = other.mu_W();

			const auto new_mu_W = _mu_W_ + o_mu_W;

			_magf_R_mean_ = (magf_R_mean()*_mu_W_ + other.magf_R_mean()*o_mu_W)/new_mu_W;

			_magf_lens_m_mean_ = (magf_m_mean()*_mu_W_ + other.magf_m_mean()*o_mu_W)/new_mu_W;
			_magf_lens_z_mean_ = (magf_z_mean()*_mu_W_ + other.magf_z_mean()*o_mu_W)/new_mu_W;
			_magf_lens_mag_mean_ = (magf_mag_mean()*_mu_W_ + other.magf_mag_mean()*o_mu_W)/new_mu_W;

			_magf_source_z_mean_ = (magf_source_z_mean()*_mu_W_ + other.magf_source_z_mean()*o_mu_W)/new_mu_W;

			_mu_hat_ = (mu_hat()*_mu_W_ + other.mu_hat()*o_mu_W)/new_mu_W;

			_mu_W_ = new_mu_W;
		}
	}

	// Add the shear data

	// If the other's shear data isn't empty
	if(!other_shear_zero)
	{

		// If this is empty, simply copy everything
		if(this_shear_zero)
		{
			_shear_pair_count_ = other.shear_pair_count();

			_shear_R_mean_ = other.shear_R_mean();
			_shear_lens_m_mean_ = other.shear_m_mean();
			_shear_lens_z_mean_ = other.shear_z_mean();
			_shear_lens_mag_mean_ = other.mag_mean();
			_shear_source_z_mean_ = other.source_z_mean();

			_delta_Sigma_t_mean_ = other.delta_Sigma_t_mean();
			_delta_Sigma_t_mean_square_ = other.delta_Sigma_t_mean_square();
			_delta_Sigma_x_mean_ = other.delta_Sigma_x_mean();
			_delta_Sigma_x_mean_square_ = other.delta_Sigma_x_mean_square();
			_shear_sum_of_weights_ = other.shear_sum_of_weights();
			_shear_sum_of_square_weights_ = other.shear_sum_of_square_weights();
		}
		else
		{
			_shear_pair_count_ += other.shear_pair_count();

			const auto o_shear_sum_of_weights = other.shear_sum_of_weights();

			const auto new_sum_of_weights = _shear_sum_of_weights_ + o_shear_sum_of_weights;

			_shear_R_mean_ = ( shear_R_mean()*_shear_sum_of_weights_ + other.shear_R_mean()*o_shear_sum_of_weights)/new_sum_of_weights;
			_shear_lens_m_mean_ = ( shear_m_mean()*_shear_sum_of_weights_ + other.shear_m_mean()*o_shear_sum_of_weights)/new_sum_of_weights;
			_shear_lens_z_mean_ = ( shear_z_mean()*_shear_sum_of_weights_ + other.shear_z_mean()*o_shear_sum_of_weights)/new_sum_of_weights;
			_shear_lens_mag_mean_ = ( shear_mag_mean()*_shear_sum_of_weights_ + other.mag_mean()*o_shear_sum_of_weights)/new_sum_of_weights;
			_shear_source_z_mean_ = ( shear_source_z_mean()*_shear_sum_of_weights_ + other.source_z_mean()*o_shear_sum_of_weights)/new_sum_of_weights;

			_delta_Sigma_t_mean_ = ( delta_Sigma_t_mean()*_shear_sum_of_weights_ + other.delta_Sigma_t_mean()*o_shear_sum_of_weights)
					/new_sum_of_weights;
			_delta_Sigma_t_mean_square_ = ( delta_Sigma_t_mean_square()*shear_sum_of_weights() +
					other.delta_Sigma_t_mean_square()*o_shear_sum_of_weights)
					/new_sum_of_weights;
			_delta_Sigma_x_mean_ = ( delta_Sigma_x_mean()*shear_sum_of_weights() + other.delta_Sigma_x_mean()*o_shear_sum_of_weights)
							/new_sum_of_weights;
			_delta_Sigma_x_mean_square_ = ( delta_Sigma_x_mean_square()*shear_sum_of_weights() +
					other.delta_Sigma_x_mean_square()*o_shear_sum_of_weights)
					/new_sum_of_weights;
			_shear_sum_of_weights_ = new_sum_of_weights;
			_shear_sum_of_square_weights_ += other.shear_sum_of_square_weights();
		}
	}

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

	brgastro::fixbad(_shear_R_mean_);
	brgastro::fixbad(_magf_R_mean_);

	brgastro::fixbad(_shear_lens_m_mean_);
	brgastro::fixbad(_magf_lens_m_mean_);
	brgastro::fixbad(_shear_lens_z_mean_);
	brgastro::fixbad(_magf_lens_z_mean_);
	brgastro::fixbad(_shear_lens_mag_mean_);
	brgastro::fixbad(_magf_lens_mag_mean_);

	brgastro::fixbad(_shear_source_z_mean_);
	brgastro::fixbad(_magf_source_z_mean_);

	brgastro::fixbad(_shear_sum_of_weights_);
	brgastro::fixbad(_shear_sum_of_square_weights_);

	brgastro::fixbad(_delta_Sigma_t_mean_);
	brgastro::fixbad(_delta_Sigma_x_mean_);
	brgastro::fixbad(_delta_Sigma_t_mean_square_);
	brgastro::fixbad(_delta_Sigma_x_mean_square_);

	brgastro::fixbad(_mu_hat_);
	brgastro::fixbad(_mu_W_);
	brgastro::fixbad(_area_);
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

BRG_UNITS pair_bin_summary::shear_sigma_crit() const
{
	#ifdef _BRG_USE_UNITS_
	BRG_UNITS result(brgastro::sigma_crit(shear_z_mean(),shear_source_z_mean()),-2,0,1,0,0,0);
	#else
	BRG_UNITS result = brgastro::sigma_crit(shear_z_mean(),shear_source_z_mean());
	#endif
	return result;
}

BRG_UNITS pair_bin_summary::magf_sigma_crit() const
{
	#ifdef _BRG_USE_UNITS_
	BRG_UNITS result(brgastro::sigma_crit(magf_z_mean(),magf_source_z_mean()),-2,0,1,0,0,0);
	#else
	BRG_UNITS result = brgastro::sigma_crit(magf_z_mean(),magf_source_z_mean());
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
	if(_shear_pair_count_<2) return std::numeric_limits<double>::max();
	BRG_UNITS result = delta_Sigma_t_std()*std::sqrt(_shear_pair_count_/(shear_effective_pair_count()*(_shear_pair_count_-1)) );
	if(isgood(result)) return result;
	return 0;
}
BRG_UNITS pair_bin_summary::delta_Sigma_x_stderr() const
{
	if(_shear_pair_count_<2) return std::numeric_limits<double>::max();
	BRG_UNITS result = delta_Sigma_x_std()*std::sqrt(_shear_pair_count_/(shear_effective_pair_count()*(_shear_pair_count_-1)) );
	if(isgood(result)) return result;
	return 0;
}

double pair_bin_summary::gamma_t_mean() const
{
	return delta_Sigma_t_mean()/shear_sigma_crit();
}
double pair_bin_summary::gamma_x_mean() const
{
	return delta_Sigma_x_mean()/shear_sigma_crit();
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
	return delta_Sigma_t_stderr()/shear_sigma_crit();
}
double pair_bin_summary::gamma_x_stderr() const
{
	return delta_Sigma_x_stderr()/shear_sigma_crit();
}
double pair_bin_summary::gamma_stderr() const
{
	return quad_add_err(gamma_t_mean(),gamma_t_stderr(),gamma_x_mean(),gamma_x_stderr());
}
double pair_bin_summary::gamma_square_stderr() const
{
	return 2*gamma_mean()*gamma_stderr();
}

BRG_UNITS pair_bin_summary::model_delta_Sigma_t(const double & MLratio) const
{
	BRG_UNITS result = _shear_model_profile(MLratio).WLsig(shear_R_mean());
	return result;
}
double pair_bin_summary::model_gamma_t(const double & MLratio) const
{
	double result = model_delta_Sigma_t(MLratio)/shear_sigma_crit();
	return result;
}

#endif // Shear

// Magnification
#if (1)

BRG_UNITS pair_bin_summary::area_per_lens() const
{
	return area()/magf_num_lenses();
}
double pair_bin_summary::mu_stderr() const
{
	// TODO Correct this for weighted lenses and pairs used for magnification if needed
	return 1./std::sqrt(mu_W());
}
double pair_bin_summary::kappa() const
{
	double gms = _magf_gamma_mean_square();
	brgastro::fixbad(gms);
	return 1.-std::sqrt(1/mu_hat()+gms);
}
double pair_bin_summary::kappa_stderr() const
{
	double gms = _magf_gamma_mean_square();
	double gserr = _magf_gamma_square_stderr();
	brgastro::fixbad(gms);
	brgastro::fixbad(gserr);
	return sqrt_err(1/mu_hat()+gms,quad_add(mu_stderr(),gserr));
}

BRG_UNITS pair_bin_summary::Sigma() const
{
	BRG_UNITS result = kappa() * magf_sigma_crit();
	return result;
}
BRG_UNITS pair_bin_summary::Sigma_stderr() const
{
	BRG_UNITS result = kappa_stderr() * magf_sigma_crit();
	return result;
}

double pair_bin_summary::model_mu(const double & MLratio) const
{
	double result = 1./(square(1-model_kappa(MLratio))+square(model_gamma_t(MLratio)));
	return result;
}
double pair_bin_summary::model_kappa(const double & MLratio) const
{
	double result = model_Sigma(MLratio)/magf_sigma_crit();
	return result;
}
BRG_UNITS pair_bin_summary::model_Sigma(const double & MLratio) const
{
	BRG_UNITS result = _magf_model_profile(MLratio).proj_dens(magf_R_mean());
	return result;
}

#endif // Magnification

#endif

} // end namespace brgastro
