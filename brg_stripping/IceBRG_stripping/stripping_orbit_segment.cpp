/**********************************************************************\
  @file stripping_orbit_segment.cpp

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

#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

#include "IceBRG_main/common.hpp"

#include "IceBRG_main/file_access/ascii_table.hpp"

#include "IceBRG_main/math/interpolator/interpolator.hpp"
#include "IceBRG_main/math/interpolator/interpolator_derivative.hpp"
#include "IceBRG_main/math/solvers/solvers.hpp"

#include "IceBRG_main/utility.hpp"

#include "gabdt.h"
#include "solve_rt_functors.h"

#include "stripping_orbit_segment.h"
#include "rm_conversions.hpp"

using namespace IceBRG;

// IceBRG::stripping_orbit_segment class method implementations
#if (1)

double IceBRG::stripping_orbit_segment::_tidal_strip_retained( const density_profile *host_group,
		const density_profile *satellite, distance_type const & r,
		velocity_type const & vr, velocity_type const & vt,
		time_type const & time_step, const long double &sum_delta_rho ) const
{
	distance_type new_rt;
	time_type inst_tangential_orbital_period, inst_full_orbital_period,
		inst_orbital_period, hm_period, stripping_period;
	double mass_frac_retained, mass_frac_lost_total;
	inst_tangential_orbital_period = 2 * pi * r / safe_d(vt);
	inst_full_orbital_period = 2 * pi * r / safe_d(quad_add(vt,vr));

	inst_orbital_period = inst_tangential_orbital_period *
			std::pow(inst_full_orbital_period/safe_d(inst_tangential_orbital_period),
				_tidal_stripping_radialness_);

	hm_period = satellite->othm();

	if(isbad(hm_period) or (hm_period<=0))
		stripping_period = inst_orbital_period;
	else
		stripping_period = inst_orbital_period *
			std::pow(hm_period/safe_d(inst_orbital_period), _tidal_stripping_deceleration_);

	new_rt = _get_rt( host_group, satellite, r, vr, vt, time_step, sum_delta_rho );

	if ( !( new_rt > 0 ) )
		mass_frac_lost_total = 0;
	else
	{
		mass_frac_lost_total = max(
				1. - satellite->enc_mass( new_rt ) / safe_d(satellite->mtot()), 0. );
	}
	mass_frac_retained = max(
			min(
					1.
							- mass_frac_lost_total * time_step / stripping_period
									* _tidal_stripping_amplification_, 1 ), 0. );

	return mass_frac_retained;
}

distance_type IceBRG::stripping_orbit_segment::_get_rt( const density_profile *host_group,
		const density_profile *satellite, distance_type const & r,
		velocity_type const & vr, velocity_type const & vt,
		time_type const & time_step, const long double &sum_delta_rho,
		const bool silent ) const
{
	flt_t omega;
	distance_type new_rt, old_rt;
	flt_t max_rt = 0;

	omega = vt / safe_d( r );

	// Check for null case
	if ( satellite->mtot() <= 0 )
	{
		return 0;
	}

	solve_rt_grid_functor rt_grid_solver( omega, satellite,
			host_group->Daccel( r ), sum_delta_rho );
	solve_rt_it_functor rt_it_solver( omega, satellite,
			host_group->Daccel( r ), sum_delta_rho );

	if ( satellite->mtot() <= 0 )
		return 0;

	old_rt = satellite->rt();

	// First, we try solving iteratively
	try
	{
	    new_rt = solve_iterate( rt_it_solver, old_rt, 1, 0.0001, 1000 );
	}
	catch(const std::exception & e)
	{
	    new_rt = 0; // signal that it failed.
	}
	if ( ( new_rt == 0 ) || ( isbad( new_rt ) ) )
	{
		// Iteratively didn't work, so we go to the grid option

		max_rt = 2 * default_tau_factor * satellite->rvir();
		old_rt = new_rt;
		try
		{
			new_rt = solve_grid( rt_grid_solver, (flt_t)0., max_rt, 100, (flt_t)0.);
		}
		catch(const std::exception &e)
		{
			if ( !silent )
				std::cerr << "WARNING: Could not solve rt:\n" << e.what() << std::endl;
			new_rt = 0; // Most likely value in the case where we can't solve it
		}
	}

	return ( max( new_rt, 0 ) );
}

IceBRG::stripping_orbit_segment::stripping_orbit_segment()
{
	_init();
}

void IceBRG::stripping_orbit_segment::_init()
{
	_current_satellite_in_use_ = false;
	_current_host_in_use_ = false;
	_host_loaded_ = false;
	_satellite_loaded_ = false;
	return clear();
}

// Swap functions
void IceBRG::stripping_orbit_segment::swap(stripping_orbit_segment &other)
{
	using std::swap;

	// Integration parameters
#if(1)
	swap(_spline_resolution_, other._spline_resolution_);
	swap(_interpolation_type_, other._interpolation_type_);
	swap(_v_0_, other._v_0_);
	swap(_r_0_, other._r_0_);
	swap(_step_length_power_, other._step_length_power_);
	swap(_step_factor_max_, other._step_factor_max_);
	swap(_step_factor_min_, other._step_factor_min_);
#endif

	// Tuning values
#if(1)
	swap(_tidal_stripping_amplification_, other._tidal_stripping_amplification_);
	swap(_tidal_stripping_deceleration_, other._tidal_stripping_deceleration_);
	swap(_tidal_stripping_radialness_, other._tidal_stripping_radialness_);
	swap(_tidal_shocking_amplification_, other._tidal_shocking_amplification_);
	swap(_tidal_shocking_persistance_, other._tidal_shocking_persistance_);
	swap(_tidal_shocking_power_, other._tidal_shocking_power_);
#endif

	swap(_rt_list_, other._rt_list_);
	swap(_rt_ratio_list_, other._rt_ratio_list_);
	swap(_delta_rho_list_, other._delta_rho_list_);
	swap(_sum_delta_rho_list_, other._sum_delta_rho_list_);
	swap(_x_data_, other._x_data_);
	swap(_y_data_, other._y_data_);
	swap(_z_data_, other._z_data_);
	swap(_vx_data_, other._vx_data_);
	swap(_vy_data_, other._vy_data_);
	swap(_vz_data_, other._vz_data_);
	swap(_t_min_natural_value_, other._t_min_natural_value_);
	swap(_t_max_natural_value_, other._t_max_natural_value_);
	swap(_satellite_parameter_data_, other._satellite_parameter_data_);
	swap(_host_parameter_data_, other._host_parameter_data_);
	swap(_num_parameters_, other._num_parameters_);
	swap(_m_ret_list_, other._m_ret_list_);
	swap(_m_vir_ret_list_, other._m_vir_ret_list_);
	swap(_satellite_parameter_unitconvs_,
			other._satellite_parameter_unitconvs_);
	swap(_satellite_output_parameters_,
			other._satellite_output_parameters_);
	swap(_host_parameter_unitconvs_, other._host_parameter_unitconvs_);
	swap(_host_output_parameters_, other._host_output_parameters_);
	swap(_gabdt_list_, other._gabdt_list_);
	swap(_sum_gabdt_list_, other._sum_gabdt_list_);
	swap(_phase_list_, other._phase_list_);
	swap(_phase_output_list_, other._phase_output_list_);
	swap(_record_full_data_, other._record_full_data_);
	swap(_x_spline_, other._x_spline_);
	swap(_y_spline_, other._y_spline_);
	swap(_z_spline_, other._z_spline_);
	swap(_vx_spline_, other._vx_spline_);
	swap(_vy_spline_, other._vy_spline_);
	swap(_vz_spline_, other._vz_spline_);
	swap(_test_mass_spline_, other._test_mass_spline_);
	swap(_host_parameter_splines_, other._host_parameter_splines_);
	swap(_init_host_ptr_, other._init_host_ptr_);
	swap(_init_satellite_ptr_, other._init_satellite_ptr_);
	swap(_init_sum_delta_rho_, other._init_sum_delta_rho_);
	swap(_init_sum_gabdt_, other._init_sum_gabdt_);
	swap(_calculated_, other._calculated_);
	swap(_bad_result_, other._bad_result_);
	swap(_current_satellite_in_use_, other._current_satellite_in_use_);
	swap(_current_satellite_ptr_, other._current_satellite_ptr_);
	swap(_current_host_in_use_, other._current_host_in_use_);
	swap(_current_host_ptr_, other._current_host_ptr_);
	swap(_host_loaded_, other._host_loaded_);
	swap(_satellite_loaded_, other._satellite_loaded_);
	swap(_using_private_init_host_, other._using_private_init_host_);
	swap(_using_private_init_satellite_,
			other._using_private_init_satellite_);
	swap(_private_tNFW_init_host_, other._private_tNFW_init_host_);
	swap(_private_tNFW_init_satellite_,
			other._private_tNFW_init_satellite_);
	swap(_evolving_host_, other._evolving_host_);
	swap(_override_t_max_, other._override_t_max_);
	swap(_override_t_min_, other._override_t_min_);
	swap(_t_max_override_val_, other._t_max_override_val_);
	swap(_t_min_override_val_, other._t_min_override_val_);
	swap(_likely_disrupted_, other._likely_disrupted_);

	if ( _using_private_init_host_ )
	{
		_init_host_ptr_ = &_private_tNFW_init_host_;
	}
	if ( other._using_private_init_host_ )
	{
		other._init_host_ptr_ = &other._private_tNFW_init_host_;
	}
	if ( _using_private_init_satellite_ )
	{
		_init_satellite_ptr_ = &_private_tNFW_init_satellite_;
	}
	if ( other._using_private_init_satellite_ )
	{
		other._init_satellite_ptr_ = &other._private_tNFW_init_satellite_;
	}
}
namespace std
{
	template <>
	void swap(IceBRG::stripping_orbit_segment &same,
			IceBRG::stripping_orbit_segment &other)
	{
		same.swap(other);
	}
}
void swap(IceBRG::stripping_orbit_segment &same,
		IceBRG::stripping_orbit_segment &other)
{
	same.swap(other);
}

IceBRG::stripping_orbit_segment::stripping_orbit_segment(
		const stripping_orbit_segment &other )
{

	// Integration parameters
#if(1)
	_spline_resolution_ = other._spline_resolution_;
	_interpolation_type_ = other._interpolation_type_;
	_v_0_ = other._v_0_;
	_r_0_ = other._r_0_;
	_step_length_power_ = other._step_length_power_;
	_step_factor_max_ = other._step_factor_max_;
	_step_factor_min_ = other._step_factor_min_;
#endif

	// Tuning values
#if(1)
	_tidal_stripping_amplification_ = other._tidal_stripping_amplification_;
	_tidal_stripping_deceleration_ = other._tidal_stripping_deceleration_;
	_tidal_stripping_radialness_ = other._tidal_stripping_radialness_;
	_tidal_shocking_amplification_ = other._tidal_shocking_amplification_;
	_tidal_shocking_persistance_ = other._tidal_shocking_persistance_;
	_tidal_shocking_power_ = other._tidal_shocking_power_;
#endif

	_rt_list_ = other._rt_list_;
	_rt_ratio_list_ = other._rt_ratio_list_;
	_delta_rho_list_ = other._delta_rho_list_;
	_sum_delta_rho_list_ = other._sum_delta_rho_list_;
	_x_data_ = other._x_data_;
	_y_data_ = other._y_data_;
	_z_data_ = other._z_data_;
	_vx_data_ = other._vx_data_;
	_vy_data_ = other._vy_data_;
	_vz_data_ = other._vz_data_;
	_t_min_natural_value_ = other._t_min_natural_value_;
	_t_max_natural_value_ = other._t_max_natural_value_;
	_satellite_parameter_data_ = other._satellite_parameter_data_;
	_host_parameter_data_ = other._host_parameter_data_;
	_num_parameters_ = other._num_parameters_;
	_m_vir_ret_list_ = other._m_vir_ret_list_;
	_satellite_parameter_unitconvs_ =
			other._satellite_parameter_unitconvs_;
	_satellite_output_parameters_ =
			other._satellite_output_parameters_;
	_host_parameter_unitconvs_ = other._host_parameter_unitconvs_;
	_host_output_parameters_ = other._host_output_parameters_;
	_gabdt_list_ = other._gabdt_list_;
	_sum_gabdt_list_ = other._sum_gabdt_list_;
	_phase_list_ = other._phase_list_;
	_phase_output_list_ = other._phase_output_list_;
	_record_full_data_ = other._record_full_data_;
	_x_spline_ = other._x_spline_;
	_y_spline_ = other._y_spline_;
	_z_spline_ = other._z_spline_;
	_vx_spline_ = other._vx_spline_;
	_vy_spline_ = other._vy_spline_;
	_vz_spline_ = other._vz_spline_;
	_test_mass_spline_ = other._test_mass_spline_;
	_host_parameter_splines_ = other._host_parameter_splines_;
	_init_host_ptr_ = other._init_host_ptr_;
	_init_satellite_ptr_ = other._init_satellite_ptr_;
	_init_sum_delta_rho_ = other._init_sum_delta_rho_;
	_init_sum_gabdt_ = other._init_sum_gabdt_;
	_calculated_ = other._calculated_;
	_bad_result_ = other._bad_result_;
	_current_satellite_in_use_ = other._current_satellite_in_use_;
	_current_host_in_use_ = other._current_host_in_use_;
	_host_loaded_ = other._host_loaded_;
	_satellite_loaded_ = other._satellite_loaded_;
	_using_private_init_host_ = other._using_private_init_host_;
	_using_private_init_satellite_ =
			other._using_private_init_satellite_;
	_private_tNFW_init_host_ = other._private_tNFW_init_host_;
	_private_tNFW_init_satellite_ =
			other._private_tNFW_init_satellite_;
	_evolving_host_ = other._evolving_host_;
	_override_t_max_ = other._override_t_max_;
	_override_t_min_ = other._override_t_min_;
	_t_max_override_val_ = other._t_max_override_val_;
	_t_min_override_val_ = other._t_min_override_val_;
	_likely_disrupted_ = other._likely_disrupted_;

	if ( _current_satellite_in_use_ )
	{
		_current_satellite_ptr_ =
				other._current_satellite_ptr_->density_profile_clone();
	}
	else
	{
		_current_satellite_ptr_ = NULL;
	}
	if ( _current_host_in_use_ )
	{
		_current_host_ptr_ =
				other._current_host_ptr_->density_profile_clone();
	}
	else
	{
		_current_host_ptr_ = NULL;
	}
	if ( _using_private_init_host_ )
	{
		_init_host_ptr_ = &_private_tNFW_init_host_;
	}
	if ( _using_private_init_satellite_ )
	{
		_init_satellite_ptr_ = &_private_tNFW_init_satellite_;
	}
}

IceBRG::stripping_orbit_segment & IceBRG::stripping_orbit_segment::operator=(
		stripping_orbit_segment other )
{
	swap(other);
	return *this;
}

IceBRG::stripping_orbit_segment::stripping_orbit_segment(
		const density_profile *init_init_host,
		const density_profile *init_init_satellite, const ssize_t init_resolution )
{
	_init();
	_init_host_ptr_ = init_init_host;
	_init_satellite_ptr_ = init_init_satellite;
	_spline_resolution_ = init_resolution;
	_current_satellite_ptr_ = _init_satellite_ptr_->density_profile_clone();
	_current_satellite_in_use_ = true;
	_current_host_ptr_ = _init_host_ptr_->density_profile_clone();
	_current_host_in_use_ = true;
	_host_loaded_ = true;
	_satellite_loaded_ = true;

}

IceBRG::stripping_orbit_segment::~stripping_orbit_segment()
{
	if ( _current_satellite_in_use_ )
	{
		delete _current_satellite_ptr_;
		_current_satellite_ptr_ = NULL;
	}
	if ( _current_host_in_use_ )
	{
		delete _current_host_ptr_;
		_current_host_ptr_ = NULL;
	}
}

IceBRG::stripping_orbit_segment *IceBRG::stripping_orbit_segment::stripping_orbit_spline_clone() const
{
	return new stripping_orbit_segment( *this );
}

void IceBRG::stripping_orbit_segment::clear()
{
	clear_calcs();

	// Integration parameters
#if(1)
	_spline_resolution_ = IceBRG::stripping_orbit::default_spline_resolution();
	_interpolation_type_ = IceBRG::stripping_orbit::default_interpolation_type();
	_v_0_ = IceBRG::stripping_orbit::default_spline_resolution();
	_r_0_ = IceBRG::stripping_orbit::default_spline_resolution();
	_step_length_power_ = IceBRG::stripping_orbit::default_spline_resolution();
	_step_factor_max_ = IceBRG::stripping_orbit::default_spline_resolution();
	_step_factor_min_ = IceBRG::stripping_orbit::default_spline_resolution();
#endif

	// Tuning values
#if(1)
	_tidal_stripping_amplification_ = IceBRG::stripping_orbit::default_tidal_stripping_amplification();
	_tidal_stripping_deceleration_ = IceBRG::stripping_orbit::default_tidal_stripping_deceleration();
	_tidal_stripping_radialness_ = IceBRG::stripping_orbit::default_tidal_stripping_radialness();
	_tidal_shocking_amplification_ = IceBRG::stripping_orbit::default_tidal_shocking_amplification();
	_tidal_shocking_persistance_ = IceBRG::stripping_orbit::default_tidal_shocking_persistance();
	_tidal_shocking_power_ = IceBRG::stripping_orbit::default_tidal_shocking_power();
#endif


	_phase_list_.clear();

	_x_data_.clear();
	_y_data_.clear();
	_z_data_.clear();
	_vx_data_.clear();
	_vy_data_.clear();
	_vz_data_.clear();

	_satellite_parameter_data_.clear();
	_host_parameter_data_.clear();
	_satellite_parameter_unitconvs_.clear();
	_satellite_output_parameters_.clear();
	_host_parameter_unitconvs_.clear();
	_host_output_parameters_.clear();

	_x_spline_.clear();
	_y_spline_.clear();
	_z_spline_.clear();
	_vx_spline_.clear();
	_vy_spline_.clear();
	_vz_spline_.clear();
	_test_mass_spline_.clear();
	_host_parameter_splines_.clear();

	_t_min_natural_value_ = std::numeric_limits<double>::max();
	_t_max_natural_value_ = ( -std::numeric_limits<double>::max() );
	_num_parameters_ = 0;
	_record_full_data_ = false;

	_init_sum_delta_rho_ = 0;
	_init_sum_gabdt_.override_zero();

	_init_host_ptr_ = NULL;
	_init_satellite_ptr_ = NULL;

	if ( _current_satellite_in_use_ )
	{
		delete _current_satellite_ptr_;
		_current_satellite_ptr_ = NULL;
	}
	if ( _current_host_in_use_ )
	{
		delete _current_host_ptr_;
		_current_host_ptr_ = NULL;
	}

	_current_satellite_ptr_ = NULL;
	_current_satellite_in_use_ = false;
	_current_host_ptr_ = NULL;
	_current_host_in_use_ = false;
	_host_loaded_ = false;
	_satellite_loaded_ = false;
	_using_private_init_host_ = false;
	_using_private_init_satellite_ = false;
	_evolving_host_ = false;
}

void IceBRG::stripping_orbit_segment::clear_calcs() const
{
	_rt_list_.clear();
	_rt_ratio_list_.clear();
	_phase_output_list_.clear();
	_m_ret_list_.clear();
	_m_vir_ret_list_.clear();
	_delta_rho_list_.clear();
	_sum_delta_rho_list_.clear();
	_gabdt_list_.clear();
	_sum_gabdt_list_.clear();

	_calculated_ = false;
	_bad_result_ = false;
	_likely_disrupted_ = false;
}

// Setting integration parameters
#if(1)
void IceBRG::stripping_orbit_segment::set_resolution( const ssize_t new_resolution,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_resolution == _spline_resolution_ )
		return;

	if ( new_resolution < 2 )
	{
		throw std::logic_error("Attempt to set resolution to value below minimum of 2.");
	}
	clear_calcs();
	_spline_resolution_ = new_resolution;
}
void IceBRG::stripping_orbit_segment::set_interpolation_type(
		const IceBRG::stripping_orbit::allowed_interpolation_type new_type,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( _interpolation_type_ == new_type )
		return;
	clear_calcs();
	_interpolation_type_ = new_type;
}
void IceBRG::stripping_orbit_segment::set_v_0( const double new_v_0,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_v_0 == _v_0_ )
		return;

	if ( new_v_0 <= 0 )
	{

		throw std::logic_error("Attempt to set v_0 to value at or below minimum of 0.");
	}
	clear_calcs();
	_v_0_ = new_v_0;
}
void IceBRG::stripping_orbit_segment::set_r_0( const double new_r_0,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_r_0 == _r_0_ )
		return;

	if ( new_r_0 <= 0 )
	{

		throw std::logic_error("Attempt to set r_0 to value at or below minimum of 0.");
	}
	clear_calcs();
	_r_0_ = new_r_0;
}
void IceBRG::stripping_orbit_segment::set_step_length_power( const double new_step_length_power,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_step_length_power == _step_length_power_ )
		return;

	clear_calcs();
	_step_length_power_ = new_step_length_power;
}
void IceBRG::stripping_orbit_segment::set_step_factor_max( const double new_step_factor_max,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_step_factor_max == _step_factor_max_ )
		return;

	if ( new_step_factor_max < 1 )
	{

		throw std::logic_error("Attempt to set step_factor_max to value below minimum of 1.");
	}
	clear_calcs();
	_step_factor_max_ = new_step_factor_max;
}
void IceBRG::stripping_orbit_segment::set_step_factor_min( const double new_step_factor_min,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_step_factor_min == _step_factor_min_ )
		return;

	if ( new_step_factor_min > 1 )
	{

		throw std::logic_error("Attempt to set step_factor_min to value above minimum of 1.");
	}
	clear_calcs();
	_step_factor_min_ = new_step_factor_min;
}
#endif

// Setting tuning parameters
#if(1)

void IceBRG::stripping_orbit_segment::set_tidal_stripping_amplification(
		const double new_tidal_stripping_amplification,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_stripping_amplification == _tidal_stripping_amplification_ )
		return;

	if ( new_tidal_stripping_amplification < 0 )
	{

		throw std::logic_error("Attempt to set tidal_stripping_amplification to value below minimum of 0.");
	}

	clear_calcs();
	_tidal_stripping_amplification_ = new_tidal_stripping_amplification;
}
void IceBRG::stripping_orbit_segment::set_tidal_stripping_deceleration(
		const double new_tidal_stripping_deceleration,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_stripping_deceleration == _tidal_stripping_deceleration_ )
		return;

	clear_calcs();
	_tidal_stripping_deceleration_ = new_tidal_stripping_deceleration;
}
void IceBRG::stripping_orbit_segment::set_tidal_stripping_radialness(
		const double new_tidal_stripping_radialness,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_stripping_radialness == _tidal_stripping_radialness_ )
		return;

	clear_calcs();
	_tidal_stripping_radialness_ = new_tidal_stripping_radialness;
}
void IceBRG::stripping_orbit_segment::set_tidal_shocking_amplification(
		const double new_tidal_shocking_amplification,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_shocking_amplification == _tidal_shocking_amplification_ )
		return;

	if ( new_tidal_shocking_amplification < 0 )
	{

		throw std::logic_error("Attempt to set tidal_shocking_amplification to value below minimum of 0.");
	}

	clear_calcs();
	_tidal_shocking_amplification_ = new_tidal_shocking_amplification;
}
void IceBRG::stripping_orbit_segment::set_tidal_shocking_persistance(
		const double new_tidal_shocking_persistance,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_shocking_persistance == _tidal_shocking_persistance_ )
		return;

	if ( new_tidal_shocking_persistance <= 0 )
	{

		throw std::logic_error("Attempt to set tidal_shocking_persistance to value at or below minimum of 0.");
	}
	clear_calcs();
	_tidal_shocking_persistance_ = new_tidal_shocking_persistance;
}
void IceBRG::stripping_orbit_segment::set_tidal_shocking_power(
		const double new_tidal_shocking_power,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_shocking_power == _tidal_shocking_power_ )
		return;

	clear_calcs();
	_tidal_shocking_power_ = new_tidal_shocking_power;
}
#endif


void IceBRG::stripping_orbit_segment::add_point( distance_type const & x,
		distance_type const & y, distance_type const & z, time_type const & t,
		const double new_mass )
{
	_calculated_ = false;
	_x_spline_.add_point( t, x );
	_y_spline_.add_point( t, y );
	_z_spline_.add_point( t, z );
	_vx_spline_.add_unknown_point( t );
	_vy_spline_.add_unknown_point( t );
	_vz_spline_.add_unknown_point( t );
	_test_mass_spline_.add_point( t, new_mass );
	if ( t < _t_min_natural_value_ )
		_t_min_natural_value_ = t;
	if ( t > _t_max_natural_value_ )
		_t_max_natural_value_ = t;
}

void IceBRG::stripping_orbit_segment::add_point( distance_type const & x,
		distance_type const & y, distance_type const & z, velocity_type const & vx,
		velocity_type const & vy, velocity_type const & vz, time_type const & t,
		const double new_test_mass )
{
	_calculated_ = false;
	_x_spline_.add_point( t, x );
	_y_spline_.add_point( t, y );
	_z_spline_.add_point( t, z );
	_vx_spline_.add_point( t, vx );
	_vy_spline_.add_point( t, vy );
	_vz_spline_.add_point( t, vz );
	_test_mass_spline_.add_point( t, new_test_mass );
	if ( t < _t_min_natural_value_ )
		_t_min_natural_value_ = t;
	if ( t > _t_max_natural_value_ )
		_t_max_natural_value_ = t;
}

void IceBRG::stripping_orbit_segment::add_x_point(
		distance_type const & x, time_type const & t )
{
	_calculated_ = false;

	_x_spline_.add_point( t, x );
	if ( t < _t_min_natural_value_ )
		_t_min_natural_value_ = t;
	if ( t > _t_max_natural_value_ )
		_t_max_natural_value_ = t;

}

void IceBRG::stripping_orbit_segment::add_y_point(
		distance_type const & y, time_type const & t )
{
	_calculated_ = false;
	_y_spline_.add_point( t, y );
	if ( t < _t_min_natural_value_ )
		_t_min_natural_value_ = t;
	if ( t > _t_max_natural_value_ )
		_t_max_natural_value_ = t;
}

void IceBRG::stripping_orbit_segment::add_z_point(
		distance_type const & z, time_type const & t )
{
	_calculated_ = false;
	_z_spline_.add_point( t, z );
	if ( t < _t_min_natural_value_ )
		_t_min_natural_value_ = t;
	if ( t > _t_max_natural_value_ )
		_t_max_natural_value_ = t;
}

void IceBRG::stripping_orbit_segment::add_vx_point(
		velocity_type const & vx, time_type const & t )
{
	_calculated_ = false;
	_vx_spline_.add_point( t, vx );
	if ( t < _t_min_natural_value_ )
		_t_min_natural_value_ = t;
	if ( t > _t_max_natural_value_ )
		_t_max_natural_value_ = t;
}

void IceBRG::stripping_orbit_segment::add_vy_point(
		velocity_type const & vy, time_type const & t )
{
	_calculated_ = false;
	_vy_spline_.add_point( t, vy );
	if ( t < _t_min_natural_value_ )
		_t_min_natural_value_ = t;
	if ( t > _t_max_natural_value_ )
		_t_max_natural_value_ = t;
}

void IceBRG::stripping_orbit_segment::add_vz_point(
		velocity_type const & vz, time_type const & t )
{
	_calculated_ = false;
	_vz_spline_.add_point( t, vz );
	if ( t < _t_min_natural_value_ )
		_t_min_natural_value_ = t;
	if ( t > _t_max_natural_value_ )
		_t_max_natural_value_ = t;
}

void IceBRG::stripping_orbit_segment::add_unknown_vx_point(
		time_type const & t )
{
	_calculated_ = false;
	_vx_spline_.add_unknown_point( t );
	if ( t < _t_min_natural_value_ )
		_t_min_natural_value_ = t;
	if ( t > _t_max_natural_value_ )
		_t_max_natural_value_ = t;
}

void IceBRG::stripping_orbit_segment::add_unknown_vy_point(
		time_type const & t )
{
	_calculated_ = false;
	_vy_spline_.add_unknown_point( t );
	if ( t < _t_min_natural_value_ )
		_t_min_natural_value_ = t;
	if ( t > _t_max_natural_value_ )
		_t_max_natural_value_ = t;
}

void IceBRG::stripping_orbit_segment::add_unknown_vz_point(
		time_type const & t )
{
	_calculated_ = false;
	_vz_spline_.add_unknown_point( t );
	if ( t < _t_min_natural_value_ )
		_t_min_natural_value_ = t;
	if ( t > _t_max_natural_value_ )
		_t_max_natural_value_ = t;
}

void IceBRG::stripping_orbit_segment::add_test_mass_point(
		const double test_mass, time_type const & t )
{
	_calculated_ = false;
	_test_mass_spline_.add_point( t, test_mass );
	if ( t < _t_min_natural_value_ )
		_t_min_natural_value_ = t;
	if ( t > _t_max_natural_value_ )
		_t_max_natural_value_ = t;
}

void IceBRG::stripping_orbit_segment::add_host_parameter_point(
		const std::vector< flt_t > & parameters, time_type const & t,
		const bool silent )
{
	const ssize_t num_parameters = ssize(parameters);

	if ( ssize(_host_parameter_splines_) == 0 )
		_host_parameter_splines_.resize( num_parameters );

	if ( ssize(_host_parameter_splines_) != num_parameters )
	{
		throw std::logic_error("All parameter lists passed to stripping_orbit_segment must have same size.");
	}

	for ( ssize_t i = 0; i < num_parameters; i++ )
	{
		_host_parameter_splines_[i].add_point( t, parameters[i] );
	}

	if ( ssize(_host_parameter_splines_[0]) >= 2 )
		_evolving_host_ = true;
	if ( t < _t_min_natural_value_ )
		_t_min_natural_value_ = t;
	if ( t > _t_max_natural_value_ )
		_t_max_natural_value_ = t;
}

void IceBRG::stripping_orbit_segment::_reserve( const ssize_t n,
		const bool silent ) const
{
	if ( n < 1 )
	{
		throw std::logic_error("Attempt to reserve length of 0 or less in stripping_orbit_segment::reserve.");
	}
	_phase_list_.reserve( n );
	_m_ret_list_.reserve( n );
	_m_vir_ret_list_.reserve( n );
	_phase_output_list_.reserve( n );
	_rt_list_.reserve( n );
	_rt_ratio_list_.reserve( n );
	_delta_rho_list_.reserve( n );
	_sum_delta_rho_list_.reserve( n );
	_gabdt_list_.reserve( n );
	_sum_gabdt_list_.reserve( n );
}

void IceBRG::stripping_orbit_segment::set_init_host(
		const density_profile *new_init_host )
{
	_init_host_ptr_ = new_init_host;
	if ( _current_host_in_use_ )
	{
		delete _current_host_ptr_;
		_current_host_ptr_ = NULL;
	}
	_current_host_ptr_ = _init_host_ptr_->density_profile_clone();
	_current_host_in_use_ = true;
	_host_loaded_ = true;
	_using_private_init_host_ = false;
	_calculated_ = false;
}

void IceBRG::stripping_orbit_segment::set_init_satellite(
		const density_profile *new_init_satellite )
{
	_init_satellite_ptr_ = new_init_satellite;
	if ( _current_satellite_in_use_ )
	{
		delete _current_satellite_ptr_;
		_current_satellite_ptr_ = NULL;
	}
	_current_satellite_ptr_ = _init_satellite_ptr_->density_profile_clone();
	_current_satellite_in_use_ = true;
	_satellite_loaded_ = true;
	_using_private_init_satellite_ = false;
	_calculated_ = false;
}

void IceBRG::stripping_orbit_segment::set_init_sum_deltarho(
		const long double &new_init_sum_deltarho )
{
	_init_sum_delta_rho_ = new_init_sum_deltarho;
}

void IceBRG::stripping_orbit_segment::set_init_sum_gabdt(
		const gabdt &new_init_sum_gabdt )
{
	_init_sum_gabdt_ = new_init_sum_gabdt;
}

void IceBRG::stripping_orbit_segment::set_tNFW_init_satellite(
		mass_type const & new_init_mvir0, const double z,
		const double new_init_c, const double new_init_tau )
{
	_using_private_init_satellite_ = true;
	_private_tNFW_init_satellite_ = tNFW_profile( new_init_mvir0, z,
			new_init_c, new_init_tau );

	// Correct mass for r200 -> rvir
	_private_tNFW_init_satellite_.set_mvir(new_init_mvir0 /
	        get_M_ratio(_private_tNFW_init_satellite_.z(),
	                    _private_tNFW_init_satellite_.c()));

	_init_satellite_ptr_ = &_private_tNFW_init_satellite_;

	if ( _current_satellite_in_use_ )
	{
		delete _current_satellite_ptr_;
		_current_satellite_ptr_ = NULL;
	}
	_current_satellite_ptr_ = _init_satellite_ptr_->density_profile_clone();
	_current_satellite_in_use_ = true;

	if ( _using_private_init_host_ )
	{
		if ( _private_tNFW_init_host_.z() != z )
			_private_tNFW_init_host_.set_z( z );
	}

	_satellite_loaded_ = true;
	_calculated_ = false;
}

void IceBRG::stripping_orbit_segment::set_tNFW_host(
		mass_type const & new_init_mvir0, const double z,
		const double new_init_c, const double new_init_tau )
{
	_using_private_init_host_ = true;
	_private_tNFW_init_host_ = tNFW_profile( new_init_mvir0, z, new_init_c,
			new_init_tau );

    // Correct mass for r200 -> rvir
	_private_tNFW_init_host_.set_mvir(new_init_mvir0 /
            get_M_ratio(_private_tNFW_init_host_.z(),
                        _private_tNFW_init_host_.c()));

	_init_host_ptr_ = &_private_tNFW_init_host_;

	if ( _using_private_init_satellite_ )
	{
		if ( _private_tNFW_init_satellite_.z() != z )
		{
			_private_tNFW_init_satellite_.set_z( z );

			if ( _current_satellite_in_use_ )
			{
				delete _current_satellite_ptr_;
				_current_satellite_ptr_ = NULL;
			}
			_current_satellite_ptr_ =
					_init_satellite_ptr_->density_profile_clone();
			_current_satellite_in_use_ = true;
		}
	}

	_host_loaded_ = true;
	_calculated_ = false;
}

void IceBRG::stripping_orbit_segment::set_t_min(
		time_type const & new_t_min )
{
	_t_min_override_val_ = new_t_min;
	_override_t_min_ = true;
}

void IceBRG::stripping_orbit_segment::set_t_max(
		time_type const & new_t_max )
{
	_t_max_override_val_ = new_t_max;
	_override_t_max_ = true;
}

void IceBRG::stripping_orbit_segment::reset_t_min()
{
	_t_min_natural_value_ = std::numeric_limits<double>::max();
	for ( ssize_t i = 0; i < ssize(_x_spline_); i++ )
	{
		if ( _x_spline_.sorted_data().at( i ).first < _t_min_natural_value_ )
			_t_min_natural_value_ = _x_spline_.sorted_data().at( i ).first;
	}
	_override_t_min_ = false;
}

void IceBRG::stripping_orbit_segment::reset_t_max()
{
	_t_max_natural_value_ = ( -std::numeric_limits<double>::max() );
	for ( ssize_t i = 0; i < ssize(_x_spline_); i++ )
	{
		if ( _x_spline_.sorted_data().at( i ).first > _t_max_natural_value_ )
			_t_max_natural_value_ = _x_spline_.sorted_data().at( i ).first;
	}
	_override_t_min_ = false;
}

void IceBRG::stripping_orbit_segment::clear_points()
{
	_x_spline_.clear();
	_y_spline_.clear();
	_z_spline_.clear();
	_vx_spline_.clear();
	_vy_spline_.clear();
	_vz_spline_.clear();
	_test_mass_spline_.clear();
	_calculated_ = false;
}
void IceBRG::stripping_orbit_segment::clear_host_parameter_points()
{
	_host_parameter_splines_.clear();
	_calculated_ = false;
}
void IceBRG::stripping_orbit_segment::clear_init_sum_deltarho()
{
	_init_sum_delta_rho_ = 0;
	_calculated_ = false;
}
void IceBRG::stripping_orbit_segment::clear_init_sum_gabdt()
{
	_init_sum_gabdt_.clear();
	_init_sum_gabdt_.override_zero();
	_calculated_ = false;
}
void IceBRG::stripping_orbit_segment::clear_init_satellite()
{
	_init_satellite_ptr_ = NULL;
	if ( _current_satellite_in_use_ )
	{
		delete _current_satellite_ptr_;
		_current_satellite_ptr_ = NULL;
	}
	_current_satellite_ptr_ = 0;
	_current_satellite_in_use_ = false;
	_satellite_loaded_ = false;
	_using_private_init_satellite_ = false;
	_calculated_ = false;
}
void IceBRG::stripping_orbit_segment::clear_init_host()
{
	_init_host_ptr_ = NULL;
	if ( _current_host_in_use_ )
	{
		delete _current_host_ptr_;
		_current_host_ptr_ = NULL;
	}
	_current_host_ptr_ = NULL;
	_current_host_in_use_ = false;
	_host_loaded_ = false;
	_using_private_init_host_ = false;
	_calculated_ = false;
}

ssize_t IceBRG::stripping_orbit_segment::length() const
{
	ssize_t result = std::numeric_limits<ssize_t>::max();
	if ( ssize(_x_spline_) < result )
		result = ssize(_x_spline_);
	if ( ssize(_y_spline_) < result )
		result = ssize(_y_spline_);
	if ( ssize(_z_spline_) < result )
		result = ssize(_z_spline_);
	if ( ssize(_vx_spline_) < result )
		result = ssize(_x_spline_);
	if ( ssize(_vy_spline_) < result )
		result = ssize(_y_spline_);
	if ( ssize(_vz_spline_) < result )
		result = ssize(_z_spline_);
	return result;
}

// Functions for determining how calc() will be called
void IceBRG::stripping_orbit_segment::set_record_full_data(
		const bool new_record_full_data ) const
{
	// Check if anything is actually changing here
	if ( new_record_full_data == _record_full_data_ )
		return;

	clear_calcs();
	_record_full_data_ = new_record_full_data;
}

// Function to calculate stripping
void IceBRG::stripping_orbit_segment::calc( const bool silent ) const
{
	distance_type r;
	velocity_type v, vt, vr;
	time_type t, t_step, t_min_to_use, t_max_to_use;
	flt_t current_deltarho;
	std::vector< flt_t > satellite_parameters;
	std::vector< flt_t > host_parameters;
	int num_host_parameters = 0;
	int counter = 0;
	double step_length_factor = 1;
	double m_ret = 1;

	if ( ( !_host_loaded_ ) || ( !_satellite_loaded_ ) )
	{
		throw std::logic_error("Host or satellite not loaded in stripping_orbit_segment::calc().");
	}

	if ( length() < 2 )
	{
		throw std::runtime_error("Too few data points to calculate stripping in stripping_orbit_segment::calc().");
	}

	// Initialisation

	if ( _override_t_min_ )
	{
		t_min_to_use = _t_min_override_val_;
	}
	else
	{
		t_min_to_use = _t_min_natural_value_;
	}

	if ( _override_t_max_ )
	{
		t_max_to_use = _t_max_override_val_;
	}
	else
	{
		t_max_to_use = _t_max_natural_value_;
	}

	if( t_max_to_use <= t_min_to_use)
	{
		if(t_max_to_use==t_min_to_use)
			throw std::runtime_error("t_max == t_min for stripping_orbit_segment!\n");
		if(!silent)
			std::cerr << "WARNING: t_max < t_min for stripping_orbit segment. They'll be swapped for\n calculation.\n";
		std::swap(t_min_to_use,t_max_to_use);
	}

	t_step = ( t_max_to_use - t_min_to_use ) / _spline_resolution_;
	gabdt current_gabdt = _init_sum_gabdt_;
	t = 0;
	r = 0;
	v = vt = vr = 0;
	current_deltarho = 0;

	clear_calcs();
	_reserve( 2 * _spline_resolution_ + 1 );

	// Setup

	// Velocity splines setup
	_vx_spline_.set_spline_ptr( &_x_spline_ );
	_vy_spline_.set_spline_ptr( &_y_spline_ );
	_vz_spline_.set_spline_ptr( &_z_spline_ );

	// Set up interpolation type
	_pass_interpolation_type();

	_m_ret_list_.push_back( 1 );
	_m_vir_ret_list_.push_back( 1 );
	_delta_rho_list_.push_back( current_deltarho );
	_sum_delta_rho_list_.push_back( _init_sum_delta_rho_ );

	_sum_gabdt_list_.push_back( current_gabdt );
	current_gabdt.set_host_ptr( _init_host_ptr_ );
	current_gabdt.override_zero();
	_gabdt_list_.push_back( current_gabdt );
	current_gabdt.override_zero();

	// Calculate radius and velocity at init position
	r = IceBRG::dist3d( _x_spline_( t_min_to_use ),
			_y_spline_( t_min_to_use ), _z_spline_( t_min_to_use ) );
	v = IceBRG::dist3d( _vx_spline_( t_min_to_use ),
			_vy_spline_( t_min_to_use ), _vz_spline_( t_min_to_use ) );
	vr = IceBRG::dot_product( _vx_spline_( t_min_to_use ),
			_vy_spline_( t_min_to_use ), _vz_spline_( t_min_to_use ),
			_x_spline_( t_min_to_use ), _y_spline_( t_min_to_use ),
			_z_spline_( t_min_to_use ) ) / r;
	if ( std::fabs( vr ) >= std::fabs( v ) )
		vt = 0;
	else
		vt = IceBRG::quad_sub( v, vr );

	// Make sure current_satellite and current_host are set up correctly
	if ( !_current_satellite_in_use_ )
	{
		_current_satellite_ptr_ =
				_init_satellite_ptr_->density_profile_clone();
		_current_satellite_in_use_ = true;
	}
	if ( !_current_host_in_use_ )
	{
		_current_host_ptr_ = _init_host_ptr_->density_profile_clone();
		_current_host_in_use_ = true;
	}

	step_length_factor = _step_length_factor(v, r);

	_rt_list_.push_back(
			_get_rt( _current_host_ptr_, _current_satellite_ptr_, r,
					vr, vt, t_step, current_deltarho ) );
	_rt_ratio_list_.push_back( _rt_list_.at( 0 ) / _rvir() );

	// Get satellite parameters at init position
	if ( _record_full_data_ )
	{
		satellite_parameters = _current_satellite_ptr_->get_parameters( );
		_satellite_parameter_data_.push_back( satellite_parameters );
		host_parameters = _current_host_ptr_->get_parameters( );
		_host_parameter_data_.push_back( host_parameters );
	}

	const mass_type init_mtot = _current_satellite_ptr_->mtot();
	const mass_type init_mvir = _current_satellite_ptr_->mvir();

	// Loop over time
	for ( t = t_min_to_use; t < t_max_to_use;
			t += t_step * step_length_factor )
	{
		if ( _record_full_data_ )
			_phase_output_list_.push_back(
					phase( _x_spline_( t ), _y_spline_( t ), _z_spline_( t ),
							_vx_spline_( t ), _vy_spline_( t ),
							_vz_spline_( t ), t ) );
		counter += 1;
		// Calculate radius and velocity at this position
		try
		{
			r = IceBRG::dist3d( _x_spline_( t ), _y_spline_( t ),
					_z_spline_( t ) );
			v = IceBRG::quad_add( _vx_spline_( t ), _vy_spline_( t ),
					_vz_spline_( t ) );
			vr = IceBRG::dot_product( _vx_spline_( t ), _vy_spline_( t ),
					_vz_spline_( t ), _x_spline_( t ), _y_spline_( t ),
					_z_spline_( t ) ) / r;
			if ( std::fabs( vr ) >= std::fabs( v ) )
				vt = 0;
			else
				vt = IceBRG::quad_sub( v, vr );
			step_length_factor = _step_length_factor(v, r);

			// Update host for this timestep if necessary
			if ( _evolving_host_ )
			{
				num_host_parameters = ssize(_host_parameter_splines_);
				host_parameters.resize( num_host_parameters );
				for ( int i = 0; i < num_host_parameters; i++ )
				{
					host_parameters.at( i ) = _host_parameter_splines_.at( i )(
							t );
				}
				_current_host_ptr_->set_parameters( host_parameters );
			}

			// Calculate effects of tidal stripping and shocks
			// Effects of stripping
			m_ret = _tidal_strip_retained( _current_host_ptr_,
					_current_satellite_ptr_, r, vr, vt,
					t_step * step_length_factor,
					_sum_delta_rho_list_.at( counter - 1 ) );
			_m_ret_list_.push_back( _m_ret_list_.at( counter - 1 ) * m_ret );

			// Calculate adjusted fraction (so numerical errors don't cause a consistent
			// offset to add up).

			double adjusted_m_ret = _m_ret_list_.back() * init_mtot/safe_d(_current_satellite_ptr_->mtot());
			if( isbad(adjusted_m_ret) or (adjusted_m_ret>1.1) or (adjusted_m_ret<0)) adjusted_m_ret = 0;

			_current_satellite_ptr_->truncate_to_fraction( adjusted_m_ret );

			_m_vir_ret_list_.push_back( _current_satellite_ptr_->mvir() /
					init_mvir);

			if(m_ret == 0)
			{
				_likely_disrupted_ = true;
				_bad_result_ = true;
				break;
			}

			// Effects of shocking

			double t_shock = r / safe_d(v);
			double t_recover = _current_satellite_ptr_->othmtot()/(2*pi);
			double x;
			double gabdt_scaling_factor;
			x = max(0.,t_shock / safe_d(t_recover) / ( 2 * pi ));

			if(isbad(t_recover) or (t_recover<=0) or isbad(x) or (x<=0))
			{
				_likely_disrupted_ = true;
				_bad_result_ = true;
				break;
			}
			else
			{
				gabdt_scaling_factor = max(0.,1
									- ( t_step * step_length_factor
											/ safe_d( t_recover * _tidal_shocking_persistance_ ) ) );
			}

			current_gabdt.set_host_ptr( _current_host_ptr_ );
			current_gabdt.set_pos( _x_spline_( t ), _y_spline_( t ),
					_z_spline_( t ) );
			current_gabdt.set_dt( t_step * step_length_factor );
			current_gabdt.calc_dv();
			_gabdt_list_.push_back( current_gabdt );
			_sum_gabdt_list_.push_back(
					_sum_gabdt_list_.back() * gabdt_scaling_factor
							+ current_gabdt );

			current_deltarho = _delta_rho( counter, x,
                                    t_step * step_length_factor );

			_delta_rho_list_.push_back( current_deltarho );

			_sum_delta_rho_list_.push_back(
							_sum_delta_rho_list_.back() + current_deltarho);

			_rt_list_.push_back(
					_get_rt( _current_host_ptr_,
							_current_satellite_ptr_, r, vr, vt,
							t_step * step_length_factor,
							_sum_delta_rho_list_.at( counter - 1 ) ) );
			_rt_ratio_list_.push_back( _rt_list_.at( counter ) / _rvir( 0 ) );

			satellite_parameters = _current_satellite_ptr_->get_parameters(  );
			_satellite_parameter_data_.push_back( satellite_parameters );

			host_parameters =_current_host_ptr_->get_parameters(  );
			_host_parameter_data_.push_back( host_parameters );

		}
		catch ( const std::exception & e )
		{
		    const char* msg = e.what();
			_calculated_ = true; // We know there's an issue, so we won't calculate again unless something changes
			_bad_result_ = true;
			throw;
		}
	}
	if ( _record_full_data_ )
		_phase_output_list_.push_back(
				phase( _x_spline_( t ), _y_spline_( t ), _z_spline_( t ),
						_vx_spline_( t ), _vy_spline_( t ), _vz_spline_( t ),
						t ) );
	_calculated_ = true;

	if((_bad_result_) || (_likely_disrupted_))
		throw std::runtime_error("Error in calculating stripping for orbit segment.\n");
}

void IceBRG::stripping_orbit_segment::set_satellite_output_parameters(
		const std::vector< bool > & new_satellite_output_parameters )
{
	assert(_init_satellite_ptr_->num_parameters()==ssize(new_satellite_output_parameters));

	_satellite_output_parameters_ = new_satellite_output_parameters;
}
void IceBRG::stripping_orbit_segment::set_satellite_parameter_unitconvs(
		const std::vector< double > &new_satellite_parameter_unitconvs )
{
	assert(_init_satellite_ptr_->num_parameters()==ssize(new_satellite_parameter_unitconvs));

	_satellite_parameter_unitconvs_ = new_satellite_parameter_unitconvs;
}

void IceBRG::stripping_orbit_segment::set_host_output_parameters(
		const std::vector< bool > & new_host_output_parameters )
{
	assert(_init_host_ptr_->num_parameters()==ssize(new_host_output_parameters));

	_host_output_parameters_ = new_host_output_parameters;
}
void IceBRG::stripping_orbit_segment::set_host_parameter_unitconvs(
		const std::vector< double > & new_host_parameter_unitconvs )
{
	assert(_init_host_ptr_->num_parameters()==ssize(new_host_parameter_unitconvs));

	_host_parameter_unitconvs_ = new_host_parameter_unitconvs;
}

void IceBRG::stripping_orbit_segment::clear_satellite_output_parameters()
{
	_satellite_output_parameters_.clear();
}
void IceBRG::stripping_orbit_segment::clear_satellite_parameter_unitconvs()
{
	_satellite_parameter_unitconvs_.clear();
}

void IceBRG::stripping_orbit_segment::clear_host_output_parameters()
{
	_host_output_parameters_.clear();
}
void IceBRG::stripping_orbit_segment::clear_host_parameter_unitconvs()
{
	_host_parameter_unitconvs_.clear();
}

void IceBRG::stripping_orbit_segment::print_full_data(
		std::ostream *out, const bool include_header,
		const double m_ret_multiplier, const double m_vir_ret_multiplier, const bool silent ) const
{
	const ssize_t num_columns_base = 18;
	ssize_t num_extra_satellite_columns = 0;
	ssize_t num_extra_host_columns = 0;
	ssize_t num_rows = 0;
	std::vector< std::string > header;
	std::vector< std::vector< std::string > > data;
	std::stringstream ss;

	// Calculate data if necessary
	if ( ( !_calculated_ ) || ( !_record_full_data_ ) )
	{
		set_record_full_data( true );
		try
		{
			calc();
		}
		catch(const std::runtime_error &e)
		{
			std::cerr << "WARNING: Printing data for orbit which had a bad calculation.\n";
		}
	}

	num_rows = ssize(_phase_output_list_);
	if(_bad_result_ || _likely_disrupted_)
	{
		// Check how many rows we can actually output
		num_rows = min(num_rows,ssize(_m_ret_list_));
		num_rows = min(num_rows,ssize(_m_vir_ret_list_));
		num_rows = min(num_rows,ssize(_rt_list_));
		num_rows = min(num_rows,ssize(_rt_ratio_list_));
		num_rows = min(num_rows,ssize(_host_parameter_data_));
		num_rows = min(num_rows,ssize(_satellite_parameter_data_));
	}
	if ( num_rows < 2 )
	{
		throw std::runtime_error("No stripping data to output in stripping_orbit_segment::print_full_data.");
	}

	// Check that vectors for density_profile params are sane.
	// If so, set up extra columns. If not, no extra columns.
	num_extra_satellite_columns = 0;
	if ( !( ( ssize(_satellite_output_parameters_) == 0 )
			&& ( ssize(_satellite_parameter_unitconvs_) == 0 ) ) )
	{
		// If output_parameter_unitconvs was never initialized, set it to all 1
		if ( ssize(_satellite_parameter_unitconvs_) == 0 )
			_satellite_parameter_unitconvs_.resize(
				ssize(_satellite_output_parameters_), 1 );

		// Check both vectors are the same size
		if ( ssize(_satellite_output_parameters_)
				!= ssize(_satellite_parameter_unitconvs_) )
		{
			if ( !silent )
				std::cerr
						<< "WARNING: satellite_output_parameters and satellite_output_parameter_unitconvs have different sizes.\n"
						<< "No extra parameters will be output.\n";
		}
		else
		{
			// Count up extra columns to print out
			try
			{
				for ( ssize_t i = 0;
						i < ssize(_satellite_output_parameters_); i++ )
				{
					if ( _satellite_output_parameters_.at( i ) )
						num_extra_satellite_columns++;
				}
			}
			catch ( const std::exception &e )
			{
				if ( !silent )
					std::cerr
							<< "WARNING: Error reading satellite_output_parameters.\n"
							<< "No extra parameters will be output.\n"
							<< e.what();
				num_extra_satellite_columns = 0;
			}
		}
	}

	// As before, for host now
	num_extra_host_columns = 0;
	if ( !( ( ssize(_host_output_parameters_) == 0 )
			&& ( ssize(_host_parameter_unitconvs_) == 0 ) ) )
	{
		// If output_parameter_unitconvs was never initialized, set it to all 1
		if ( ssize(_host_parameter_unitconvs_) == 0 )
			_host_parameter_unitconvs_.resize( ssize(_host_output_parameters_),
					1 );

		// Check both vectors are the same size
		if ( ssize(_host_output_parameters_)
				!= ssize(_host_parameter_unitconvs_) )
		{
			if ( !silent )
				std::cerr
						<< "WARNING: host_output_parameters and host_output_parameter_unitconvs have different sizes.\n"
						<< "No extra parameters will be output.\n";
		}
		else
		{
			// Count up extra columns to print out
			try
			{
				for ( ssize_t i = 0; i < ssize(_host_output_parameters_);
						i++ )
				{
					if ( _host_output_parameters_.at( i ) )
						num_extra_host_columns++;
				}
			}
			catch ( const std::exception &e )
			{
				if ( !silent )
					std::cerr
							<< "WARNING: Error reading host_output_parameters.\n"
							<< "No extra parameters will be output.\n"
							<< e.what();
				num_extra_host_columns = 0;
			}
		}
	}

	int num_columns = num_columns_base + num_extra_satellite_columns
			+ num_extra_host_columns;

	header.resize( num_columns );
	make_vector_default( data, num_columns, num_rows );

	header[0] = "#";
	header[1] = "t";
	header[2] = "x";
	header[3] = "y";
	header[4] = "z";
	header[5] = "vx";
	header[6] = "vy";
	header[7] = "vz";
	header[8] = "d";
	header[9] = "v";
	header[10] = "m_ret";
	header[11] = "m_frac_lost";
	header[12] = "m_vir_ret";
	header[13] = "m_vir_frac_lost";
	header[14] = "comp_m_ret";
	header[15] = "comp_m_frac_lost";
	header[16] = "rt";
	header[17] = "rt/rvir";

	if ( num_extra_satellite_columns > 0 )
	{
		int extra_column_counter = 0;
		std::vector< std::string > parameter_names( 0 );

		parameter_names = _current_satellite_ptr_->get_parameter_names(  );

		for ( ssize_t i = 0; i < ssize(_satellite_output_parameters_);
				i++ )
		{
			if ( _satellite_output_parameters_.at( i ) )
			{
				try
				{
					ss.str( "" );
					ss << "Satellite_" << parameter_names.at( i );
					header[num_columns_base + extra_column_counter] = ss.str();
				}
				catch ( const std::out_of_range & )
				{
					ss.str( "" );
					ss << "Unknown_parameter_" << ( extra_column_counter + 1 );
					header[num_columns_base + extra_column_counter] = ss.str();
				}
				extra_column_counter++;
			}
		}
	}

	if ( num_extra_host_columns > 0 )
	{
		int extra_column_counter = 0;
		std::vector< std::string > parameter_names( 0 );

		parameter_names = _current_host_ptr_->get_parameter_names(  );

		for ( ssize_t i = 0; i < ssize(_host_output_parameters_); i++ )
		{
			if ( _host_output_parameters_.at( i ) )
			{
				try
				{
					ss.str( "" );
					ss << "Host_" << parameter_names.at( i );
					header[num_columns_base + num_extra_satellite_columns
							+ extra_column_counter] = ss.str();
				}
				catch ( const std::out_of_range & )
				{
					ss.str( "" );
					ss << "Unknown_parameter_"
							<< ( num_extra_satellite_columns
									+ extra_column_counter + 1 );
					header[num_columns_base + num_extra_satellite_columns
							+ extra_column_counter] = ss.str();
				}
				extra_column_counter++;
			}
		}
	}

	time_type t_min_to_use, t_max_to_use;

	if ( _override_t_min_ )
	{
		t_min_to_use = _t_min_override_val_;
	}
	else
	{
		t_min_to_use = _t_min_natural_value_;
	}

	if ( _override_t_max_ )
	{
		t_max_to_use = _t_max_override_val_;
	}
	else
	{
		t_max_to_use = _t_max_natural_value_;
	}

	if ( t_max_to_use <= t_min_to_use )
	{
		if(t_max_to_use == t_min_to_use)
			throw std::logic_error("t_max == t_min for stripping_orbit_segment.");
		std::cerr << "WARNING: t_max < t_min for stripping_orbit_segment; they'll be swapped for printing.\n";
		std::swap(t_max_to_use,t_min_to_use);
	}

	for ( ssize_t i = 0; i < num_rows; i++ )
	{
		data[0][i] = "";
		ss.str( "" );
		ss << _phase_output_list_.at( i ).t * unitconv::stoGyr;
		data[1][i] = ss.str();
		ss.str( "" );
		ss << _phase_output_list_.at( i ).x * unitconv::mtokpc;
		data[2][i] = ss.str();
		ss.str( "" );
		ss << _phase_output_list_.at( i ).y * unitconv::mtokpc;
		data[3][i] = ss.str();
		ss.str( "" );
		ss << _phase_output_list_.at( i ).z * unitconv::mtokpc;
		data[4][i] = ss.str();
		ss.str( "" );
		ss
				<< _phase_output_list_.at( i ).vx * unitconv::mtokpc
						/ unitconv::stoGyr;
		data[5][i] = ss.str();
		ss.str( "" );
		ss
				<< _phase_output_list_.at( i ).vy * unitconv::mtokpc
						/ unitconv::stoGyr;
		data[6][i] = ss.str();
		ss.str( "" );
		ss
				<< _phase_output_list_.at( i ).vz * unitconv::mtokpc
						/ unitconv::stoGyr;
		data[7][i] = ss.str();
		ss.str( "" );
		ss
				<< IceBRG::dist3d( _phase_output_list_.at( i ).x,
						_phase_output_list_.at( i ).y,
						_phase_output_list_.at( i ).z ) * unitconv::mtokpc;
		data[8][i] = ss.str();
		ss.str( "" );
		ss
				<< quad_add( _vx_spline_( _phase_output_list_.at( i ).t ),
						_vy_spline_( _phase_output_list_.at( i ).t ),
						_vz_spline_( _phase_output_list_.at( i ).t ) )
						* unitconv::mtokpc / unitconv::stoGyr;
		data[9][i] = ss.str();
		ss.str( "" );
		ss << _m_ret_list_.at( i ) * m_ret_multiplier;
		data[10][i] = ss.str();
		ss.str( "" );
		if ( i > 0 )
			ss << ( 1 - _m_ret_list_.at( i ) / _m_ret_list_.at( i - 1 ) )
							/ ( _phase_output_list_.at( i ).t
									- _phase_output_list_.at( i - 1 ).t );
		else
			ss << ( 1 - _m_ret_list_.at( i+1 ) / _m_ret_list_.at( i ) )
							/ ( _phase_output_list_.at( i+1 ).t
									- _phase_output_list_.at( i ).t );
		data[11][i] = ss.str();
		ss.str( "" );
		ss << _m_vir_ret_list_.at( i ) * m_vir_ret_multiplier;
		data[12][i] = ss.str();
		ss.str( "" );
		if ( i > 0 )
			ss << ( 1 - _m_vir_ret_list_.at( i ) / _m_vir_ret_list_.at( i - 1 ) )
							/ ( _phase_output_list_.at( i ).t
									- _phase_output_list_.at( i - 1 ).t );
		else
			ss << ( 1 - _m_vir_ret_list_.at( i+1 ) / _m_vir_ret_list_.at( i ) )
							/ ( _phase_output_list_.at( i+1 ).t
									- _phase_output_list_.at( i ).t );
		data[13][i] = ss.str();
		ss.str( "" );
		ss << _test_mass_spline_( _phase_output_list_.at( i ).t );
		data[14][i] = ss.str();
		ss.str( "" );
		if ( i > 0 )
			ss
					<< ( 1
							- _test_mass_spline_(
									_phase_output_list_.at( i ).t )
									/ safe_d(
											_test_mass_spline_(
													_phase_output_list_.at(
															i - 1 ).t ) ) )
							/ ( _phase_output_list_.at( i ).t
									- _phase_output_list_.at( i - 1 ).t );
		else
			ss << 0;
		data[15][i] = ss.str();
		ss.str( "" );
		ss << _rt_list_.at( i ) * unitconv::mtokpc;
		data[16][i] = ss.str();
		ss.str( "" );
		ss << _rt_ratio_list_.at( i );
		data[17][i] = ss.str();

		if ( num_extra_satellite_columns > 0 )
		{
			int extra_column_counter = 0;
			for ( ssize_t j = 0; j < ssize(_satellite_output_parameters_);
					j++ )
			{
				if ( _satellite_output_parameters_.at( j ) )
				{
					try
					{
						ss.str( "" );
						ss
								<< ( _satellite_parameter_data_.at( i ).at( j )
										/ safe_d(
												_satellite_parameter_unitconvs_.at(
														j ) ) );
						data[num_columns_base + extra_column_counter][i] =
								ss.str();
					}
					catch ( const std::out_of_range & )
					{
						ss.str( "" );
						ss << std::numeric_limits<double>::min();
						data[num_columns_base + extra_column_counter][i] =
								ss.str();
					}
					extra_column_counter++;
				}
			}
		}

		if ( num_extra_host_columns > 0 )
		{
			int extra_column_counter = 0;
			for ( ssize_t j = 0; j < ssize(_host_output_parameters_);
					j++ )
			{
				if ( _host_output_parameters_.at( j ) )
				{
					try
					{
						ss.str( "" );
						ss
								<< ( _host_parameter_data_.at( i ).at( j )
										/ safe_d(
												_host_parameter_unitconvs_.at(
														j ) ) );
						data[num_columns_base + num_extra_satellite_columns
								+ extra_column_counter][i] = ss.str();
					}
					catch ( const std::out_of_range & )
					{
						ss.str( "" );
						ss << std::numeric_limits<double>::min();
						data[num_columns_base + num_extra_satellite_columns
								+ extra_column_counter][i] = ss.str();
					}
					extra_column_counter++;
				}
			}
		}
	}

	if ( include_header )
		print_table( *out, data, header, silent );
	else
		print_table( *out, data, std::vector<std::string>(), silent );
}

flt_t IceBRG::stripping_orbit_segment::_delta_rho(
		const int index, const double x, time_type const & t_step,
		const bool silent ) const
{
	if ( index < 1 )
		return 0;
	double result = ( 1. / 6. / Gc * std::pow( 1 + x, _tidal_shocking_power_ )
			* _tidal_shocking_amplification_
			* ( 2 * ( _gabdt_list_[index] * _sum_gabdt_list_[index - 1] )
					+ _gabdt_list_[index] * _gabdt_list_[index] ) );
	if ( isbad( result ) )
	{
		if(!silent)
		{
			std::cerr << "ERROR: Could not calculate delta_rho in stripping_orbit_segment::_delta_rho.\n"
					<< "x = " << x << std::endl
					<< "_tidal_shocking_power_ = " << _tidal_shocking_power_ << std::endl
					<< "_tidal_shocking_amplification_ = " << _tidal_shocking_amplification_ << std::endl
					<< "gabdt product with sum = " << (_gabdt_list_[index] * _sum_gabdt_list_[index - 1]) << std::endl
					<< "gabdt product with self = " << (_gabdt_list_[index] * _gabdt_list_[index]) << std::endl;


			std::cerr.flush();
		}
		throw std::runtime_error("ERROR: Could not calculate delta_rho in stripping_orbit_segment::_delta_rho.\n");
	}
	return max( result, 0. );
}

const double IceBRG::stripping_orbit_segment::_step_length_factor( velocity_type const &  v, distance_type const &  r ) const
{
	/* Method to determine the appropriate factor for the step length.
	 * It uses the host's virial velocity and radius so that parts of
	 * the orbit with either high velocity or low distance from the host
	 * will take small steps (and the inverse as well).
	 */
	velocity_type v_0;
	distance_type r_0;

	// Determine proper scaling from the virial velocities, using the
	// most appropriate value depending on what's been loaded.
	if(_current_host_in_use_)
	{
		v_0 = _current_host_ptr_->vvir();
		r_0 = _current_host_ptr_->rvir();
	}
	else
	{
		if(_host_loaded_)
		{
			v_0 = _init_host_ptr_->vvir();
			r_0 = _init_host_ptr_->rvir();
		}
		else
		{
			v_0 = _v_0_;
			r_0 = _r_0_;
		}
	}

	// Figure out what the step factor would be from each of v and r
	// alone.

	// Low-v -> large step factor (and the inverse), bounded by max and min
	double step_length_factor_v = bound( _step_factor_min_,
			std::pow( std::fabs(v_0 / safe_d( v )), _step_length_power_ ), _step_factor_max_);

	// Low-r -> small step factor (and the inverse), bounded by max and min
	double step_length_factor_r = bound( _step_factor_min_,
			std::pow( std::fabs(r / safe_d( r_0 )), _step_length_power_ ), _step_factor_max_);

	if(isbad(step_length_factor_r))
	{
		std::cerr << "WARNING: Bad step_length_factor_r: " << step_length_factor_r
				<< "r = " << r
				<< "r_0 = " << r_0 << std::endl;
		if(isgood(step_length_factor_v)) return step_length_factor_v;
	}

	if(isbad(step_length_factor_v))
	{
		std::cerr << "WARNING: Bad step_length_factor_v: " << step_length_factor_r
				<< "v = " << r
				<< "v_0 = " << r_0 << std::endl;
		if(isgood(step_length_factor_r))
			return step_length_factor_v;
		else
			return 1;
	}

	// Use the smaller of the two calculated factors
	return min(step_length_factor_v,step_length_factor_r);
}

distance_type IceBRG::stripping_orbit_segment::_rvir(
		const int index ) const
{
	return _current_satellite_ptr_->rvir();
}
void IceBRG::stripping_orbit_segment::_pass_interpolation_type() const
{
	IceBRG::interpolator::allowed_interpolation_type type_to_pass = IceBRG::interpolator::SPLINE;
	if(_interpolation_type_ == IceBRG::stripping_orbit::UNSET)
	{
		if(length() * std::sqrt(step_factor_min()) > _spline_resolution_) type_to_pass = IceBRG::interpolator::LINEAR;
		if(length() * step_factor_min() > _spline_resolution_) type_to_pass = IceBRG::interpolator::LOWER;
	}
	else
	{
		if(_interpolation_type_ == IceBRG::stripping_orbit::LINEAR)
			type_to_pass = IceBRG::interpolator::LINEAR;
		if(_interpolation_type_ == IceBRG::stripping_orbit::SPLINE)
			type_to_pass = IceBRG::interpolator::SPLINE;
		if(_interpolation_type_ == IceBRG::stripping_orbit::UPPER)
			type_to_pass = IceBRG::interpolator::UPPER;
		if(_interpolation_type_ == IceBRG::stripping_orbit::LOWER)
			type_to_pass = IceBRG::interpolator::LOWER;
	}

	_x_spline_.set_interpolation_type(type_to_pass);
	_y_spline_.set_interpolation_type(type_to_pass);
	_z_spline_.set_interpolation_type(type_to_pass);
	_vx_spline_.set_interpolation_type(type_to_pass);
	_vy_spline_.set_interpolation_type(type_to_pass);
	_vz_spline_.set_interpolation_type(type_to_pass);
	for(ssize_t i=0; i<ssize(_host_parameter_splines_); i++)
	{
		_host_parameter_splines_[i].set_interpolation_type(type_to_pass);
	}
}

int IceBRG::stripping_orbit_segment::get_final_m_ret(
mass_type & m_ret, const bool silent ) const
{
	try
	{
		m_ret = final_m_ret();
	}
	catch(const std::runtime_error &e)
	{
		return 1;
	}
	return 0;
}
int IceBRG::stripping_orbit_segment::get_final_frac_m_ret( double & frac_m_ret,
		const bool silent ) const
{
	try
	{
		frac_m_ret = final_frac_m_ret();
	}
	catch(const std::runtime_error &e)
	{
		return 1;
	}
	return 0;
}
int IceBRG::stripping_orbit_segment::get_m_ret_points(
		std::vector< std::pair<double,double> > & m_ret_points_output,
		const bool silent ) const
{
	try
	{
		m_ret_points_output = m_ret_points();
	}
	catch(const std::runtime_error &e)
	{
		return 1;
	}
	return 0;
}

int IceBRG::stripping_orbit_segment::get_final_m_vir_ret(
mass_type & m_ret, const bool silent ) const
{
	try
	{
		m_ret = final_m_vir_ret();
	}
	catch(const std::runtime_error &e)
	{
		return 1;
	}
	return 0;
}
int IceBRG::stripping_orbit_segment::get_final_frac_m_vir_ret( double & frac_m_ret,
		const bool silent ) const
{
	try
	{
		frac_m_ret = final_frac_m_vir_ret();
	}
	catch(const std::runtime_error &e)
	{
		return 1;
	}
	return 0;
}
int IceBRG::stripping_orbit_segment::get_m_vir_ret_points(
		std::vector< std::pair<double,double> > & m_ret_points,
		const bool silent ) const
{
	try
	{
		m_ret_points = m_vir_ret_points();
	}
	catch(const std::runtime_error &e)
	{
		return 1;
	}
	return 0;
}

int IceBRG::stripping_orbit_segment::get_final_sum_deltarho(
long double & sum_deltarho, const bool silent ) const
{
	try
	{
		sum_deltarho = long_final_sum_deltarho();
	}
	catch(const std::runtime_error &e)
	{
		return 1;
	}
	return 0;
}

int IceBRG::stripping_orbit_segment::get_final_sum_deltarho(
		flt_t & sum_deltarho, const bool silent ) const
{
	try
	{
		sum_deltarho = final_sum_deltarho();
	}
	catch(const std::runtime_error &e)
	{
		return 1;
	}
	return 0;
}

int IceBRG::stripping_orbit_segment::get_final_sum_gabdt(
		gabdt & sum_gabdt, const bool silent ) const
{
	try
	{
		sum_gabdt = final_sum_gabdt();
	}
	catch(const std::runtime_error &e)
	{
		return 1;
	}
	return 0;
}

int IceBRG::stripping_orbit_segment::clone_final_satellite(
		density_profile * & final_satellite_clone, const bool silent ) const
{
	if(!_satellite_loaded_)
	{
		final_satellite_clone = new tNFW_profile();
		return 1;
	}
	try
	{
		final_satellite_clone = final_satellite()->density_profile_clone();
	}
	catch(const std::runtime_error &e)
	{
		final_satellite_clone = _init_satellite_ptr_->density_profile_clone();
		return 1;
	}
	return 0;
}

int IceBRG::stripping_orbit_segment::clone_final_host(
		density_profile * & final_host_clone, const bool silent ) const
{
	if(!_host_loaded_)
	{
		final_host_clone = new tNFW_profile();
		return 1;
	}
	try
	{
		final_host_clone = final_host()->density_profile_clone();
	}
	catch(const std::runtime_error &e)
	{
		final_host_clone = _init_host_ptr_->density_profile_clone();
		return 1;
	}
	return 0;
}

#if (1) // Get final data (throws exception on error)

mass_type IceBRG::stripping_orbit_segment::final_m_ret() const
{
	if(!_calculated_) calc();

	return _init_satellite_ptr_->mtot()*_m_ret_list_.back();
}
double IceBRG::stripping_orbit_segment::final_frac_m_ret() const
{
	if(!_calculated_) calc();

	return _m_ret_list_.back();
}
std::vector< std::pair<double,double> > IceBRG::stripping_orbit_segment::m_ret_points() const
{
	if(!_calculated_) calc();

	std::vector< std::pair<double,double> > result;
	for(ssize_t i = 0; i<ssize(_phase_output_list_); i++)
	{
		result.push_back( std::make_pair(
				_phase_output_list_[i].t,
				_m_ret_list_.at(i)));
	}

	return result;
}
mass_type IceBRG::stripping_orbit_segment::final_m_vir_ret() const
{
	if(!_calculated_) calc();

	return _init_satellite_ptr_->mvir()*_m_vir_ret_list_.back();
}
double IceBRG::stripping_orbit_segment::final_frac_m_vir_ret() const
{
	if(!_calculated_) calc();

	return _m_vir_ret_list_.back();
}
std::vector< std::pair<double,double> > IceBRG::stripping_orbit_segment::m_vir_ret_points() const
{
	if(!_calculated_) calc();

	std::vector< std::pair<double,double> > result;
	for(ssize_t i = 0; i<ssize(_phase_output_list_); i++)
	{
		result.push_back( std::make_pair(
				_phase_output_list_[i].t,
				_m_vir_ret_list_.at(i)));
	}

	return result;
}

flt_t IceBRG::stripping_orbit_segment::final_sum_deltarho() const
{
	if(!_calculated_) calc();

#ifdef _BRG_USE_UNITS_
	return IceBRG::unit_obj(_sum_delta_rho_list_.back(),-3,0,1,0,0);
#else
	return _sum_delta_rho_list_.back();
#endif
}

long double IceBRG::stripping_orbit_segment::long_final_sum_deltarho() const
{
	if(!_calculated_) calc();

	return _sum_delta_rho_list_.back();
}

IceBRG::gabdt IceBRG::stripping_orbit_segment::final_sum_gabdt() const
{
	if(!_calculated_) calc();

	return _sum_gabdt_list_.back();
}

bool IceBRG::stripping_orbit_segment::likely_disrupted() const
{
	if(!_calculated_)
	{
		try
		{
			calc();
		}
		catch(const std::exception &e)
		{
			// Do nothing on exception here
		}
	}
	return _likely_disrupted_;
}

const IceBRG::density_profile * IceBRG::stripping_orbit_segment::final_satellite() const
{
	if ( !_calculated_ )
	{
		if(_satellite_loaded_)
		{
			try
			{
				calc();
			} catch (const std::exception &e)
			{
			}
		}
		else
		{
			throw std::logic_error("Attempt to call stripping_orbit::final_satellite() without init_satellite assigned.\n");
		}
	}
	if( _bad_result_ || _likely_disrupted_ )
		return _init_satellite_ptr_; // Sanest option in this case
	else
		return _current_satellite_ptr_;
}

const IceBRG::density_profile * IceBRG::stripping_orbit_segment::final_host() const
{
	if ( !_calculated_ )
	{
		if(_host_loaded_)
		{
			try
			{
				calc();
			} catch (const std::exception &e)
			{
			}
		}
		else
		{
			throw std::logic_error("Attempt to call stripping_orbit::final_host() without init_host assigned.\n");
		}
	}
	if( _bad_result_ || _likely_disrupted_ )
		return _init_host_ptr_; // Sanest option in this case
	else
		return _current_host_ptr_;
}

#endif

#endif // end IceBRG::stripping_orbit_segment class function definitions
