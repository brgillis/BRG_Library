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

#include <vector>

#include "brg/global.h"

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
	BRG_DISTANCE _R_mean_;
	BRG_MASS _m_mean_;
	double _z_mean_;
	double _mag_mean_;

	BRG_UNITS _delta_Sigma_t_mean_, _delta_Sigma_x_mean_;
	BRG_UNITS _delta_Sigma_t_mean_square_, _delta_Sigma_x_mean_square_;

#endif

public:

	// Constructors and destructor
#if(1)
	pair_bin_summary( CONST_BRG_DISTANCE_REF init_R_min=0, CONST_BRG_DISTANCE_REF init_R_max=0,
			CONST_BRG_MASS_REF init_m_min=0, CONST_BRG_MASS_REF init_m_max=0,
			double init_z_min=0, double init_z_max=0,
			double init_mag_min=0, double init_mag_max=0 );
	pair_bin_summary( const pair_bin & bin);
	virtual ~pair_bin_summary()
	{
	}
#endif

	// Adding and clearing data
#if(1)

	virtual void clear();

#endif

	// Count
	virtual size_t count() const
	{
		return _count_;
	}

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
		return delta_Sigma_t_std()/std::sqrt(_count_-1);
	}
	virtual BRG_UNITS delta_Sigma_x_stderr() const
	{
		if(_count_<2) return std::numeric_limits<double>::max();
		return delta_Sigma_x_std()/std::sqrt(_count_-1);
	}

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
