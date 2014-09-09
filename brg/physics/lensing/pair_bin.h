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

#include "brg/global.h"

#include "brg/physics/units/unit_obj.h"
#include "brg/vector/summary_functions.hpp"

namespace brgastro {

/**
 *
 */
class pair_bin {
private:

	// Bin limits
#if(1)

	BRG_DISTANCE _R_min_, _R_max_;
	BRG_MASS _m_min_, _m_max_;
	double _z_min_, _z_max_;
	double _mag_min_, _mag_max_;

#endif // Bin limits

	// Pair data
#if(1)

	std::vector<BRG_DISTANCE> _R_values_;
	std::vector<BRG_MASS> _m_values_;
	std::vector<double> _z_values_;
	std::vector<double> _mag_values_;
	std::vector<BRG_UNITS> _delta_Sigma_t_values_;
	std::vector<BRG_UNITS> _delta_Sigma_x_values_;

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

	// Limits accessors
#if(1)

	CONST_BRG_DISTANCE_REF R_min() const
	{
		return _R_min_;
	}
	CONST_BRG_DISTANCE_REF R_max() const
	{
		return _R_max_;
	}

	CONST_BRG_MASS_REF m_min() const
	{
		return _R_min_;
	}
	CONST_BRG_MASS_REF m_max() const
	{
		return _R_max_;
	}

	double z_min() const
	{
		return _z_min_;
	}
	double z_max() const
	{
		return _z_max_;
	}

	double mag_min() const
	{
		return _z_min_;
	}
	double mag_max() const
	{
		return _z_max_;
	}

#endif

	// Accessors to pair values
#if (1)

	std::vector<BRG_DISTANCE> R_values()
	{
		return _R_values_;
	}
	std::vector<BRG_MASS> m_values()
	{
		return _m_values_;
	}
	std::vector<double> z_values()
	{
		return _z_values_;
	}
	std::vector<double> mag_values()
	{
		return _mag_values_;
	}
	std::vector<BRG_UNITS> delta_Sigma_t_values()
	{
		return _delta_Sigma_t_values_;
	}
	std::vector<BRG_UNITS> delta_Sigma_x_values()
	{
		return _delta_Sigma_x_values_;
	}

#endif

	// Calculations on stored values
#if (1)

	BRG_UNITS delta_Sigma_t_mean();
	BRG_UNITS delta_Sigma_x_mean();

	BRG_UNITS delta_Sigma_t_std();
	BRG_UNITS delta_Sigma_x_std();

	BRG_UNITS delta_Sigma_t_stderr();
	BRG_UNITS delta_Sigma_x_stderr();

#endif
};

} // end namespace brgastro

#endif // _BRG_PAIR_BIN_H_INCLUDED_
