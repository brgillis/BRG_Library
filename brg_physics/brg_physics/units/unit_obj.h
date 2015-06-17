/**********************************************************************\
  @file unit_obj.h

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

// body file: brg_physics/units/unit_obj.cpp

/**********************************************************************\
 unit_obj.h
 -----------

 brgastro::unit_obj declarations: A class for a flt_type with units. This
 section is only included if the header directive _BRG_USE_UNITS_ is defined.
 This can be changed in the file global.h. If this is the case, the
 file units.cpp must be included and compiled with the project. This
 declares the following classes:

 brgastro::unit_obj - Double with generic units

 brgastro::unit_distance    - Child of unit_obj, with units fixed to
 distance
 brgastro::unit_time        - " fixed to time
 brgastro::unit_mass        - " fixed to mass
 brgastro::unit_temperature - " fixed to temp
 brgastro::unit_angle       - " fixed to angle
 brgastro::unit_charge      - " fixed to charge
 brgastro::unit_velocity    - " fixed to distance/time


 The class internally stores the value of the variable in default units
 and the powers of each unit type in a size-6 vector. All operations are
 performed in the default unit set, except for "set" and "get" functions
 which specify alternate units to use.

 I've found that using this class is most useful for debugging, as unit
 mismatches raise an obvious flag when something is wrong. As use of this
 class does slow the program down (order of ~50% longer for some programs
 I've run), it might be recommended to disable its use once debugging is
 complete, by toggling the _BRG_USE_UNITS_ flag in the brg_global.h file.

 As programmed, the class allows only basic units of distance, time,
 mass, temperature, angle, and charge. This can be expanded if necessary
 or convenient, but be aware of possible backwards-compatibility issues
 from changing the function formats.


 The class defaults to kms units for operations with unitless values. For
 instance, the expression:
 unit_distance d = 2;
 will assign d a value of 2 metres. The full default unit set is:

 Distance:    Metres    (m)
 Time:        Seconds   (s)
 Mass:        Kilograms (kg)
 Temperature: Kelvins   (K)
 Angle:       Radians   (rad)
 Charge:      Coulombs  (C)


 Functions for unit_objs:


 Key things to remember for these functions:
 -When in doubt, use default units (see list above), and expect
 unit_objs to behave like doubles.
 -For any unit conversions, use the "specified-to-default" format
 (unitconv::kmtom, not unitconv::mtokm). Do not raise the unit
 conversion factor to any power when you're giving it as a
 function argument, except for when using the
 set_value(new_val, conv_factor) function.
 -Values for angle and charge units/unit powers can always be left
 off if not needed; they'll take default values for all functions.


 brgastro::unit_obj x;

 Initializes x with a value of zero and unitless.


 brgastro::unit_obj x(const flt_type init_val=0,
 const flt_type d_units_power=0, const flt_type t_units_power=0,
 const flt_type m_units_power=0, const flt_type T_units_power=0,
 const flt_type a_units_power=0, const flt_type c_units_power=0);

 Initializes x with a value of init_val and specified unit powers. For
 instance, brgastro::unit_obj x(4,2) will initializes x with as 4 m^2.


 brgastro::unit_obj x(const flt_type init_val,
 const flt_type d_units, const float d_units_power,
 const flt_type t_units, const float t_units_power,
 const flt_type m_units, const float m_units_power,
 const flt_type T_units, const float T_units_power,
 const flt_type a_units=0, const float a_units_power=0,
 const flt_type c_units=0, const float c_units_power=0);

 Initializes x with a value of init_val in the specified units, with the
 specified unit powers. For instance,
 brgastro::unit_obj x( 4, unitconv::kmtom, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0 )
 will initialize x with a value of 4 km^2, which is internally stored as
 4e6 m^2. Always use the "kmtom"-like form of the unit conversion, with
 the default unit on the right, for this and all in all other class
 functions which use specified units.


 brgastro::unit_obj x(const unit_obj& other_unit_obj,
 const bool maintain_unit_fix=false);

 Initializes x as a copy of other_unit_obj. If maintain_unit_fix is set
 to true, and other_unit_obj is a value with fixed units (eg. it was
 declared as a unit_distance), x will have fixed units as well. Generally
 this is not necessary, and this value can be ignored for any explicit
 declarations.


 Set functions:


 const int_type brgastro::unit_obj::reset(const float d_units_power=0,
 const float t_units_power=0, const float m_units_power=0,
 const float T_units_power=0, const float a_units_power=0,
 const float c_units_power=0);

 Resets the value and optionally alters the unit powers, otherwise
 setting them to zero. Always returns a value of 0, but
 can be defined otherwise if necessary.


 const int_type brgastro::unit_obj::set_unit_powers(
 const float d_units_power=0, const float t_units_power=0,
 const float m_units_power=0, const float T_units_power=0,
 const float a_units_power=0, const float c_units_power=0);
 const int_type brgastro::unit_obj::set_unit_powers(
 std::vector<float> new_unitpowers);

 Set the unit powers of the variable, but DO NOT affect the stored value.
 Not recommended to be used, but available if necessary. Either function
 will return 1 if the variable's unit powers are fixed, and the latter
 will return 1 if the vector does not have length equal to the number
 of unit types.


 const int_type brgastro::unit_obj::set_value(const flt_type new_val,
 const flt_type conv_factor=1);

 Set the value of the variable in default units, or in the units
 specified by the conv_factor parameter. For instance,
 x.set_value(4, unitconv::kmtom/squarew(unitconv::hrtos)) will assign
 x a value of 4 km/hr^2, but it WILL NOT change the unit powers of
 x correspondingly. Use this only if you are sure x already has
 the correct unit powers. Always returns a value of 0.


 const int_type brgastro::unit_obj::set_value(const flt_type new_val,
 const flt_type d_units, const flt_type t_units, const flt_type m_units,
 const flt_type T_units, const flt_type a_units=1,
 const flt_type c_units=1);

 As above, but units are set individually. DO NOT raise units to
 the appropriate power; the function determines it automatically
 from the stored unit powers. For instance,
 x.set_value(4, unitconv::kmtom, unitconv::hrtos, 1, 1, 1, 1)
 will assign x a value of 4 km/hr if x has unit powers of distance/time,
 but it will assign x a value of 4 km/hr^2 if x has unit powers of
 distance/time^2.


 const int_type brgastro::unit_obj::set(const flt_type new_val,
 const std::vector<float> new_unitpowers);
 const int_type brgastro::unit_obj::set(const flt_type new_val,
 const float d_units_power, const float t_units_power,
 const float m_units_power, const float T_units_power,
 const float a_units_power=0, const float c_units_power=0);

 Sets the value in default units and changes the unit powers of the
 variable.


 const int_type brgastro::unit_obj::set(const flt_type new_val,
 const flt_type d_units, const float d_units_power,
 const flt_type t_units, const float t_units_power,
 const flt_type m_units, const float m_units_power,
 const flt_type T_units, const float T_units_power,
 const flt_type a_units=1, const float a_units_power=0,
 const flt_type c_units=1, const float c_units_power=0);

 Sets the value in the specified units and changes the unit powers of the
 variable. DO NOT raise units to the appropriate power; the function
 determines it automatically from the specified unit powers.


 Get functions:


 const flt_type brgastro::unit_obj::get_value() const;

 Returns the value in default units.


 const flt_type brgastro::unit_obj::get_value(const flt_type d_units,
 const flt_type t_units, const flt_type m_units, const flt_type T_units,
 const flt_type a_units=1, const flt_type C_units=1) const;

 Returns the value in specified units. DO NOT raise units to
 the appropriate power; the function determines it automatically
 from the stored unit powers.


 const flt_type brgastro::unit_obj::get_value(
 const flt_type d_units, const float d_units_power,
 const flt_type t_units, const float t_units_power,
 const flt_type m_units, const float m_units_power,
 const flt_type T_units, const float T_units_power,
 const flt_type a_units=1, const float a_units_power=0,
 const flt_type c_units=1, const float c_units_power=0) const;

 Returns the value in specified units, using specified unit powers
 instead of the stored powers. Typically should not be needed.


 const std::vector<float> brgastro::unit_obj::get_unit_powers() const;

 Returns a vector of the unit powers of the variable.


 const std::string brgastro::unit_obj::get_string() const;

 Returns a string writing out the variable's value in default units,
 along with the powers of those units.


 Specific get functions:
 These functions automatically adjust for the unit powers. For instance,
 if x has a value of 100 m^2, x.km() will return a value of 0.0001 (km^2).

 Get distance in units of...
 const flt_type brgastro::unit_obj::m() const;
 const flt_type brgastro::unit_obj::mm() const;
 const flt_type brgastro::unit_obj::um() const;
 const flt_type brgastro::unit_obj::nm() const;
 const flt_type brgastro::unit_obj::cm() const;
 const flt_type brgastro::unit_obj::angstrom() const;
 const flt_type brgastro::unit_obj::km() const;
 const flt_type brgastro::unit_obj::ltyr() const;
 const flt_type brgastro::unit_obj::pc() const;
 const flt_type brgastro::unit_obj::kpc() const;
 const flt_type brgastro::unit_obj::Mpc() const;
 const flt_type brgastro::unit_obj::mi() const;
 const flt_type brgastro::unit_obj::Mmi() const;
 const flt_type brgastro::unit_obj::ft() const;
 const flt_type brgastro::unit_obj::yd() const;

 Get time in units of...
 const flt_type brgastro::unit_obj::s() const;
 const flt_type brgastro::unit_obj::ms() const;
 const flt_type brgastro::unit_obj::cs() const;
 const flt_type brgastro::unit_obj::ns() const;
 const flt_type brgastro::unit_obj::us() const;
 const flt_type brgastro::unit_obj::min() const;
 const flt_type brgastro::unit_obj::hr() const;
 const flt_type brgastro::unit_obj::day() const;
 const flt_type brgastro::unit_obj::week() const;
 const flt_type brgastro::unit_obj::month() const;
 const flt_type brgastro::unit_obj::yr() const;
 const flt_type brgastro::unit_obj::kyr() const;
 const flt_type brgastro::unit_obj::Myr() const;
 const flt_type brgastro::unit_obj::Gyr() const;

 Get mass in units of...
 const flt_type brgastro::unit_obj::kg() const;
 const flt_type brgastro::unit_obj::gm() const;
 const flt_type brgastro::unit_obj::Mearth() const;
 const flt_type brgastro::unit_obj::Msun() const;
 const flt_type brgastro::unit_obj::ttMsun() const; // 10^10 Msun

 Get temperature in units of...
 Note: These automatically adjust for the different
 zero-values of the different temperature unit scales
 if the value is a pure temperature. If scaling only is
 preferred, simply use only the K() and degR() functions.
 const flt_type brgastro::unit_obj::K() const ;
 const flt_type brgastro::unit_obj::degC() const ;
 const flt_type brgastro::unit_obj::degF() const;
 const flt_type brgastro::unit_obj::degR() const;

 Get angle in units of...
 const flt_type brgastro::unit_obj::rad() const;
 const flt_type brgastro::unit_obj::deg() const;
 const flt_type brgastro::unit_obj::amin() const;
 const flt_type brgastro::unit_obj::asec() const;

 Get charge in units of...
 const flt_type brgastro::unit_obj::C() const;
 const flt_type brgastro::unit_obj::esu() const;

 Get other values in units of...
 Note: These functions do not automatically adjust for the variable's
 unit powers, since they involve more than one unit type.

 const flt_type brgastro::unit_obj::mps() const; // For velocity
 const flt_type brgastro::unit_obj::kmps() const; // ''
 const flt_type brgastro::unit_obj::c() const; // ''
 const flt_type brgastro::unit_obj::miphr() const; // ''
 const flt_type brgastro::unit_obj::kpcpGyr() const; // ''

 const flt_type brgastro::unit_obj::kmpspGyr() const; // For acceleration
 const flt_type brgastro::unit_obj::kpcpGyr2() const; // ''


 Other functions:

 const int_type round_powers();

 Rounds the unit powers. Precision is determined by the
 power_round_precision value. By default, it rounds so that no unit
 power will have a denominator greater than 12. Use this if round-off
 error is a concern after raising a variable to a fractional power.

 \**********************************************************************/

#ifndef _BRG_UNIT_OBJ_H_INCLUDED_
#define _BRG_UNIT_OBJ_H_INCLUDED_

#ifdef _BRG_USE_UNITS_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "brg/common.h"

#include "brg_physics/units/unit_conversions.hpp"

#define DIST_UNIT_INDEX 0
#define TIME_UNIT_INDEX 1
#define MASS_UNIT_INDEX 2
#define TEMP_UNIT_INDEX 3
#define ANGL_UNIT_INDEX 4
#define CHRG_UNIT_INDEX 5
#define NUM_UNIT_TYPES 6
#define ABS_ZERO_C -273.15
#define ABS_ZERO_F -459.67
// Power round precision - Set to zero to enforce no rounding. Higher is finer rounding
#define UNIT_POWER_ROUND_PRECISION_FACTOR 27720 // Lowest number that's evenly divisible by every integer 1 through 12

namespace brgastro
{
	class unit_obj;
}

namespace std
{

// Some specialisations to add to the std:: namespace
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const int_type rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const short_int_type rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const long int_type rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const int_type rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const int_type rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const long_int_type rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const flt_type rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const long_flt_type rhs);
brgastro::unit_obj pow(const brgastro::unit_obj &lhs, const float rhs);
brgastro::unit_obj sqrt(const brgastro::unit_obj &obj);
brgastro::unit_obj fabs(const brgastro::unit_obj &obj);

}

namespace brgastro
{

	class unit_obj
	{
	protected:

		// Default units: m, s, kg, K, rad, C
		std::vector<float> _unit_powers_;// Distance, Time, Mass, Temperature, Angle
		flt_type _value_;// Value in default units

		bool _fix_unit_powers_;// Used by derived classes to keep unit powers from being altered

	public:

		//------------------------------
		// Function prototypes (unit_obj)
		//------------------------------

		// Constructors
		unit_obj(flt_type init_val=0, float d_units_power=0, float t_units_power=0, float m_units_power=0,
				float T_units_power=0, float a_units_power=0, float c_units_power=0);
		unit_obj(flt_type init_val,
				flt_type d_units, float d_units_power, flt_type t_units, float t_units_power,
				flt_type m_units, float m_units_power, flt_type T_units, float T_units_power,
				 flt_type a_units=1, float a_units_power=0, flt_type c_units=1, float c_units_power=0);

		// Copy constructor
		unit_obj(const unit_obj& other_unit_obj, bool maintain_unit_fix=false);

		// Virtual Destructor (in case a derived class might need this functionality)
		virtual ~unit_obj();

		// Reset function
		void reset(float d_units_power=0, float t_units_power=0, float m_units_power=0,
				float T_units_power=0, float a_units_power=0, float c_units_power=0);

		// Set functions
		void set_unit_powers(float d_units_power=0, float t_units_power=0, float m_units_power=0,
				float T_units_power=0, float a_units_power=0, float c_units_power=0);// Use with care
		void set_unit_powers( std::vector<float> new_unitpowers);// Use with care
		void set_value(flt_type new_val, flt_type conv_factor=1);
		void set_value(flt_type new_val, flt_type d_units, flt_type t_units, flt_type m_units, flt_type T_units,
				flt_type a_units=1, flt_type c_units=1);
		void set(flt_type new_val, const std::vector<float> & new_unitpowers);
		void set(flt_type new_val, float d_units_power, float t_units_power, float m_units_power,
				float T_units_power, float a_units_power=0, float c_units_power=0);
		void set(flt_type new_val, flt_type d_units, float d_units_power, flt_type t_units, float t_units_power,
				flt_type m_units, float m_units_power, flt_type T_units, float T_units_power, flt_type a_units=1,
				float a_units_power=0, flt_type c_units=1, float c_units_power=0);

		// Get functions
		flt_type get_value() const;
		flt_type get_value(flt_type d_units, flt_type t_units, flt_type m_units, flt_type T_units, flt_type a_units=1,
				flt_type C_units=1) const;
		flt_type get_value(flt_type d_units, float d_units_power, flt_type t_units, float t_units_power, flt_type m_units,
				float m_units_power, flt_type T_units, float T_units_power, flt_type a_units=1,
				float a_units_power=0, flt_type c_units=1, float c_units_power=0) const;// Use with care
		std::vector<float> get_unit_powers() const;
		std::string get_string() const;

		// Misc functions
		void round_powers();

		// Get distance in units of...
		flt_type m() const;
		flt_type mm() const;
		flt_type um() const;
		flt_type nm() const;
		flt_type cm() const;
		flt_type angstrom() const;
		flt_type km() const;
		flt_type ltyr() const;
		flt_type pc() const;
		flt_type kpc() const;
		flt_type Mpc() const;
		flt_type mi() const;
		flt_type Mmi() const;
		flt_type ft() const;
		flt_type yd() const;

		// Get time in units of...
		flt_type s() const;
		flt_type ms() const;
		flt_type cs() const;
		flt_type ns() const;
		flt_type us() const;
		flt_type min() const;
		flt_type hr() const;
		flt_type day() const;
		flt_type week() const;
		flt_type month() const;
		flt_type yr() const;
		flt_type kyr() const;
		flt_type Myr() const;
		flt_type Gyr() const;

		// Get mass in units of...
		flt_type kg() const;
		flt_type gm() const;
		flt_type Mearth() const;
		flt_type Msun() const;
		flt_type ttMsun() const;

		// Get temperature in units of...
		flt_type K() const;
		flt_type degC() const;
		flt_type degF() const;
		flt_type degR() const;

		// Get angle in units of...
		flt_type rad() const;
		flt_type deg() const;
		flt_type amin() const;
		flt_type asec() const;

		// Get charge in units of...
		flt_type C() const;
		flt_type esu() const;

		// Get other values in units of...

		flt_type mps() const;// For velocity
		flt_type kmps() const;// ''
		flt_type c() const;// ''
		flt_type miphr() const;// ''
		flt_type kpcpGyr() const;// ''

		flt_type kmpspGyr() const;// For acceleration
		flt_type kpcpGyr2() const;// ''

		// Operator overloading
		unit_obj & operator=(const unit_obj & );
		unit_obj & operator=(int_type);
		unit_obj & operator=(short_int_type);
		unit_obj & operator=(long int_type);
		unit_obj & operator=(int_type);
		unit_obj & operator=(int_type);
		unit_obj & operator=(long_int_type);
		unit_obj & operator=(flt_type);
		unit_obj & operator=(long_flt_type);
		unit_obj & operator=(float);
		unit_obj operator+(const unit_obj & ) const;
		unit_obj operator+(int_type) const;
		unit_obj operator+(short_int_type) const;
		unit_obj operator+(long int_type) const;
		unit_obj operator+(int_type) const;
		unit_obj operator+(int_type) const;
		unit_obj operator+(long_int_type) const;
		unit_obj operator+(flt_type) const;
		unit_obj operator+(long_flt_type) const;
		unit_obj operator+(float) const;
		unit_obj operator-(const unit_obj & ) const;
		unit_obj operator-(int_type) const;
		unit_obj operator-(short_int_type) const;
		unit_obj operator-(long int_type) const;
		unit_obj operator-(int_type) const;
		unit_obj operator-(int_type) const;
		unit_obj operator-(long_int_type) const;
		unit_obj operator-(flt_type) const;
		unit_obj operator-(long_flt_type) const;
		unit_obj operator-(float) const;
		unit_obj operator*(const unit_obj & ) const;
		unit_obj operator*(int_type) const;
		unit_obj operator*(short_int_type) const;
		unit_obj operator*(long int_type) const;
		unit_obj operator*(int_type) const;
		unit_obj operator*(int_type) const;
		unit_obj operator*(long_int_type) const;
		unit_obj operator*(flt_type) const;
		unit_obj operator*(long_flt_type) const;
		unit_obj operator*(float) const;
		unit_obj operator/(const unit_obj & ) const;
		unit_obj operator/(int_type) const;
		unit_obj operator/(short_int_type) const;
		unit_obj operator/(long int_type) const;
		unit_obj operator/(int_type) const;
		unit_obj operator/(int_type) const;
		unit_obj operator/(long_int_type) const;
		unit_obj operator/(flt_type) const;
		unit_obj operator/(long_flt_type) const;
		unit_obj operator/(float) const;
		unit_obj & operator+=(const unit_obj & );
		unit_obj & operator+=(int_type);
		unit_obj & operator+=(short_int_type);
		unit_obj & operator+=(long int_type);
		unit_obj & operator+=(int_type);
		unit_obj & operator+=(int_type);
		unit_obj & operator+=(long_int_type);
		unit_obj & operator+=(flt_type);
		unit_obj & operator+=(long_flt_type);
		unit_obj & operator+=(float);
		unit_obj & operator-=(const unit_obj & );
		unit_obj & operator-=(int_type);
		unit_obj & operator-=(short_int_type);
		unit_obj & operator-=(long int_type);
		unit_obj & operator-=(int_type);
		unit_obj & operator-=(int_type);
		unit_obj & operator-=(long_int_type);
		unit_obj & operator-=(flt_type);
		unit_obj & operator-=(long_flt_type);
		unit_obj & operator-=(float);
		unit_obj & operator*=(const unit_obj & );
		unit_obj & operator*=(int_type);
		unit_obj & operator*=(short_int_type);
		unit_obj & operator*=(long int_type);
		unit_obj & operator*=(int_type);
		unit_obj & operator*=(int_type);
		unit_obj & operator*=(long_int_type);
		unit_obj & operator*=(flt_type);
		unit_obj & operator*=(long_flt_type);
		unit_obj & operator*=(float);
		unit_obj & operator/=(const unit_obj & );
		unit_obj & operator/=(int_type);
		unit_obj & operator/=(short_int_type);
		unit_obj & operator/=(long int_type);
		unit_obj & operator/=(int_type);
		unit_obj & operator/=(int_type);
		unit_obj & operator/=(long_int_type);
		unit_obj & operator/=(flt_type);
		unit_obj & operator/=(long_flt_type);
		unit_obj & operator/=(float);
		unit_obj & operator++();
		unit_obj operator++(int_type);
		unit_obj & operator--();
		unit_obj operator--(int_type);
		unit_obj operator-() const;
		bool operator<(const unit_obj & ) const;
		bool operator<(int_type) const;
		bool operator<(short_int_type) const;
		bool operator<(long int_type) const;
		bool operator<(int_type) const;
		bool operator<(int_type) const;
		bool operator<(long_int_type) const;
		bool operator<(flt_type) const;
		bool operator<(long_flt_type) const;
		bool operator<(float) const;
		bool operator>(const unit_obj & ) const;
		bool operator>(int_type) const;
		bool operator>(short_int_type) const;
		bool operator>(long int_type) const;
		bool operator>(int_type) const;
		bool operator>(int_type) const;
		bool operator>(long_int_type) const;
		bool operator>(flt_type) const;
		bool operator>(long_flt_type) const;
		bool operator>(float) const;
		bool operator==(const unit_obj & ) const;
		bool operator==(int_type) const;
		bool operator==(short_int_type) const;
		bool operator==(long int_type) const;
		bool operator==(int_type) const;
		bool operator==(int_type) const;
		bool operator==(long_int_type) const;
		bool operator==(flt_type) const;
		bool operator==(long_flt_type) const;
		bool operator==(float) const;
		bool operator<=(const unit_obj & ) const;
		bool operator<=(int_type) const;
		bool operator<=(short_int_type) const;
		bool operator<=(long int_type) const;
		bool operator<=(int_type) const;
		bool operator<=(int_type) const;
		bool operator<=(long_int_type) const;
		bool operator<=(flt_type) const;
		bool operator<=(long_flt_type) const;
		bool operator<=(float) const;
		bool operator>=(const unit_obj & ) const;
		bool operator>=(int_type) const;
		bool operator>=(short_int_type) const;
		bool operator>=(long int_type) const;
		bool operator>=(int_type) const;
		bool operator>=(int_type) const;
		bool operator>=(long_int_type) const;
		bool operator>=(flt_type) const;
		bool operator>=(long_flt_type) const;
		bool operator>=(float ) const;
		bool operator!=(const unit_obj & ) const;
		bool operator!=(int_type) const;
		bool operator!=(short_int_type) const;
		bool operator!=(long int_type) const;
		bool operator!=(int_type) const;
		bool operator!=(int_type) const;
		bool operator!=(long_int_type) const;
		bool operator!=(flt_type) const;
		bool operator!=(long_flt_type) const;
		bool operator!=(float) const;
		operator flt_type() const;

		// Friends

		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, int_type rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, short_int_type rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, long int_type rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, int_type rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, int_type rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, long_int_type rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, flt_type rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, long_flt_type rhs);
		friend brgastro::unit_obj std::pow(const brgastro::unit_obj &lhs, float rhs);
		friend brgastro::unit_obj std::sqrt(const brgastro::unit_obj &obj);
		friend brgastro::unit_obj std::fabs(const brgastro::unit_obj &obj);

	}; // unit_obj class

// Distance class
	class unit_distance: public unit_obj
	{

	public:
		unit_distance();
		unit_distance(flt_type init_val, flt_type conv_factor=1);
		unit_distance(const unit_obj &other_unit_obj);

	}; // class distance

	class unit_time: public unit_obj
	{

	public:
		unit_time();
		unit_time(flt_type init_val, flt_type conv_factor=1);
		unit_time(const unit_obj &other_unit_obj);

	}; // class time

	class unit_velocity: public unit_obj
	{
	public:
		unit_velocity();
		unit_velocity(flt_type init_val, flt_type conv_factor=1);
		unit_velocity(const unit_obj &other_unit_obj);

	}; // class unit_velocity

	class unit_mass: public unit_obj
	{
	public:
		unit_mass();
		unit_mass(flt_type init_val, flt_type conv_factor=1);
		unit_mass(const unit_obj &other_unit_obj);

	}; // class mass

	class unit_temperature: public unit_obj
	{
	public:
		unit_temperature();
		unit_temperature(flt_type init_val, flt_type conv_factor=1);
		unit_temperature(const unit_obj &other_unit_obj);

		flt_type K() const;
		flt_type C() const;
		flt_type R() const;
		flt_type F() const;

	}; // class unit_temperature

	class unit_angle: public unit_obj
	{

	public:
		unit_angle();
		unit_angle(flt_type init_val, flt_type conv_factor=1);
		unit_angle(const unit_obj &other_unit_obj);

	}; // class unit_angle

	class unit_charge: public unit_obj
	{

	public:
		unit_charge();
		unit_charge(flt_type init_val, flt_type conv_factor=1);
		unit_charge(const unit_obj &other_unit_obj);

	}; // class unit_charge

	// brgastro function overloads for unit_objs
#if (1)

	inline void set_zero( brgastro::unit_obj & v1)
	{
		v1 = unit_obj(0);
	}
	inline void set_zero( brgastro::unit_distance & v1)
	{
		v1 = unit_obj(0);
	}
	inline void set_zero( brgastro::unit_time & v1)
	{
		v1 = unit_obj(0);
	}
	inline void set_zero( brgastro::unit_velocity & v1)
	{
		v1 = unit_obj(0);
	}
	inline void set_zero( brgastro::unit_mass & v1)
	{
		v1 = unit_obj(0);
	}
	inline void set_zero( brgastro::unit_temperature & v1)
	{
		v1 = unit_obj(0);
	}
	inline void set_zero( brgastro::unit_angle & v1)
	{
		v1 = unit_obj(0);
	}
	inline void set_zero( brgastro::unit_charge & v1)
	{
		v1 = unit_obj(0);
	}

	inline bool isinf( unit_obj val )
	{
		return std::fabs( val.get_value() ) > std::numeric_limits<flt_type>::max();
	}

	inline brgastro::unit_obj square( const brgastro::unit_obj & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj square( const brgastro::unit_distance & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj square( const brgastro::unit_time & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj square( const brgastro::unit_velocity & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj square( const brgastro::unit_mass & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj square( const brgastro::unit_temperature & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj square( const brgastro::unit_angle & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj square( const brgastro::unit_charge & v1)
	{
		return v1*v1;
	}
	inline brgastro::unit_obj cube( const brgastro::unit_obj & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj cube( const brgastro::unit_distance & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj cube( const brgastro::unit_time & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj cube( const brgastro::unit_velocity & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj cube( const brgastro::unit_mass & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj cube( const brgastro::unit_temperature & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj cube( const brgastro::unit_angle & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj cube( const brgastro::unit_charge & v1)
	{
		return v1*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_obj & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_distance & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_time & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_velocity & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_mass & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_temperature & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_angle & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj quart( const brgastro::unit_charge & v1)
	{
		return square(v1)*square(v1);
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_obj & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_distance & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_time & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_velocity & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_mass & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_temperature & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_angle & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inverse( const brgastro::unit_charge & v1)
	{
		return 1/v1;
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_obj & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_distance & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_time & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_velocity & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_mass & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_temperature & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_angle & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_square( const brgastro::unit_charge & v1)
	{
		return inverse(square(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_obj & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_distance & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_time & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_velocity & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_mass & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_temperature & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_angle & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_cube( const brgastro::unit_charge & v1)
	{
		return inverse(cube(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_obj & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_distance & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_time & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_velocity & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_mass & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_temperature & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_angle & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj inv_quart( const brgastro::unit_charge & v1)
	{
		return inverse(quart(v1));
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_obj & v1, int_type p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_distance & v1, int_type p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_time & v1, int_type p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_velocity & v1, int_type p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_mass & v1, int_type p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_temperature & v1, int_type p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_angle & v1, int_type p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
	inline brgastro::unit_obj ipow( const brgastro::unit_charge & v1, int_type p)
	{
		return ipow((brgastro::unit_obj)v1,p);
	}
#endif

} // namespace brgastro

// Overloaded operators relating to unit_objs

brgastro::unit_obj operator+(int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(short_int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(long int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(long_int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(flt_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(long_flt_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator+(float lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(short_int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(long int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(long_int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(flt_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(long_flt_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator-(float lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(short_int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(long int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(long_int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(flt_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(long_flt_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator*(float lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(short_int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(long int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(long_int_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(flt_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(long_flt_type lhs, const brgastro::unit_obj &rhs);
brgastro::unit_obj operator/(float lhs, const brgastro::unit_obj &rhs);

bool operator<(int_type lhs, const brgastro::unit_obj &rhs);
bool operator<(short_int_type lhs, const brgastro::unit_obj &rhs);
bool operator<(long int_type lhs, const brgastro::unit_obj &rhs);
bool operator<(int_type lhs, const brgastro::unit_obj &rhs);
bool operator<(int_type lhs, const brgastro::unit_obj &rhs);
bool operator<(long_int_type lhs, const brgastro::unit_obj &rhs);
bool operator<(flt_type lhs, const brgastro::unit_obj &rhs);
bool operator<(long_flt_type lhs, const brgastro::unit_obj &rhs);
bool operator<(float lhs, const brgastro::unit_obj &rhs);
bool operator>(int_type lhs, const brgastro::unit_obj &rhs);
bool operator>(short_int_type lhs, const brgastro::unit_obj &rhs);
bool operator>(long int_type lhs, const brgastro::unit_obj &rhs);
bool operator>(int_type lhs, const brgastro::unit_obj &rhs);
bool operator>(int_type lhs, const brgastro::unit_obj &rhs);
bool operator>(long_int_type lhs, const brgastro::unit_obj &rhs);
bool operator>(flt_type lhs, const brgastro::unit_obj &rhs);
bool operator>(long_flt_type lhs, const brgastro::unit_obj &rhs);
bool operator>(float lhs, const brgastro::unit_obj &rhs);
bool operator<=(int_type lhs, const brgastro::unit_obj &rhs);
bool operator<=(short_int_type lhs, const brgastro::unit_obj &rhs);
bool operator<=(long int_type lhs, const brgastro::unit_obj &rhs);
bool operator<=(int_type lhs, const brgastro::unit_obj &rhs);
bool operator<=(int_type lhs, const brgastro::unit_obj &rhs);
bool operator<=(long_int_type lhs, const brgastro::unit_obj &rhs);
bool operator<=(flt_type lhs, const brgastro::unit_obj &rhs);
bool operator<=(long_flt_type lhs, const brgastro::unit_obj &rhs);
bool operator<=(float lhs, const brgastro::unit_obj &rhs);
bool operator>=(int_type lhs, const brgastro::unit_obj &rhs);
bool operator>=(short_int_type lhs, const brgastro::unit_obj &rhs);
bool operator>=(long int_type lhs, const brgastro::unit_obj &rhs);
bool operator>=(int_type lhs, const brgastro::unit_obj &rhs);
bool operator>=(int_type lhs, const brgastro::unit_obj &rhs);
bool operator>=(long_int_type lhs, const brgastro::unit_obj &rhs);
bool operator>=(flt_type lhs, const brgastro::unit_obj &rhs);
bool operator>=(long_flt_type lhs, const brgastro::unit_obj &rhs);
bool operator>=(float lhs, const brgastro::unit_obj &rhs);
bool operator==(int_type lhs, const brgastro::unit_obj &rhs);
bool operator==(short_int_type lhs, const brgastro::unit_obj &rhs);
bool operator==(long int_type lhs, const brgastro::unit_obj &rhs);
bool operator==(int_type lhs, const brgastro::unit_obj &rhs);
bool operator==(int_type lhs, const brgastro::unit_obj &rhs);
bool operator==(long_int_type lhs, const brgastro::unit_obj &rhs);
bool operator==(flt_type lhs, const brgastro::unit_obj &rhs);
bool operator==(long_flt_type lhs, const brgastro::unit_obj &rhs);
bool operator==(float lhs, const brgastro::unit_obj &rhs);
bool operator!=(int_type lhs, const brgastro::unit_obj &rhs);
bool operator!=(short_int_type lhs, const brgastro::unit_obj &rhs);
bool operator!=(long int_type lhs, const brgastro::unit_obj &rhs);
bool operator!=(int_type lhs, const brgastro::unit_obj &rhs);
bool operator!=(int_type lhs, const brgastro::unit_obj &rhs);
bool operator!=(long_int_type lhs, const brgastro::unit_obj &rhs);
bool operator!=(flt_type lhs, const brgastro::unit_obj &rhs);
bool operator!=(long_flt_type lhs, const brgastro::unit_obj &rhs);
bool operator!=(float lhs, const brgastro::unit_obj &rhs);

std::ostream & operator<<(std::ostream &out, brgastro::unit_obj &obj);

#endif

#endif // _BRG_UNIT_OBJ_H_INCLUDED_
