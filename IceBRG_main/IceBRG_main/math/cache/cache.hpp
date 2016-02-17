/**********************************************************************\
  @file cache.hpp

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

#ifndef _BRG_CACHE_HPP_INCLUDED_
#define _BRG_CACHE_HPP_INCLUDED_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <sstream>
#include "IceBRG_main/common.h"

#include "IceBRG_main/Eigen.hpp"
#include "IceBRG_main/error_handling.h"
#include "IceBRG_main/file_access/ascii_table.hpp"
#include "IceBRG_main/file_access/open_file.hpp"
#include "IceBRG_main/file_access/trim_comments.hpp"
#include "IceBRG_main/math/misc_math.hpp"
#include "IceBRG_main/math/safe_math.hpp"

#include "IceBRG_main/units/units.hpp"

// Macro definitions

// SPCP: "Static Polymorphic Const Pointer"
// The constness isn't actually enforced, but this is for the reader's understanding
#define SPCP(name) static_cast<const name*>(this)

// SPP: "Static Polymorphic Pointer"
#define SPP(name) static_cast<name*>(this)

#define DECLARE_BRG_CACHE(class_name,name_base,Tin,Tout) \
class class_name : public IceBRG::brg_cache<class_name,Tin,Tout> \
{ \
private: \
 \
    static Tin _min_1_, _max_1_, _step_1_; \
	static IceBRG::ssize_t _resolution_1_; \
	static IceBRG::array_t<Tout> _results_; \
 \
	static IceBRG::str_t _file_name_; \
	static IceBRG::str_t _header_string_; \
 \
	static bool _loaded_, _initialised_; \
	static IceBRG::short_int_t _is_monotonic_; \
 \
	static int_t _sig_digits_; \
 \
	friend class IceBRG::brg_cache<class_name,Tin,Tout>; \
 \
protected: \
 \
	IceBRG::str_t _name_base() const \
	{ \
		return #name_base; \
	} \
 \
	Tout _calculate( Tin const & in_param ) const; \
 \
	void _load_cache_dependencies() const; \
}; \

#define DEFINE_BRG_CACHE(class_name,Tin,Tout, \
		init_min,init_max,init_step, \
		calc_method, \
		dependency_loading) \
	 \
	Tin class_name::_min_1_ = init_min; \
	Tin class_name::_max_1_ = init_max; \
	Tin class_name::_step_1_ = init_step; \
	IceBRG::ssize_t class_name::_resolution_1_ = 0; \
	 \
	bool class_name::_loaded_ = false; \
	bool class_name::_initialised_ = false; \
	 \
	IceBRG::short_int_t class_name::_is_monotonic_ = 0; \
	IceBRG::int_t class_name::_sig_digits_ = 8; \
	IceBRG::str_t class_name::_file_name_ = ""; \
	IceBRG::str_t class_name::_header_string_ = ""; \
	 \
	IceBRG::array_t<Tout> class_name::_results_; \
	 \
	Tout class_name::_calculate( Tin const & in_param ) const \
	{ \
		calc_method \
	} \
	void class_name::_load_cache_dependencies() const \
	{ \
		dependency_loading \
	}

namespace IceBRG
{

template<typename name, typename Tin=flt_t, typename Tout=flt_t>
class brg_cache
{
private:

	// Private variables
#if (1)

    static Tin _min_1_, _max_1_, _step_1_;
	static IceBRG::ssize_t _resolution_1_;
	static IceBRG::array_t<Tout> _results_;

	static IceBRG::str_t _file_name_;
	static IceBRG::str_t _header_string_;

	static bool _loaded_, _initialised_;
	static IceBRG::short_int_t _is_monotonic_;

	static int_t _sig_digits_;

#endif // Private variables

	// Private methods
#if (1)
	void _init() const
	{
		// We check for initialisation twice due to the critical section here.
		// It's expensive to enter, and we don't want to do anything inside it more than once,
		// so we check whether we need to both once outside it and once inside it.
		if(SPCP(name)->_initialised_) return;

		#ifdef _OPENMP
		#pragma omp critical(init_brg_cache)
		#endif
		if(!SPCP(name)->_initialised_)
		{
			SPCP(name)->_resolution_1_ = (ssize_t) max( ( ( SPCP(name)->_max_1_ - SPCP(name)->_min_1_ ) / safe_d(SPCP(name)->_step_1_)) + 1, 2);
			SPCP(name)->_file_name_ = SPCP(name)->_name_base() + "_cache.dat";
			SPCP(name)->_header_string_ = "# " + SPCP(name)->_name_base() + "_cache v1.0";

			SPCP(name)->_initialised_ = true;
		}
	}
	bool _critical_load() const
	{
		bool bad_result = false;
		#ifdef _OPENMP
		#pragma omp critical(load_brg_cache)
		#endif
		{
			try
			{
				SPCP(name)->_load();
			}
			catch(const std::exception &e)
			{
				handle_error_message(e.what());
				bad_result = true;
			}
		}
		return bad_result;
	}
	void _load() const
	{
		if ( SPCP(name)->_loaded_ ) return;

		std::ifstream in_file;
		std::string file_data;
		bool need_to_calc = false;
		ssize_t i;
		int_t loop_counter = 0;

		do
		{
			if ( loop_counter >= 2 )
			{
				throw std::runtime_error("Infinite loop detected trying to load " + SPCP(name)->_file_name_
						+ " in IceBRG::brg_cache.\n");
			}
			else
			{
				loop_counter++;
			}
			need_to_calc = false;

			try
			{
				open_file_input( in_file, SPCP(name)->_file_name_ );
			}
			catch(const std::exception &e)
			{
				need_to_calc = true;
				SPCP(name)->_calc();
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

			// Check that it has the right header
			getline( in_file, file_data );
			if ( file_data.compare( SPCP(name)->_header_string_ ) )
			{
				need_to_calc = true;
				SPCP(name)->_calc();
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

			// Trim out any other commented lines
			trim_comments_all_at_top( in_file );

			// Load range parameters;
			flt_t v_min, v_max, v_step;
			if ( !( in_file >> v_min >> v_max >> v_step ) )
			{
				need_to_calc = true;
				SPCP(name)->_calc();
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}
			SPCP(name)->_min_1_ = units_cast<Tin>(v_min);
			SPCP(name)->_max_1_ = units_cast<Tin>(v_max);
			SPCP(name)->_step_1_ = units_cast<Tin>(v_step);

			// Set up data
			SPCP(name)->_resolution_1_ = (ssize_t) max( ( ( SPCP(name)->_max_1_ - SPCP(name)->_min_1_ ) /
					safe_d(SPCP(name)->_step_1_)) + 1, 2);
			SPCP(name)->_results_.resize(SPCP(name)->_resolution_1_ );

			// Read in data

			flt_t temp_data;
			flt_t last_data=0;

			i = 0;
			SPCP(name)->_is_monotonic_ = 0;
			while ( ( !in_file.eof() ) && ( i < SPCP(name)->_resolution_1_ ) )
			{
				in_file >> temp_data;
				SPCP(name)->_results_[i] = units_cast<Tout>(temp_data);
				if(i==1)
				{
					// First monotonic check, so we don't compare to its past values
					if(temp_data > last_data)
					{
						SPCP(name)->_is_monotonic_ = 1;
					}
					else if(temp_data < last_data)
					{
						SPCP(name)->_is_monotonic_ = -1;
					}
					else
					{
						SPCP(name)->_is_monotonic_ = 0;
					}
				}
				else if(i>1)
				{
					// Check for monotonic increase/decrease
					if(temp_data > last_data)
					{
						if(SPCP(name)->_is_monotonic_ != 1)
							SPCP(name)->_is_monotonic_ = 0;
					}
					else if(temp_data < last_data)
					{
						if(SPCP(name)->_is_monotonic_ != -1)
							SPCP(name)->_is_monotonic_ = 0;
					}
					else
					{
						SPCP(name)->_is_monotonic_ = 0;
					}
				}
				last_data = temp_data;
				i++;
			}

			// Check that it was all read properly
			if ( i < SPCP(name)->_resolution_1_ )
			{
				need_to_calc = true;
				SPCP(name)->_calc();
				SPCP(name)->_unload();
				continue;
			}

		} while ( need_to_calc );

		// Finish up
		in_file.close();
		in_file.clear();
		SPCP(name)->_loaded_ = true;
	}
	void _unload() const
	{
		SPCP(name)->_loaded_ = false;
		set_zero(SPCP(name)->_results_);
	}
	void _calc() const
	{

		// Test that range is sane
		if ( ( SPCP(name)->_max_1_ <= SPCP(name)->_min_1_ ) || ( SPCP(name)->_step_1_ <= 0 ) )
		{
			throw std::runtime_error("ERROR: Bad range passed to brg_cache::_calc() for " +
					SPCP(name)->_name_base() + "\n");
		}

		// Print a message that we're generating the cache
		handle_notification("Generating " + SPCP(name)->_file_name_ + ". This may take some time.");

		// Set up data
		SPCP(name)->_resolution_1_ = (ssize_t) max( ( ( SPCP(name)->_max_1_ - SPCP(name)->_min_1_ ) /
				safe_d(SPCP(name)->_step_1_)) + 1, 2);
		SPCP(name)->_results_.resize(SPCP(name)->_resolution_1_ );

		// Calculate data
		bool bad_result = false;

		#ifdef _OPENMP
		#pragma omp parallel for schedule(dynamic)
		#endif
		for ( ssize_t i = 0; i < SPCP(name)->_resolution_1_; i++ )
		{
			Tout result = 0;
			Tin x = SPCP(name)->_min_1_ + i*SPCP(name)->_step_1_;
			try
			{
				result = SPCP(name)->_calculate(x);
			}
			catch(const std::exception &e)
			{
				handle_error_message(e.what());
				bad_result = true;
			}
			SPCP(name)->_results_[i] = result;
		}

		if(bad_result) throw std::runtime_error("One or more calculations failed in generating cache " +
				SPCP(name)->_name_base());
		SPCP(name)->_loaded_ = true;

		// Print a message that we've finished generating the cache
		handle_notification("Finished generating " + SPCP(name)->_file_name_ + "!");
	}
	void _output() const
	{

		std::ofstream out_file;
		std::string file_data;

		if ( !SPCP(name)->_loaded_ )
		{
			SPCP(name)->_calc();
		}

		open_file_output( out_file, SPCP(name)->_file_name_ );

		// Output header
		out_file << SPCP(name)->_header_string_ << "\n#\n";

		// Set number of significant digits
		out_file.precision(SPCP(name)->_sig_digits_);

		// Output range
		out_file << value_of(SPCP(name)->_min_1_) << "\t" << value_of(SPCP(name)->_max_1_) << "\t"
				<< value_of(SPCP(name)->_step_1_) << "\n";

		// Output data
		for ( ssize_t i = 0; i < SPCP(name)->_resolution_1_; i++ )
		{
			out_file << value_of(SPCP(name)->_results_[i]) << "\n";
		}

		out_file.close();
		out_file.clear();
	}
#endif // Private methods

protected:

	// Protected methods
	// These are made protected instead of private so base classes can overload them
#if (1)

	/// Long calculation function, which is used to generate the cache; must be overloaded by each
	/// child.
	Tout _calculate(Tin const & x) const;

	/// The default name (without extension) for the cache file; should be unique for each cache.
	str_t _name_base() const;

	/// This function should be overloaded to call each cache of the same dimensionality as
	/// this cache, which this depends upon in calculation. This is necessary in order to avoid critical
	/// sections of the same name being called recursively.
	void _load_cache_dependencies() const
	{
	}

#endif // Protected methods

public:

	// Public methods
#if (1)

	/**
	 * Set the name of the cache file to use.
	 *
	 * @param new_name The name of the cache file to use
	 */
	void set_file_name( str_t const & new_name )
	{
		SPP(name)->_file_name_ = new_name;
		if ( SPCP(name)->_loaded_ )
		{
			SPCP(name)->_unload();
		}
	} // set_file_name()

	/**
	 * Set the range of the independent parameter for wish you want values to be
	 * cached.
	 *
	 * @param new_min The new minimum value
	 * @param new_max The new maximum value
	 * @param new_step The number of points at which to cache the results
	 */
	void set_range( Tin const & new_min, Tin const & new_max,
			Tin const & new_step)
	{
		// First we try to load, so we can see if there are any changes from
		// the existing cache
		if ( !SPCP(name)->_loaded_ )
			SPCP(name)->_critical_load();

		// Go through variables, check if any are actually changed. If so, recalculate cache
		if ( ( SPCP(name)->_min_1_ != new_min ) || ( SPCP(name)->_max_1_ != new_max )
				|| ( SPCP(name)->_step_1_ != new_step ) )
		{
			SPP(name)->_min_1_ = new_min;
			SPP(name)->_max_1_ = new_max;
			SPP(name)->_step_1_ = new_step;

			SPCP(name)->_unload();
			SPCP(name)->_calc();
		}
	} // void set_range()

	/**
	 * Set the precision you wish the values stored in the cache to have.
	 *
	 * @param new_precision The desired precision
	 */
	void set_precision( ssize_t const & new_precision)
	{
		if ( new_precision > 0 )
		{
			SPP(name)->_sig_digits_ = min( new_precision, DBL_MAX_PRECISION );
		}
		else
		{
			throw std::runtime_error("Precision for dfa_cache must be > 0.\n");
		}
	} // void set_precision()

	/**
	 * Print the cached input and output values to an output stream.
	 *
	 * @param out The output stream you wish to print the cached values to.
	 */
	template<typename otype>
	void print( otype & out ) const
	{
		// Load if necessary
		if ( !SPCP(name)->_loaded_ )
		{
			SPCP(name)->_critical_load(SPCP(name)->_min_1_);
		}

		// Fill up header
		vector_t< str_t > header(3);
		header[0] = "#";
		header[1] = "x_1";
		header[2] = "y";

		// Fill up data
		vector_t< vector_t<str_t> > data(3);
		std::stringstream ss;
		for(ssize_t i_1=0; i_1<SPCP(name)->_resolution_1_; ++i_1)
		{
			data[0].push_back("");
			ss.str("");
			ss << value_of(SPCP(name)->_min_1_ + i_1*SPCP(name)->_step_1_);
			data[1].push_back(ss.str());
			ss.str("");
			ss << value_of(SPCP(name)->_results_[i_1]);
			data[2].push_back(ss.str());
		}

		print_table(out,data,header);
	}

	/**
	 * Get the result of the cached function for a given value.
	 *
	 * @param init_x The value for which you desired the cached result.
	 * @return The cached result for the input value.
	 */
	Tout get( const Tin & x ) const
	{
		Tin xlo, xhi;
		ssize_t x_i; // Lower nearby array point
		Tout result;

		// Load if necessary
		if ( !SPCP(name)->_loaded_ )
		{
			// Load any caches we depend upon before the critical section
			SPCP(name)->_load_cache_dependencies();

			if ( SPCP(name)->_critical_load() )
			{
				throw std::runtime_error("ERROR: Could neither load " + SPCP(name)->_file_name_ + " nor calculate in brg_cache::get()\n");
			}
		}

		x_i = (ssize_t)bound(0,
				( ( x - SPCP(name)->_min_1_ ) / SPCP(name)->_step_1_ ),
				SPCP(name)->_resolution_1_ - 2 );

		xlo = SPCP(name)->_min_1_ + SPCP(name)->_step_1_ * x_i;
		xhi = SPCP(name)->_min_1_ + SPCP(name)->_step_1_ * ( x_i + 1 );

		result = ( ( x - xlo ) * SPCP(name)->_results_[x_i + 1] + ( xhi - x ) * SPCP(name)->_results_[x_i] )
				/ SPCP(name)->_step_1_;

		return result;

	} // get()

	/**
	 * Get the independent parameter "x" from the corresponding dependent parameter "y". This only
	 * works if the function is monotonic increasing or monotonic decreasing.
	 *
	 * @param y The independent parameter of the cached function
	 * @return The corresponding independent parameter of the cached function
	 */
	Tin inverse_get( Tout const & y ) const
	{
		// Check if it's possible to do an inverse get
		if((SPCP(name)->_is_monotonic_!=1)&&((SPCP(name)->_is_monotonic_!=-1)))
		{
			// Not a monotonic function. Inverse get isn't possible
			str_t err = "ERROR: Attempt to use inverse_get in cache for " + SPCP(name)->_file_name_ +
					" for function which isn't monotonic.\n";
			throw std::runtime_error(err);
		}


		Tin xlo, xhi;
		Tout ylo, yhi;
		Tin result = 0;

		if ( !SPCP(name)->_loaded_ )
		{
			SPCP(name)->_critical_load();
		}
		if ( result == -1 )
		{
			str_t err = "ERROR: Could neither load " + SPCP(name)->_file_name_ +
					" nor calculate in brg_cache::inverse_get()\n";
			throw std::runtime_error(err);
		}

		if(SPCP(name)->_is_monotonic_==1)
		{

			for ( ssize_t x_i = 0; x_i < SPCP(name)->_resolution_1_ - 1; x_i++ )
			{
				// Loop through till we find the proper y or reach the end
				yhi = SPCP(name)->_results_[x_i];
				if ( ( yhi > y ) || (x_i >= SPCP(name)->_resolution_1_ - 2) )
				{
					ylo = SPCP(name)->_results_[x_i + 1];

					xlo = SPCP(name)->_min_1_ + SPCP(name)->_step_1_ * x_i;
					xhi = SPCP(name)->_min_1_ + SPCP(name)->_step_1_ * ( x_i + 1 );
					result = xlo + ( xhi - xlo ) * ( y - ylo ) / safe_d( yhi - ylo );
					break;
				}
			}
		} // if(_is_monotonic_==1)
		else
		{

			for ( ssize_t x_i = 0; x_i < SPCP(name)->_resolution_1_ - 1; x_i++ )
			{
				// Loop through till we find the proper y or reach the end
				ylo = SPCP(name)->_results_[x_i];
				if ( ( ylo < y ) || (x_i >= SPCP(name)->_resolution_1_ - 2) )
				{
					yhi = SPCP(name)->_results_[x_i + 1];

					xlo = SPCP(name)->_min_1_ + SPCP(name)->_step_1_ * x_i;
					xhi = SPCP(name)->_min_1_ + SPCP(name)->_step_1_ * ( x_i + 1 );
					result = xlo + ( xhi - xlo ) * ( y - ylo ) / safe_d( yhi - ylo );
					break;
				}
			}

		} // _is_monotonic == -1

		return result;

	}

	/// Load the cache, calculating if necessary
	void load() const
	{
		_load();
	}

	/// Reload the cache, calculating if necessary.
	void reload() const
	{
		_unload();
		_load();
	}

	/// Recalculate function. Call if you want to overwrite a cache when something's changed in the code
	/// (for instance, the _calculate() function has been altered)
	void recalc() const
	{
		SPCP(name)->_unload();
		SPCP(name)->_calc();
		SPCP(name)->_output();
	}

	// Constructor
	brg_cache()
	{
		SPP(name)->_init();
	}

	// Deconstructor
	virtual ~brg_cache()
	{
	}

#endif // Public methods

}; // class brg_cache

} // namespace IceBRG

// Undef macros
#undef SPP
#undef SPCP

#endif // __BRG_CACHE_HPP_INCLUDED__
