/*
 * cache_nd.hpp
 *
 *  Created on: 25 Apr 2014
 *      Author: brg
 */

#ifndef _BRG_CACHE_ND_HPP_INCLUDED_
#define _BRG_CACHE_ND_HPP_INCLUDED_

#ifndef BRG_CACHE_ND_NAME_SIZE
#define BRG_CACHE_ND_NAME_SIZE 9 // Needs an end character, so will only actually allow 8 chars
#endif

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <exception>

#include "brg/global.h"

#include "brg/file_functions.h"
#include "brg/physics/units/unit_obj.h"
#include "brg/vector/vector.hpp"
#include "brg/vector/elementwise_functions.hpp"

// Macro definitions

// SPCP: "Static Polymorphic Const Pointer"
// The constness isn't actually enforced, but this is for the reader's understanding
#define SPCP(name) static_cast<const name*>(this)

// SPP: "Static Polymorphic Pointer"
#define SPP(name) static_cast<name*>(this)

#define DECLARE_BRG_CACHE_ND_STATIC_VARS()		       \
	static brgastro::vector<double> _mins_, _maxes_, _steps_;      \
	static brgastro::vector<unsigned int> _resolutions_;           \
	static brgastro::vector<double> _results_;                     \
											                       \
	static std::string _file_name_;                                \
	static unsigned int _version_number_;                          \
										                           \
	static bool _loaded_, _initialised_;                           \
											                       \
	static unsigned short int _num_dim_;

// Be careful when using this not to use the default constructor for init_steps, which would result in
// divide-by-zero errors
#define DEFINE_BRG_CACHE_ND_STATIC_VARS(class_name,init_mins,init_maxes,init_steps,init_num_dim) \
	brgastro::vector<double> brgastro::class_name::_mins_ = init_mins;	                         \
	brgastro::vector<double> brgastro::class_name::_maxes_ = init_maxes;                         \
	brgastro::vector<double> brgastro::class_name::_steps_ = init_steps;                         \
	bool brgastro::class_name::_loaded_ = false;							                     \
	bool brgastro::class_name::_initialised_ = false;					                         \
	brgastro::vector<unsigned int> brgastro::class_name::_resolutions_ =                         \
		max( (((brgastro::class_name::_maxes_-brgastro::class_name::_mins_) /                    \
				safe_d(brgastro::class_name::_steps_))+1), 1);                                   \
	std::string brgastro::class_name::_file_name_ = "";					     	                 \
	unsigned int brgastro::class_name::_version_number_ = 1;		     		                 \
	brgastro::vector<double> brgastro::class_name::_results_;                                    \
	unsigned short int brgastro::class_name::_num_dim_ = init_num_dim;

namespace brgastro
{

template<typename name>
class brg_cache_nd
{
private:

	// Private variables
#if (1)

	DECLARE_BRG_CACHE_ND_STATIC_VARS();

#endif // Private variables

	// Private methods
#if (1)
	void _init() const throw()
	{
		// We check for initialisation twice due to the critical section here.
		// It's expensive to enter, and we don't want to do anything inside it more than once,
		// so we check whether we need to both once outside it and once inside it.
		if(SPCP(name)->_initialised_) return;

		#pragma omp critical(init_brg_cache_nd)
		if(!SPCP(name)->_initialised_)
		{
			SPCP(name)->_resolutions_ = max( (((SPCP(name)->_maxes_-SPCP(name)->_mins_) / safe_d(SPCP(name)->_steps_))+1), 1);
			SPCP(name)->_file_name_ = SPCP(name)->_name_base() + "_cache.bin";
			SPCP(name)->_version_number_ = 1; // This should be changed when there are changes to this code

			SPCP(name)->_initialised_ = true;
		}
	}
	void _load( const bool silent = false ) const
	{
		std::ifstream in_file;
		std::string file_data;
		bool need_to_calc = false;
		int loop_counter = 0;

		if ( SPCP(name)->_loaded_ )
			return;

		do
		{
			if ( loop_counter >= 2 )
			{
				throw std::runtime_error("Infinite loop detected trying to load " + SPCP(name)->_file_name_ + " in brgastro::brg_cache_nd.\n");
			}
			else
			{
				loop_counter++;
			}
			need_to_calc = false;

			if ( open_bin_file( in_file, SPCP(name)->_file_name_, true ) )
			{
				need_to_calc = true;
				SPCP(name)->_calc( silent );
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

			// Check that it has the right name and version

			char file_name[BRG_CACHE_ND_NAME_SIZE];
			unsigned int file_version = std::numeric_limits<unsigned short int>::max();

			in_file.read(file_name,BRG_CACHE_ND_NAME_SIZE);
			in_file.read((char *)&file_version,sizeof(file_version));

			if( (!in_file) || (((std::string)file_name) != SPCP(name)->_name_base()) ||
					(file_version != SPCP(name)->_version_number_) )
			{
				need_to_calc = true;
				SPCP(name)->_calc( silent );
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

			// Load range parameters;
			SPCP(name)->_mins_.resize(SPCP(name)->_num_dim_,0);
			SPCP(name)->_maxes_.resize(SPCP(name)->_num_dim_,0);
			SPCP(name)->_steps_.resize(SPCP(name)->_num_dim_,0);
			for(unsigned int i = 0; i < SPCP(name)->_num_dim_; i++)
			{
				in_file.read((char *)&(SPCP(name)->_mins_[i]),sizeof(SPCP(name)->_mins_[i]));
				in_file.read((char *)&(SPCP(name)->_maxes_[i]),sizeof(SPCP(name)->_maxes_[i]));
				in_file.read((char *)&(SPCP(name)->_steps_[i]),sizeof(SPCP(name)->_steps_[i]));
			}
			if ( !(in_file) )
			{
				need_to_calc = true;
				SPCP(name)->_calc( silent );
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

			// Set up data
			SPCP(name)->_resolutions_ = max( (((SPCP(name)->_maxes_-SPCP(name)->_mins_) / safe_d(SPCP(name)->_steps_))+1.), 1.);
			SPCP(name)->_results_.reshape(SPCP(name)->_resolutions_.v(),0);

			// Read in data

			// Initialize
			unsigned int i = 0;
			const std::streamsize size = sizeof(SPCP(name)->_results_[0]); // Store the size for speed
			brgastro::vector<unsigned int> position(SPCP(name)->_num_dim_,0);

			while ( ( !in_file.eof() ) && ( i < product(SPCP(name)->_resolutions_) ) && (in_file) )
			{
				in_file.read((char *)&(SPCP(name)->_results_(position)),size);
				i++;

				for(unsigned int d=0; d<SPCP(name)->_num_dim_; d++)
				{
					position[d]++;
					if(position[d] != SPCP(name)->_resolutions_[d])
						break;
					position[d] = 0;
					// If we get here, we'll go on to increase the next index by 1
				}
			}

			// Check that it was all read properly
			if ( (i != product(SPCP(name)->_resolutions_)) || (!in_file) )
			{
				need_to_calc = true;
				SPCP(name)->_calc( silent );
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

		} while ( need_to_calc );

		// Finish up
		in_file.close();
		in_file.clear();
		SPCP(name)->_loaded_ = true;
	}
	void _unload() const throw()
	{
		SPCP(name)->_loaded_ = false;
		SPCP(name)->_results_.clear();
	}
	void _calc( const bool silent = false ) const
	{
		// Test that range is sane
		for(unsigned int i = 0; i < SPCP(name)->_num_dim_; i++)
		{
			try {
				if ( ( SPCP(name)->_maxes_.at(i) <= SPCP(name)->_mins_.at(i) ) || ( SPCP(name)->_steps_.at(i) <= 0 ) )
				{
					throw std::runtime_error("ERROR: Bad range passed to brg_cache::_calc() for " + static_cast<const name*>(this)->_name_base() + "\n");
				}
			} catch (std::exception &e) {
				throw std::runtime_error("ERROR: Bad range passed to brg_cache::_calc() for " + static_cast<const name*>(this)->_name_base() + "\n");
			}
		}

		// Print a message that we're generating the cache if not silent
		if(!silent)
		{
			std::cout << "Generating " << SPCP(name)->_file_name_ << ". This may take some time.\n";
		}

		// Set up data
		SPCP(name)->_resolutions_.resize(SPCP(name)->_num_dim_);
		SPCP(name)->_resolutions_ = max( (((SPCP(name)->_maxes_-SPCP(name)->_mins_) / safe_d(SPCP(name)->_steps_))+1.), 1.);
		SPCP(name)->_results_.reshape(SPCP(name)->_resolutions_.v() );

		brgastro::vector<unsigned int> position(SPCP(name)->_num_dim_,0);
		brgastro::vector<double> x(SPCP(name)->_num_dim_,0);
		for ( unsigned int i = 0; i < SPCP(name)->_results_.size(); i++ )
		{
			x = SPCP(name)->_mins_ + SPCP(name)->_steps_*position;
			double result = 0;
			SPCP(name)->_results_(position) = SPCP(name)->_calculate(x);

			for(unsigned int d=0; d<SPCP(name)->_num_dim_; d++)
			{
				position[d]++;
				if(position[d] != SPCP(name)->_resolutions_[d])
					break;
				position[d] = 0;
				// If we get here, we'll go on to increase the next index by 1
			}
		}

		// Print a message that we've finished generating the cache if not silent.
		if(!silent)
		{
			std::cout << "Finished generating " << SPCP(name)->_file_name_ << "!\n";
		}

		SPCP(name)->_loaded_ = true;
	}
	void _output( const bool silent = false ) const
	{

		std::ofstream out_file;
		std::string file_data;

		if ( !SPCP(name)->_loaded_ )
		{
			SPCP(name)->_calc( silent );;
		}

		open_bin_file( out_file, SPCP(name)->_file_name_, true );

		// Output name and version

		std::string file_name = SPCP(name)->_name_base();
		unsigned int file_version = SPCP(name)->_version_number_;

		out_file.write(file_name.c_str(),BRG_CACHE_ND_NAME_SIZE);
		out_file.write((char *)&file_version,sizeof(file_version));

		// Output range parameters
		for(unsigned int i = 0; i < SPCP(name)->_num_dim_; i++)
		{
			out_file.write((char *)&(SPCP(name)->_mins_[i]),sizeof(SPCP(name)->_mins_[i]));
			out_file.write((char *)&(SPCP(name)->_maxes_[i]),sizeof(SPCP(name)->_maxes_[i]));
			out_file.write((char *)&(SPCP(name)->_steps_[i]),sizeof(SPCP(name)->_steps_[i]));
		}

		// Output data

		// Initialize
		unsigned int i = 0;
		const std::streamsize size = sizeof(SPCP(name)->_results_[0]);
		brgastro::vector<unsigned int> position(SPCP(name)->_num_dim_,0);

		while ( i < product(SPCP(name)->_resolutions_) )
		{
			out_file.write((char *)&(SPCP(name)->_results_(position)),size);
			i++;

			for(unsigned int d=0; d<SPCP(name)->_num_dim_; d++)
			{
				position[d]++;
				if(position[d] != SPCP(name)->_resolutions_[d])
					break;
				position[d] = 0;
				// If we get here, we'll go on to increase the next index by 1
			}
		}

		out_file.close();
		out_file.clear();
	}
#endif // Private methods

protected:

	// Protected methods
	// These are made protected instead of private so base classes can overload them
#if (1)

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	virtual const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0);
	}
	virtual const brgastro::unit_obj _inverse_units() const throw()
	{
		return brgastro::unit_obj(0);
	}

#endif // _BRG_USE_UNITS_

	// Long calculation function, which is used to generate the cache
	virtual const double _calculate(const brgastro::vector<double> & x) const
	{
		return 0;
	}

	// This function should be overloaded to provide a unique name for this cache
	virtual const std::string _name_base() const throw()
	{
		char name_base[BRG_CACHE_ND_NAME_SIZE] = "";
		return name_base;
	}

	// This function should be overloaded to call each cache of the same dimensionality
	// this cache depends upon in calculation. This is necessary in order to avoid critical
	// sections of the same name being called recursively.
	virtual void _load_cache_dependencies() const
	{
	}

#endif // Protected methods

public:

	// Public methods
#if (1)

	void set_file_name( const std::string new_name )
	{
		if(!SPCP(name)->_initialised_) SPP(name)->_init();
		SPP(name)->_file_name_ = new_name;
		if ( SPCP(name)->_loaded_ )
		{
			SPCP(name)->_unload();
		}
	} // void set_file_name()

	void set_range( const brgastro::vector<double> & new_mins, const brgastro::vector<double> & new_maxes,
			const brgastro::vector<double> & new_steps, const bool silent = false )
	{
		if(!SPCP(name)->_initialised_) SPP(name)->_init();

		// First we try to load, so we can see if there are any changes from
		// the existing cache
		if ( !SPCP(name)->_loaded_ )
			SPCP(name)->_load( true );

		// Check sizes of passed vectors
		if( (new_mins.size() != SPCP(name)->_num_dim_) || (new_maxes.size() != SPCP(name)->_num_dim_) ||
				(new_steps.size() != SPCP(name)->_num_dim_) )
		{
			throw std::runtime_error("ERROR: Incorrect sizes of vectors passed to set_range.\n");
		}

		// Go through variables, check if any are actually changed. If so, recalculate cache
		for(unsigned int i = 0; i < SPCP(name)->_num_dim_; i++)
		{
			if ( ( SPCP(name)->_mins_.at(i) != new_mins.at(i) ) || ( SPCP(name)->_maxes_.at(i) != new_maxes.at(i) )
					|| ( SPCP(name)->_steps_.at(i) != new_steps.at(i) ) )
			{
				SPP(name)->_mins_ = new_mins;
				SPP(name)->_maxes_ = new_maxes;
				SPP(name)->_steps_ = new_steps;

				SPCP(name)->_unload();
				SPCP(name)->_calc( silent );
				break;
			}
		}
	} // void set_range()

	const BRG_UNITS get( const brgastro::vector<double> & x, const bool silent = false ) const
	{

		brgastro::vector<double> xlo, xhi;
		brgastro::vector<unsigned int> x_i; // Lower nearby array points
#ifdef _BRG_USE_UNITS_
		BRG_UNITS result = SPCP(name)->_units(); // Ensure the result has the proper units
		result = 0;
#else
		double result = 0;
#endif

		if(!SPCP(name)->_initialised_) SPCP(name)->_init();

		// Load if necessary
		if ( !SPCP(name)->_loaded_ )
		{
			// Load any caches we depend upon before the critical section
			_load_cache_dependencies();

			// Critical section here, since we can't load multiple times simultaneously
			#pragma omp critical(load_brg_cache_nd)
			{
				try
				{
					SPCP(name)->_load( silent );
				}
				catch(const std::exception &e)
				{
					result = -1;
				}
			}
			if ( result == -1 )
			{
				throw std::runtime_error("ERROR: Could neither load " + SPCP(name)->_file_name_ + " nor calculate in brg_cache_nd::get()\n");
			}
		}


		x_i = bound(0,
				( ( x - SPCP(name)->_min_1_ ) / SPCP(name)->_step_1_ ),
				SPCP(name)->_resolution_1_ - 2 );

		xlo = SPCP(name)->_mins_ + SPCP(name)->_steps_ * x_i;
		xhi = SPCP(name)->_mins_ + SPCP(name)->_steps_ * ( x_i + 1 );

		unsigned int num_surrounding_points = 2;
		for(int i=0; i<SPCP(name)->_num_dim_-1; ++i ) num_surrounding_points*=2;

		result = 0;
		double total_weight = 0;
		brgastro::vector<unsigned int> position(SPCP(name)->_num_dim_,0);

		for(unsigned int j=0; j < num_surrounding_points; j++)
		{
			double weight = 1;
			unsigned int divisor = 1;
			for(unsigned int i=0; i < SPCP(name)->_num_dim_; i++)
			{
				if(divisible(j/divisor,2))
				{
					position[i] = x_i[i]+1;
					weight *= x[i]-xlo[i];
				}
				else
				{
					position[i] = x_i[i];
					weight *= xhi[i]-x[i];
				}
				divisor *= 2;
			}
			result += SPCP(name)->_results_.at(position) * weight;
			total_weight += weight;
		}

		result /= safe_d(total_weight);

		return result;

	} // get()

	// Recalculate function. Call if you want to overwrite a cache when something's changed in the code
	// (for instance, the _calculate() function has been altered)
	void recalc( const bool silent = false ) const
	{
		SPCP(name)->_unload();
		SPCP(name)->_calc(silent);
		SPCP(name)->_output(silent);
		SPCP(name)->_unload();
		SPCP(name)->_load(silent);
	}

	// Constructor
	brg_cache_nd() throw()
	{
		if(!SPCP(name)->_initialised_) SPP(name)->_init();
	}

	// Deconstructor
	virtual ~brg_cache_nd()
	{
	}

#endif // Public methods

}; // class brg_cache

} // namespace brgastro

// Undef macros
#undef SPP
#undef SPCP

#endif // __BRG_CACHE_ND_HPP_INCLUDED__
