/*
 * cache_4d.hpp
 *
 *  Created on: 25 Apr 2014
 *      Author: brg
 */

#ifndef __BRG_CACHE_4D_HPP_INCLUDED__
#define __BRG_CACHE_4D_HPP_INCLUDED__

#ifndef BRG_CACHE_ND_NAME_SIZE
#define BRG_CACHE_ND_NAME_SIZE 9 // Needs an end character, so will only actually allow 8 chars
#endif

#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "../brg_global.h"

#include "../brg_file_functions.h"
#include "../brg_units.h"

// Macro definitions

// SPCP: "Static Polymorphic Const Pointer"
// The constness isn't actually enforced, but this is for the reader's understanding
#define SPCP(name) static_cast<const name*>(this)

// SPP: "Static Polymorphic Pointer"
#define SPP(name) static_cast<name*>(this)

#define DECLARE_BRG_CACHE_4D_STATIC_VARS()		               \
	static double _min_1_, _max_1_, _step_1_;                  \
	static unsigned int _resolution_1_;                        \
	static double _min_2_, _max_2_, _step_2_;                  \
	static unsigned int _resolution_2_;                        \
	static double _min_3_, _max_3_, _step_3_;                  \
	static unsigned int _resolution_3_;                        \
	static double _min_4_, _max_4_, _step_4_;                  \
	static unsigned int _resolution_4_;                        \
	static std::vector< std::vector< std::vector< std::vector< double > > > > _results_;     \
											                   \
	static std::string _file_name_;                            \
	static unsigned int _version_number_;                      \
											                   \
	static bool _loaded_, _initialised_;

#define DEFINE_BRG_CACHE_4D_STATIC_VARS(class_name,init_min_1,init_max_1,init_step_1, \
                                        init_min_2, init_max_2, init_step_2,\
                                        init_min_3, init_max_3, init_step_3,\
                                        init_min_4, init_max_4, init_step_4)\
	double brgastro::class_name::_min_1_ = init_min_1;						\
	double brgastro::class_name::_max_1_ = init_max_1; 						\
	double brgastro::class_name::_step_1_ = init_step_1;					\
	unsigned int brgastro::class_name::_resolution_1_ = 0;					\
	double brgastro::class_name::_min_2_ = init_min_2;						\
	double brgastro::class_name::_max_2_ = init_max_2; 						\
	double brgastro::class_name::_step_2_ = init_step_2;					\
	unsigned int brgastro::class_name::_resolution_2_ = 0;					\
	double brgastro::class_name::_min_3_ = init_min_3;						\
	double brgastro::class_name::_max_3_ = init_max_3; 						\
	double brgastro::class_name::_step_3_ = init_step_3;					\
	unsigned int brgastro::class_name::_resolution_3_ = 0;					\
	double brgastro::class_name::_min_4_ = init_min_4;						\
	double brgastro::class_name::_max_4_ = init_max_4; 						\
	double brgastro::class_name::_step_4_ = init_step_4;					\
	unsigned int brgastro::class_name::_resolution_4_ = 0;					\
	bool brgastro::class_name::_loaded_ = false;							\
	bool brgastro::class_name::_initialised_ = false;						\
	std::string brgastro::class_name::_file_name_ = "";						\
	unsigned int brgastro::class_name::_version_number_ = 2;		        \
	std::vector< std::vector< std::vector< std::vector< double > > > > brgastro::class_name::_results_;

namespace brgastro
{

template<typename name>
class brg_cache_4d
{
private:

	// Private variables
#if (1)

	DECLARE_BRG_CACHE_4D_STATIC_VARS();

#endif // Private variables

	// Private methods
#if (1)
	const int _init() const throw()
	{
		if(SPCP(name)->_initialised_) return 0;

		SPCP(name)->_resolution_1_ = (unsigned int) max( ( ( SPCP(name)->_max_1_ - SPCP(name)->_min_1_ ) / safe_d(SPCP(name)->_step_1_)) + 1, 2);
		SPCP(name)->_resolution_2_ = (unsigned int) max( ( ( SPCP(name)->_max_2_ - SPCP(name)->_min_2_ ) / safe_d(SPCP(name)->_step_2_)) + 1, 2);
		SPCP(name)->_resolution_3_ = (unsigned int) max( ( ( SPCP(name)->_max_3_ - SPCP(name)->_min_3_ ) / safe_d(SPCP(name)->_step_3_)) + 1, 2);
		SPCP(name)->_resolution_4_ = (unsigned int) max( ( ( SPCP(name)->_max_4_ - SPCP(name)->_min_4_ ) / safe_d(SPCP(name)->_step_4_)) + 1, 2);
		SPCP(name)->_file_name_ = SPCP(name)->_name_base() + "_cache.bin";
		SPCP(name)->_version_number_ = 2; // This should be changed when there are changes to this code

		SPCP(name)->_initialised_ = true;

		return 0;
	}
	const int _load( const bool silent = false ) const
	{
		std::ifstream in_file;
		std::string file_data;
		bool need_to_calc = false;
		int loop_counter = 0;

		if ( SPCP(name)->_loaded_ )
			return 0;

		do
		{
			if ( loop_counter >= 2 )
			{
				if ( !silent )
					std::cerr << "ERROR: infinite loop detected trying to load " + SPCP(name)->_file_name_ + " in brgastro::brg_cache_2.\n";
				return INFINITE_LOOP_ERROR;
			}
			else
			{
				loop_counter++;
			}
			need_to_calc = false;

			if ( open_bin_file( in_file, SPCP(name)->_file_name_, true ) )
			{
				need_to_calc = true;
				if ( SPCP(name)->_calc( silent ) )
					return UNSPECIFIED_ERROR;
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
				if ( SPCP(name)->_calc( silent ) )
					return UNSPECIFIED_ERROR;
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}
			// Load range parameters;
			in_file.read((char *)&(SPCP(name)->_min_1_),sizeof(SPCP(name)->_min_1_));
			in_file.read((char *)&(SPCP(name)->_max_1_),sizeof(SPCP(name)->_max_1_));
			in_file.read((char *)&(SPCP(name)->_step_1_),sizeof(SPCP(name)->_step_1_));
			in_file.read((char *)&(SPCP(name)->_min_2_),sizeof(SPCP(name)->_min_2_));
			in_file.read((char *)&(SPCP(name)->_max_2_),sizeof(SPCP(name)->_max_2_));
			in_file.read((char *)&(SPCP(name)->_step_2_),sizeof(SPCP(name)->_step_2_));
			in_file.read((char *)&(SPCP(name)->_min_3_),sizeof(SPCP(name)->_min_3_));
			in_file.read((char *)&(SPCP(name)->_max_3_),sizeof(SPCP(name)->_max_3_));
			in_file.read((char *)&(SPCP(name)->_step_3_),sizeof(SPCP(name)->_step_3_));
			in_file.read((char *)&(SPCP(name)->_min_4_),sizeof(SPCP(name)->_min_4_));
			in_file.read((char *)&(SPCP(name)->_max_4_),sizeof(SPCP(name)->_max_4_));
			in_file.read((char *)&(SPCP(name)->_step_4_),sizeof(SPCP(name)->_step_4_));

			// Set up data
			SPCP(name)->_resolution_1_ = (unsigned int) max( ( ( SPCP(name)->_max_1_ - SPCP(name)->_min_1_ ) / safe_d(SPCP(name)->_step_1_)) + 1, 2);
			SPCP(name)->_resolution_2_ = (unsigned int) max( ( ( SPCP(name)->_max_2_ - SPCP(name)->_min_2_ ) / safe_d(SPCP(name)->_step_2_)) + 1, 2);
			SPCP(name)->_resolution_3_ = (unsigned int) max( ( ( SPCP(name)->_max_3_ - SPCP(name)->_min_3_ ) / safe_d(SPCP(name)->_step_3_)) + 1, 2);
			SPCP(name)->_resolution_4_ = (unsigned int) max( ( ( SPCP(name)->_max_4_ - SPCP(name)->_min_4_ ) / safe_d(SPCP(name)->_step_4_)) + 1, 2);
			make_array4d( SPCP(name)->_results_, SPCP(name)->_resolution_1_, SPCP(name)->_resolution_2_,
					SPCP(name)->_resolution_3_, SPCP(name)->_resolution_4_ );

			// Read in data

			// Initialise
			const std::streamsize size = sizeof(SPCP(name)->_results_[0][0][0][0]); // Store the size
			unsigned int i_1=0, i_2=0, i_3=0, i_4=0;

			while ( ( !in_file.eof() ) && ( i_4 < SPCP(name)->_resolution_4_ )
					&& (in_file) )
			{
				in_file.read((char *)&(SPCP(name)->_results_[i_1][i_2][i_3][i_4]),size);

				++i_1;
				if(i_1==SPCP(name)->_resolution_1_)
				{
					++i_2;
					if(i_2==SPCP(name)->_resolution_2_)
					{
						++i_3;
						if(i_3==SPCP(name)->_resolution_3_)
						{
							++i_4;
							i_3=0;
						}
						i_2=0;
					}
					i_1=0;
				}
			}

			// Check that it was all read properly
			if ( (i_4 != SPCP(name)->_resolution_4_) || (i_3 != 0) || (i_2 != 0) || (i_1 != 0)
					|| (!in_file) )
			{
				need_to_calc = true;
				if ( SPCP(name)->_calc( silent ) )
					return UNSPECIFIED_ERROR;
				SPCP(name)->_output();
				SPCP(name)->_unload();
				continue;
			}

		} while ( need_to_calc );

		// Finish up
		in_file.close();
		in_file.clear();
		SPCP(name)->_loaded_ = true;
		return 0;
	}
	const int _unload() const throw()
	{
		SPCP(name)->_loaded_ = false;
		SPCP(name)->_results_.clear();
		return 0;
	}
	const int _calc( const bool silent = false ) const
	{

		// Test that range is sane
		if ( ( SPCP(name)->_max_1_ <= SPCP(name)->_min_1_ ) || ( SPCP(name)->_step_1_ <= 0 ) ||
				 ( SPCP(name)->_max_2_ <= SPCP(name)->_min_2_ ) || ( SPCP(name)->_step_2_ <= 0 ) ||
				 ( SPCP(name)->_max_3_ <= SPCP(name)->_min_3_ ) || ( SPCP(name)->_step_3_ <= 0 ) ||
				 ( SPCP(name)->_max_4_ <= SPCP(name)->_min_4_ ) || ( SPCP(name)->_step_4_ <= 0 ) )
		{
			if ( !silent )
				std::cerr
						<< "ERROR: Bad range passed to brg_cache_4d::_calc() for " + SPCP(name)->_name_base() + "\n";
			return INVALID_ARGUMENTS_ERROR;
		}

		// Print a message that we're generating the cache if not silent
		if(!silent)
		{
			std::cout << "Generating " << SPCP(name)->_file_name_ << ". This may take some time.\n";
		}

		// Set up data
		SPCP(name)->_resolution_1_ = (unsigned int) max( ( ( SPCP(name)->_max_1_ - SPCP(name)->_min_1_ ) / safe_d(SPCP(name)->_step_1_)) + 1, 2);
		SPCP(name)->_resolution_2_ = (unsigned int) max( ( ( SPCP(name)->_max_2_ - SPCP(name)->_min_2_ ) / safe_d(SPCP(name)->_step_2_)) + 1, 2);
		SPCP(name)->_resolution_3_ = (unsigned int) max( ( ( SPCP(name)->_max_3_ - SPCP(name)->_min_3_ ) / safe_d(SPCP(name)->_step_3_)) + 1, 2);
		SPCP(name)->_resolution_4_ = (unsigned int) max( ( ( SPCP(name)->_max_4_ - SPCP(name)->_min_4_ ) / safe_d(SPCP(name)->_step_4_)) + 1, 2);
		if ( make_array4d( SPCP(name)->_results_, SPCP(name)->_resolution_1_, SPCP(name)->_resolution_2_,
				SPCP(name)->_resolution_3_, SPCP(name)->_resolution_4_ ) )
			return 1;

		// Calculate data
		bool bad_result = false;
		#pragma omp parallel for
		for ( unsigned int i_1 = 0; i_1 < SPCP(name)->_resolution_1_; ++i_1 )
		{
			double x_1 = SPCP(name)->_min_1_ + i_1*SPCP(name)->_step_1_;
			for( unsigned int i_2 = 0; i_2 < SPCP(name)->_resolution_2_; ++i_2)
			{
				double x_2 = SPCP(name)->_min_2_ + i_2*SPCP(name)->_step_2_;
				for( unsigned int i_3 = 0; i_3 < SPCP(name)->_resolution_3_; ++i_3)
				{
					double x_3 = SPCP(name)->_min_3_ + i_3*SPCP(name)->_step_3_;
					for( unsigned int i_4 = 0; i_4 < SPCP(name)->_resolution_4_; ++i_4)
					{
						double x_4 = SPCP(name)->_min_4_ + i_4*SPCP(name)->_step_4_;
						double result = 0;
						if(SPCP(name)->_calculate(x_1, x_2, x_3, x_4, result))
						{
							bad_result = true;
						}
						SPCP(name)->_results_[i_1][i_2][i_3][i_4] = result;
					}
				}
			}
		}

		if(bad_result) return UNSPECIFIED_ERROR;
		SPCP(name)->_loaded_ = true;

		// Print a message that we've finished generating the cache if not silent.
		if(!silent)
		{
			std::cout << "Finished generating " << SPCP(name)->_file_name_ << "!\n";
		}

		return 0;
	}
	const int _output( const bool silent = false ) const
	{

		std::ofstream out_file;
		std::string file_data;

		if ( !SPCP(name)->_loaded_ )
		{
			if ( SPCP(name)->_calc( silent ) )
				return LOWER_LEVEL_ERROR;
		}

		if ( open_bin_file( out_file, SPCP(name)->_file_name_, true ) )
			return LOWER_LEVEL_ERROR;

		// Output name and version

		std::string file_name = SPCP(name)->_name_base();
		unsigned int file_version = SPCP(name)->_version_number_;

		out_file.write(file_name.c_str(),BRG_CACHE_ND_NAME_SIZE);
		out_file.write((char *)&file_version,sizeof(file_version));

		// Output range parameters
		out_file.write((char *)&(SPCP(name)->_min_1_),sizeof(SPCP(name)->_min_1_));
		out_file.write((char *)&(SPCP(name)->_max_1_),sizeof(SPCP(name)->_max_1_));
		out_file.write((char *)&(SPCP(name)->_step_1_),sizeof(SPCP(name)->_step_1_));
		out_file.write((char *)&(SPCP(name)->_min_2_),sizeof(SPCP(name)->_min_2_));
		out_file.write((char *)&(SPCP(name)->_max_2_),sizeof(SPCP(name)->_max_2_));
		out_file.write((char *)&(SPCP(name)->_step_2_),sizeof(SPCP(name)->_step_2_));
		out_file.write((char *)&(SPCP(name)->_min_3_),sizeof(SPCP(name)->_min_3_));
		out_file.write((char *)&(SPCP(name)->_max_3_),sizeof(SPCP(name)->_max_3_));
		out_file.write((char *)&(SPCP(name)->_step_3_),sizeof(SPCP(name)->_step_3_));
		out_file.write((char *)&(SPCP(name)->_min_4_),sizeof(SPCP(name)->_min_4_));
		out_file.write((char *)&(SPCP(name)->_max_4_),sizeof(SPCP(name)->_max_4_));
		out_file.write((char *)&(SPCP(name)->_step_4_),sizeof(SPCP(name)->_step_4_));

		// Output data

		// Initialize
		const std::streamsize size = sizeof(SPCP(name)->_results_[0][0][0][0]);
		unsigned int i_1=0, i_2=0, i_3=0, i_4=0;

		while ( i_4<SPCP(name)->_resolution_4_ )
		{
			out_file.write((char *)&(SPCP(name)->_results_[i_1][i_2][i_3][i_4]),size);

			++i_1;
			if(i_1==(SPCP(name)->_resolution_1_))
			{
				++i_2;
				if(i_2==(SPCP(name)->_resolution_2_))
				{
					++i_3;
					if(i_3==(SPCP(name)->_resolution_3_))
					{
						++i_4;
						i_3 = 0;
					}
					i_2 = 0;
				}
				i_1 = 0;
			}
		}

		out_file.close();
		out_file.clear();

		return 0;
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

#endif // _BRG_USE_UNITS_

	const std::string _name_base() const throw()
	{
		char name_base[BRG_CACHE_ND_NAME_SIZE] = "";
		return name_base;
	}

#endif // Protected methods

public:

	// Public methods
#if (1)

	const int set_file_name( const std::string new_name )
	{
		if(!SPCP(name)->_initialised_) SPP(name)->_init();
		SPP(name)->_file_name_ = new_name;
		if ( SPCP(name)->_loaded_ )
		{
			return SPCP(name)->_unload();
		}
		return 0;
	} // const int set_file_name()

	const int set_range( const double new_min_1, const double new_max_1, const double new_step_1,
	         const double new_min_2, const double new_max_2, const double new_step_2,
 	         const double new_min_3, const double new_max_3, const double new_step_3,
 	         const double new_min_4, const double new_max_4, const double new_step_4,
			 const bool silent = false )
	{
		if(!SPCP(name)->_initialised_) SPP(name)->_init();

		// First we try to load, so we can see if there are any changes from
		// the existing cache
		if ( !SPCP(name)->_loaded_ )
			SPCP(name)->_load( true );

		// Go through variables, check if any are actually changed. If so, recalculate cache
		if ( ( SPCP(name)->_min_1_ != new_min_1 ) || ( SPCP(name)->_max_1_ != new_max_1 )
				|| ( SPCP(name)->_step_1_ != new_step_1 ) ||
			 ( SPCP(name)->_min_2_ != new_min_2 ) || ( SPCP(name)->_max_2_ != new_max_2 )
				|| ( SPCP(name)->_step_2_ != new_step_2 ) ||
			 ( SPCP(name)->_min_3_ != new_min_3 ) || ( SPCP(name)->_max_3_ != new_max_3 )
				|| ( SPCP(name)->_step_3_ != new_step_3 ) ||
			 ( SPCP(name)->_min_3_ != new_min_4 ) || ( SPCP(name)->_max_4_ != new_max_4 )
				|| ( SPCP(name)->_step_4_ != new_step_4 ))
		{
			SPP(name)->_min_1_ = new_min_1;
			SPP(name)->_max_1_ = new_max_1;
			SPP(name)->_step_1_ = new_step_1;
			SPP(name)->_min_2_ = new_min_2;
			SPP(name)->_max_2_ = new_max_2;
			SPP(name)->_step_2_ = new_step_2;
			SPP(name)->_min_3_ = new_min_3;
			SPP(name)->_max_3_ = new_max_3;
			SPP(name)->_step_3_ = new_step_3;
			SPP(name)->_min_4_ = new_min_4;
			SPP(name)->_max_4_ = new_max_4;
			SPP(name)->_step_4_ = new_step_4;

			if ( SPCP(name)->_unload() )
				return errorNOS( silent );
			if ( SPCP(name)->_calc( silent ) )
				return UNSPECIFIED_ERROR;
		}
		return 0;
	} // const int set_range()

	const int print( std::ostream & out, const bool silent = false ) const
	{

		if(!SPCP(name)->_initialised_) SPCP(name)->_init();

		// Load if necessary
		if ( !SPCP(name)->_loaded_ )
		{
			bool good_load=true;;
			// Critical section here, since we can't load multiple times simultaneously
			#pragma omp critical(load_brg_cache)
			{
				if ( SPCP(name)->_load( silent ) )
				{
					good_load = false;
				}
			}
			if ( !good_load )
			{
				std::string err = "ERROR: Could neither load " + SPCP(name)->_file_name_ + " nor calculate in brg_cache_4d::print()\n";
				if ( !silent )
					std::cerr << err;
				throw std::runtime_error(err);
			}
		}

		// Fill up header
		std::vector< std::string > header(6);
		header[0] = "#";
		header[1] = "x_1";
		header[2] = "x_2";
		header[3] = "x_3";
		header[4] = "x_4";
		header[5] = "y";

		// Fill up data
		std::vector< std::vector<std::string> > data(6);
		std::stringstream ss;
		for(unsigned int i_1=0; i_1<SPCP(name)->_resolution_1_; ++i_1)
		{
			for(unsigned int i_2=0; i_2<SPCP(name)->_resolution_2_; ++i_2)
			{
				for(unsigned int i_3=0; i_3<SPCP(name)->_resolution_3_; ++i_3)
				{
					for(unsigned int i_4=0; i_4<SPCP(name)->_resolution_4_; ++i_4)
					{
						data[0].push_back("");
						ss.str("");
						ss << SPCP(name)->_min_1_ + i_1*SPCP(name)->_step_1_;
						data[1].push_back(ss.str());
						ss.str("");
						ss << SPCP(name)->_min_2_ + i_2*SPCP(name)->_step_2_;
						data[2].push_back(ss.str());
						ss.str("");
						ss << SPCP(name)->_min_3_ + i_3*SPCP(name)->_step_3_;
						data[3].push_back(ss.str());
						ss.str("");
						ss << SPCP(name)->_min_4_ + i_4*SPCP(name)->_step_4_;
						data[4].push_back(ss.str());
						ss.str("");
						ss << SPCP(name)->_results_[i_1][i_2][i_3][i_4];
						data[5].push_back(ss.str());
					}
				}
			}
		}

		return print_table(out,data.size(),data[0].size(),header,data,false,silent);
	}

	const BRG_UNITS get( const double x_1, const double x_2, const double x_3, const double x_4,
			const bool silent = false ) const
	{

		double xlo_1, xhi_1;
		unsigned int xi_1; // Lower nearby array point
		double xlo_2, xhi_2;
		unsigned int xi_2; // Lower nearby array point
		double xlo_3, xhi_3;
		unsigned int xi_3; // Lower nearby array point
		double xlo_4, xhi_4;
		unsigned int xi_4; // Lower nearby array point
#ifdef _BRG_USE_UNITS_
		BRG_UNITS result = SPCP(name)->_units(); // Ensure the result has the proper units
		result = 0;
#else
		double result = 0;
#endif
		double total_weight = 0;

		if(!SPCP(name)->_initialised_) SPCP(name)->_init();

		// Load if necessary
		if ( !SPCP(name)->_loaded_ )
		{
			// Critical section here, since we can't load multiple times simultaneously
			#pragma omp critical(load_brg_cache)
			{
				if ( SPCP(name)->_load( silent ) )
				{
					result = -1;
				}
			}
			if ( result == -1 )
			{
				std::string err = "ERROR: Could neither load " + SPCP(name)->_file_name_ + " nor calculate in brg_cache_4d::get()\n";
				if ( !silent )
					std::cerr << err;
				throw std::runtime_error(err);
			}
		}

		xi_1 = (unsigned int)( ( x_1 - SPCP(name)->_min_1_ ) / SPCP(name)->_step_1_ );
		xi_1 = bound(0, xi_1, SPCP(name)->_resolution_1_ - 2 );
		xi_2 = (unsigned int)( ( x_2 - SPCP(name)->_min_2_ ) / SPCP(name)->_step_2_ );
		xi_2 = bound(0, xi_2, SPCP(name)->_resolution_2_ - 2 );
		xi_3 = (unsigned int)( ( x_3 - SPCP(name)->_min_3_ ) / SPCP(name)->_step_3_ );
		xi_3 = bound(0, xi_3, SPCP(name)->_resolution_3_ - 2 );
		xi_4 = (unsigned int)( ( x_4 - SPCP(name)->_min_4_ ) / SPCP(name)->_step_4_ );
		xi_4 = bound(0, xi_4, SPCP(name)->_resolution_4_ - 2 );

		xlo_1 = SPCP(name)->_min_1_ + SPCP(name)->_step_1_ * xi_1;
		xhi_1 = SPCP(name)->_min_1_ + SPCP(name)->_step_1_ * ( xi_1 + 1 );
		xlo_2 = SPCP(name)->_min_2_ + SPCP(name)->_step_2_ * xi_2;
		xhi_2 = SPCP(name)->_min_2_ + SPCP(name)->_step_2_ * ( xi_2 + 1 );
		xlo_3 = SPCP(name)->_min_3_ + SPCP(name)->_step_3_ * xi_3;
		xhi_3 = SPCP(name)->_min_3_ + SPCP(name)->_step_3_ * ( xi_3 + 1 );
		xlo_4 = SPCP(name)->_min_4_ + SPCP(name)->_step_4_ * xi_4;
		xhi_4 = SPCP(name)->_min_4_ + SPCP(name)->_step_4_ * ( xi_4 + 1 );

		total_weight = (xhi_1-xlo_1)*(xhi_2-xlo_2)*(xhi_3-xlo_3)*(xhi_4-xlo_4);

		result += SPCP(name)->_results_[xi_1][xi_2][xi_3][xi_4] *
				(xhi_1-x_1)*(xhi_2-x_2)*(xhi_3-x_3)*(xhi_4-x_4);
		result += SPCP(name)->_results_[xi_1+1][xi_2][xi_3][xi_4] *
				(x_1-xlo_1)*(xhi_2-x_2)*(xhi_3-x_3)*(xhi_4-x_4);

		result += SPCP(name)->_results_[xi_1][xi_2+1][xi_3][xi_4] *
				(xhi_1-x_1)*(x_2-xlo_2)*(xhi_3-x_3)*(xhi_4-x_4);
		result += SPCP(name)->_results_[xi_1+1][xi_2+1][xi_3][xi_4] *
				(x_1-xlo_1)*(x_2-xlo_2)*(xhi_3-x_3)*(xhi_4-x_4);


		result += SPCP(name)->_results_[xi_1][xi_2][xi_3+1][xi_4] *
				(xhi_1-x_1)*(xhi_2-x_2)*(x_3-xlo_3)*(xhi_4-x_4);
		result += SPCP(name)->_results_[xi_1+1][xi_2][xi_3+1][xi_4] *
				(x_1-xlo_1)*(xhi_2-x_2)*(x_3-xlo_3)*(xhi_4-x_4);

		result += SPCP(name)->_results_[xi_1][xi_2+1][xi_3+1][xi_4] *
				(xhi_1-x_1)*(x_2-xlo_2)*(x_3-xlo_3)*(xhi_4-x_4);
		result += SPCP(name)->_results_[xi_1+1][xi_2+1][xi_3+1][xi_4] *
				(x_1-xlo_1)*(x_2-xlo_2)*(x_3-xlo_3)*(xhi_4-x_4);




		result += SPCP(name)->_results_[xi_1][xi_2][xi_3][xi_4+1] *
				(xhi_1-x_1)*(xhi_2-x_2)*(xhi_3-x_3)*(x_4-xlo_4);
		result += SPCP(name)->_results_[xi_1+1][xi_2][xi_3][xi_4+1] *
				(x_1-xlo_1)*(xhi_2-x_2)*(xhi_3-x_3)*(x_4-xlo_4);

		result += SPCP(name)->_results_[xi_1][xi_2+1][xi_3][xi_4+1] *
				(xhi_1-x_1)*(x_2-xlo_2)*(xhi_3-x_3)*(x_4-xlo_4);
		result += SPCP(name)->_results_[xi_1+1][xi_2+1][xi_3][xi_4+1] *
				(x_1-xlo_1)*(x_2-xlo_2)*(xhi_3-x_3)*(x_4-xlo_4);


		result += SPCP(name)->_results_[xi_1][xi_2][xi_3+1][xi_4+1] *
				(xhi_1-x_1)*(xhi_2-x_2)*(x_3-xlo_3)*(x_4-xlo_4);
		result += SPCP(name)->_results_[xi_1+1][xi_2][xi_3+1][xi_4+1] *
				(x_1-xlo_1)*(xhi_2-x_2)*(x_3-xlo_3)*(x_4-xlo_4);

		result += SPCP(name)->_results_[xi_1][xi_2+1][xi_3+1][xi_4+1] *
				(xhi_1-x_1)*(x_2-xlo_2)*(x_3-xlo_3)*(x_4-xlo_4);
		result += SPCP(name)->_results_[xi_1+1][xi_2+1][xi_3+1][xi_4+1] *
				(x_1-xlo_1)*(x_2-xlo_2)*(x_3-xlo_3)*(x_4-xlo_4);

		result /= safe_d(total_weight);

		return result;

	} // get()

	// Recalculate function. Call if you want to overwrite a cache when something's changed in the code
	// (for instance, the _calculate() function has been altered)
	const int recalc( const bool silent = false ) const
	{
		SPCP(name)->_unload();
		if(SPCP(name)->_calc(silent)) return LOWER_LEVEL_ERROR;
		if(SPCP(name)->_output(silent)) return LOWER_LEVEL_ERROR;
		SPCP(name)->_unload();
		if(SPCP(name)->_load(silent)) return LOWER_LEVEL_ERROR;
		return 0;
	}

	// Constructor
	brg_cache_4d() throw()
	{
	}

	// Deconstructor
	virtual ~brg_cache_4d()
	{
	}

#endif // Public methods

}; // class brg_cache

} // namespace brgastro

// Undef macros
#undef SPP
#undef SPCP

#endif // __BRG_CACHE_2D_HPP_INCLUDED__
