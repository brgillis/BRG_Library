/**       @file tNFW_caches.cpp
 *
 *     Project: brg
 *        Path: /brg/density_profile/tNFW_caches.cpp
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../brg_global.h"

#include "../brg_astro.h"
#include "../brg_units.h"
#include "tNFW_profile.h"

#include "tNFW_caches.h"

// Initialisation for brgastro::tNFW_sig_cache
#if (1)
bool brgastro::tNFW_sig_cache::loaded = false;
std::string brgastro::tNFW_sig_cache::file_name = "brgastro_NFW_sig_cache.dat";
double brgastro::tNFW_sig_cache::halo_z_min = 0.1; // These values will only be used if the cache needs to be created and they aren't overridden
double brgastro::tNFW_sig_cache::halo_z_max = 1.3;
double brgastro::tNFW_sig_cache::halo_z_step = 0.1;
int brgastro::tNFW_sig_cache::halo_z_res = (int)( ( tNFW_sig_cache::halo_z_max
		- tNFW_sig_cache::halo_z_min ) / tNFW_sig_cache::halo_z_step ) + 1;
double brgastro::tNFW_sig_cache::halo_m_min = 10;
double brgastro::tNFW_sig_cache::halo_m_max = 16;
double brgastro::tNFW_sig_cache::halo_m_step = 0.01;
int brgastro::tNFW_sig_cache::halo_m_res = (int)( ( tNFW_sig_cache::halo_m_max
		- tNFW_sig_cache::halo_m_min ) / tNFW_sig_cache::halo_m_step ) + 1;
double brgastro::tNFW_sig_cache::r_min = 1 * unitconv::kpctom;
double brgastro::tNFW_sig_cache::r_max = 5000 * unitconv::kpctom;
double brgastro::tNFW_sig_cache::r_step = 1 * unitconv::kpctom;
int brgastro::tNFW_sig_cache::r_res = (int)( ( r_max - r_min ) / r_step ) + 1;
std::string brgastro::tNFW_sig_cache::header_string = "# NFW_sig v1.0";
int brgastro::tNFW_sig_cache::sig_digits = 8;
std::vector< std::vector< std::vector< double > > > brgastro::tNFW_sig_cache::signal;
#endif // End initialization for brgastro::tNFW_sig_cache

// Initialization for brgastro::tNFW_offset_sig_cache
#if (1)
bool brgastro::tNFW_offset_sig_cache::loaded = false;
std::string brgastro::tNFW_offset_sig_cache::file_name =
		"brgastro_NFW_offset_sig_cache.dat";
double brgastro::tNFW_offset_sig_cache::halo_z_min = 0.1; // These values will only be used if the cache needs to be created and they aren't overridden
double brgastro::tNFW_offset_sig_cache::halo_z_max = 1.3;
double brgastro::tNFW_offset_sig_cache::halo_z_step = 0.1;
int brgastro::tNFW_offset_sig_cache::halo_z_res =
		(int)( ( tNFW_offset_sig_cache::halo_z_max
				- tNFW_offset_sig_cache::halo_z_min )
				/ tNFW_offset_sig_cache::halo_z_step ) + 1;
double brgastro::tNFW_offset_sig_cache::halo_m_min = 10;
double brgastro::tNFW_offset_sig_cache::halo_m_max = 16;
double brgastro::tNFW_offset_sig_cache::halo_m_step = 0.01;
int brgastro::tNFW_offset_sig_cache::halo_m_res =
		(int)( ( tNFW_offset_sig_cache::halo_m_max
				- tNFW_offset_sig_cache::halo_m_min )
				/ tNFW_offset_sig_cache::halo_m_step ) + 1;
double brgastro::tNFW_offset_sig_cache::r_min = 1 * unitconv::kpctom;
double brgastro::tNFW_offset_sig_cache::r_max = 5000 * unitconv::kpctom;
double brgastro::tNFW_offset_sig_cache::r_step = 1 * unitconv::kpctom;
int brgastro::tNFW_offset_sig_cache::r_res =
		(int)( ( tNFW_offset_sig_cache::r_max - tNFW_offset_sig_cache::r_min )
				/ tNFW_offset_sig_cache::r_step ) + 1;
double brgastro::tNFW_offset_sig_cache::offset_r_min = -1;
double brgastro::tNFW_offset_sig_cache::offset_r_max = 4;
double brgastro::tNFW_offset_sig_cache::offset_r_step = 0.1;
int brgastro::tNFW_offset_sig_cache::offset_r_res =
		(int)( ( tNFW_offset_sig_cache::offset_r_max
				- tNFW_offset_sig_cache::offset_r_min )
				/ tNFW_offset_sig_cache::offset_r_step ) + 1;
std::string brgastro::tNFW_offset_sig_cache::header_string =
		"# NFW_offset_sig v1.0";
int brgastro::tNFW_offset_sig_cache::sig_digits = 8;
std::vector< std::vector< std::vector< std::vector< double > > > > brgastro::tNFW_offset_sig_cache::signal;
#endif // End initialization for brgastro::tNFW_offset_sig_cache

// Initialization for brgastro::tNFW_group_sig_cache
#if (1)
bool brgastro::tNFW_group_sig_cache::loaded = false;
std::string brgastro::tNFW_group_sig_cache::file_name =
		"brgastro::tNFW_group_sig_cache.dat";
double brgastro::tNFW_group_sig_cache::halo_z_min = 0.1; // These values will only be used if the cache needs to be created and they aren't overridden
double brgastro::tNFW_group_sig_cache::halo_z_max = 1.3;
double brgastro::tNFW_group_sig_cache::halo_z_step = 0.1;
int brgastro::tNFW_group_sig_cache::halo_z_res =
		(int)( ( tNFW_group_sig_cache::halo_z_max
				- tNFW_group_sig_cache::halo_z_min )
				/ tNFW_group_sig_cache::halo_z_step ) + 1;
double brgastro::tNFW_group_sig_cache::halo_m_min = 10;
double brgastro::tNFW_group_sig_cache::halo_m_max = 16;
double brgastro::tNFW_group_sig_cache::halo_m_step = 0.01;
int brgastro::tNFW_group_sig_cache::halo_m_res =
		(int)( ( tNFW_group_sig_cache::halo_m_max
				- tNFW_group_sig_cache::halo_m_min )
				/ tNFW_group_sig_cache::halo_m_step ) + 1;
double brgastro::tNFW_group_sig_cache::r_min = 1 * unitconv::kpctom;
double brgastro::tNFW_group_sig_cache::r_max = 5000 * unitconv::kpctom;
double brgastro::tNFW_group_sig_cache::r_step = 1 * unitconv::kpctom;
int brgastro::tNFW_group_sig_cache::r_res =
		(int)( ( tNFW_group_sig_cache::r_max - tNFW_group_sig_cache::r_min )
				/ tNFW_group_sig_cache::r_step ) + 1;
double brgastro::tNFW_group_sig_cache::group_c_min = 2;
double brgastro::tNFW_group_sig_cache::group_c_max = 8;
double brgastro::tNFW_group_sig_cache::group_c_step = 1;
int brgastro::tNFW_group_sig_cache::group_c_res =
		(int)( ( tNFW_group_sig_cache::group_c_max
				- tNFW_group_sig_cache::group_c_min )
				/ tNFW_group_sig_cache::group_c_step ) + 1;
std::string brgastro::tNFW_group_sig_cache::header_string =
		"# NFW_group_sig v1.0";
int brgastro::tNFW_group_sig_cache::sig_digits = 8;
std::vector< std::vector< std::vector< std::vector< double > > > > brgastro::tNFW_group_sig_cache::signal;
#endif // End initialization for brgastro::tNFW_group_sig_cache



// brgastro::dfa_cache class methods
#if (1)
const int brgastro::dfa_cache::_calculate( const double in_params, double & out_params ) const
{
	try
	{
		out_params = brgastro::integrate_add( 0, in_params );
	}
	catch( const std::exception &e)
	{
		return LOWER_LEVEL_ERROR;
	}

	return 0;
}

#endif // end brgastro::dfa_cache functions

// brgastro::add_cache class methods
const int brgastro::add_cache::_calculate( const brgastro::vector<double> in_params,
		double  & out_params ) const
{
	try
	{
		out_params = brgastro::integrate_add(in_params.at(0),in_params.at(1));
	}
	catch(const std::out_of_range &e)
	{
		return INVALID_ARGUMENTS_ERROR;
	}
	return 0;
}

const int brgastro::tfa_cache::_calculate( const double in_params, double & out_params ) const
{
	try
	{
		out_params = -brgastro::integrate_ltd( 0, brgastro::zfa( in_params ) ) / c;
	}
	catch(const std::exception &e)
	{
		std::cerr << "ERROR: Could not calculate cache for " << _name_base() << "\n"
				<< "Exception: " << e.what() << "\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

// brgastro::tNFW_sig_cache class methods
#if (1)

const int brgastro::tNFW_sig_cache::load( const bool silent )
{
	std::ifstream in_file;
	std::string file_data;
	bool need_to_calc = false;
	int loop_counter = 0;
	double temp_data;
	int i, j, k;

	if ( loaded )
		return 0;

	do
	{
		if ( loop_counter >= 2 )
		{
			if ( !silent )
				std::cerr
						<< "ERROR: infinite loop detected in brgastro::tNFW_sig_cache.\n";
			return INFINITE_LOOP_ERROR;
		}
		else
		{
			loop_counter++;
		}
		need_to_calc = false;

		if ( open_file( in_file, file_name ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Check that it has the right header
		getline( in_file, file_data );
		if ( file_data.compare( header_string ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Trim out any other commented lines
		if ( trim_comments_all_at_top( in_file ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Load range parameters;
		if ( !( in_file >> halo_z_min >> halo_z_max >> halo_z_step
				>> halo_m_min >> halo_m_max >> halo_m_step >> r_min >> r_max
				>> r_step ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Set up data
		halo_z_res = (int)( ( halo_z_max - halo_z_min ) / halo_z_step ) + 1;
		halo_m_res = (int)( ( halo_m_max - halo_m_min ) / halo_m_step ) + 1;
		r_res = (int)( ( r_max - r_min ) / r_step ) + 1;
		if ( make_array3d( signal, halo_z_res, halo_m_res, r_res ) )
			return 1;

		// Read in data
		i = j = k = 0;
		while ( !in_file.eof() )
		{
			in_file >> temp_data;
			signal[i][j][k] = temp_data;
			k++;
			if ( k >= r_res )
			{
				k = 0;
				j++;
				if ( j >= halo_m_res )
				{
					j = 0;
					i++;
					if ( i >= halo_z_res )
						break;
				}
			}
		}

		// Check that it was all read properly
		if ( i < halo_z_res )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

	} while ( need_to_calc );

	// Finish up
	in_file.close();
	in_file.clear();
	loaded = true;
	return 0;
}

const int brgastro::tNFW_sig_cache::unload()
{
	loaded = false;
	signal.clear();
	return 0;
}

const int brgastro::tNFW_sig_cache::calc( const bool silent )
{
	int i, j, k;
	double mass;

	// Test that range is sane
	if ( ( halo_z_max <= halo_z_min ) || ( halo_z_step <= 0 )
			|| ( halo_m_max <= halo_m_min ) || ( halo_m_step <= 0 )
			|| ( r_max <= r_min ) || ( r_step <= 0 ) )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Bad range passed to tNFW_sig_cache::calc(silent)\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	// Set up data
	halo_z_res = (int)( ( halo_z_max - halo_z_min ) / halo_z_step ) + 1;
	halo_m_res = (int)( ( halo_m_max - halo_m_min ) / halo_z_step ) + 1;
	r_res = (int)( ( r_max - r_min ) / r_step ) + 1;
	if ( make_array3d( signal, halo_z_res, halo_m_res, r_res ) )
		return 1;

	i = j = k = 0;

	for ( double z = halo_z_min; z < halo_z_max; z += halo_z_step )
	{
		for ( double m = halo_m_min; m < halo_m_max; m += halo_m_step )
		{
			mass = std::pow( 10, m ) * unitconv::Msuntokg;
			brgastro::tNFW_profile profile(mass,z);
			for ( double r = r_min; r < r_max; r += r_step )
			{

				signal[i][j][k] = profile.deltasigma( r );
				k++;
			}
			k = 0;
			j++;
		}
		j = k = 0;
		i++;
	}

	return 0;
}

const int brgastro::tNFW_sig_cache::output( const bool silent )
{
	std::ofstream out_file;
	std::string file_data;

	if ( !loaded )
	{
		if ( calc( silent ) )
			return 1;
	}

	if ( open_file( out_file, file_name ) )
		return 1;

	// Output header
	out_file << header_string << "\n#\n";

	// Output range
	out_file << halo_z_min << "\t" << halo_z_max << "\t" << halo_z_step << "\n"
			<< halo_m_min << "\t" << halo_m_max << "\t" << halo_m_step << "\n"
			<< r_min << "\t" << r_max << "\t" << r_step << "\n";

	// Output data
	for ( int i = 0; i < halo_z_res; i++ )
	{
		for ( int j = 0; j < halo_m_res; j++ )
		{
			for ( int k = 0; k < r_res; k++ )
			{
				if ( !( out_file << signal[i][j][k] << "\t" ) )
					return errorNOS();
			}
			out_file << "\n";
		}
	}

	out_file.close();
	out_file.clear();

	return 0;

}

const int brgastro::tNFW_sig_cache::set_file_name( const std::string new_name )
{
	file_name = new_name;
	if ( loaded )
	{
		return unload();
	}
	return 0;
}

const int brgastro::tNFW_sig_cache::set_range( const double new_halo_z_min,
		const double new_halo_z_max, const double new_halo_z_step,
		const double new_halo_m_min, const double new_halo_m_max,
		const double new_halo_m_step, const double new_r_min,
		const double new_r_max, const double new_r_step, const bool silent )
{
	// First we try to load
	if ( !loaded )
		load( silent );

	// Go through variables, check if any are actually changed. If so, recalculate cache
	if ( ( halo_z_min != new_halo_z_min ) || ( halo_z_max != new_halo_z_max )
			|| ( halo_z_step != new_halo_z_step )
			|| ( halo_m_min != new_halo_m_min )
			|| ( halo_m_max != new_halo_m_max )
			|| ( halo_m_step != new_halo_m_step )
			|| ( halo_z_min != new_r_min ) || ( halo_z_max != new_r_max )
			|| ( halo_z_step != new_r_step ) )
	{
		halo_z_min = new_halo_z_min;
		halo_z_max = new_halo_z_max;
		halo_z_step = new_halo_z_step;
		halo_m_min = new_halo_m_min;
		halo_m_max = new_halo_m_max;
		halo_m_step = new_halo_m_step;
		r_min = new_r_min;
		r_max = new_r_max;
		r_step = new_r_step;

		if ( unload() )
			return errorNOS( silent );
		if ( calc( silent ) )
			return 1;
	}
	return 0;
}

const int brgastro::tNFW_sig_cache::set_precision( const int new_precision,
		const bool silent )
{
	if ( new_precision > 0 )
	{
		sig_digits = min( new_precision, DBL_MAX_PRECISION );
		return 0;
	}
	else
	{
		if ( !silent )
			std::cerr << "ERROR: Precision for tNFW_sig_cache must be > 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
}

const BRG_UNITS brgastro::tNFW_sig_cache::get( const double z,
		const BRG_MASS m, const BRG_DISTANCE r, const bool silent )
{
	double zlo, zhi, mlo, mhi, rlo, rhi;
	int z_i, m_i, r_i; // Lower nearby array points
	BRG_DISTANCE rwlo, rwhi;
	double mwlo, mwhi, zwlo, zwhi;
	double lm = log10( m * unitconv::kgtoMsun );
	double result = 0;

	if ( !loaded )
	{
		#pragma omp critical(tNFW_sig_load)
		{
			if ( load( silent ) )
			{
				result = -1;
			}
		}

		if ( result == -1 )
		{
			if ( !silent )
				std::cerr
						<< "ERROR: Could neither load nor calculate tNFW_sig_cache!\n";
			return result;
		}
	}

	z_i = (int)( ( z - halo_z_min ) / halo_z_step );
	z_i = max( z_i, 0 );
	z_i = min( z_i, halo_z_res - 2 );

	zlo = halo_z_min + halo_z_step * z_i;
	zhi = halo_z_min + halo_z_step * ( z_i + 1 );
	zwlo = z - zlo;
	zwhi = zhi - z;

	m_i = (int)( ( lm - halo_m_min ) / halo_m_step );
	m_i = max( m_i, 0 );
	m_i = min( m_i, halo_m_res - 2 );

	mlo = halo_m_min + halo_m_step * m_i;
	mhi = halo_m_min + halo_m_step * ( m_i + 1 );
	mwlo = lm - mlo;
	mwhi = mhi - lm;

	r_i = (int)( ( r - r_min ) / r_step );
	r_i = max( r_i, 0 );
	r_i = min( r_i, r_res - 2 );

	rlo = r_min + r_step * r_i;
	rhi = r_min + r_step * ( r_i + 1 );
	rwlo = r - rlo;
	rwhi = rhi - r;

	double totweight = halo_z_step * halo_m_step * r_step;

#ifdef _BRG_USE_UNITS_
	result = BRG_UNITS(0,-2,0,1,0,0,0); // To get the right units
#else
	result = 0;
#endif

	result += signal[z_i][m_i][r_i] * zwhi * mwhi * rwhi;
	result += signal[z_i][m_i][r_i + 1] * zwhi * mwhi * rwlo;
	result += signal[z_i][m_i + 1][r_i] * zwhi * mwlo * rwhi;
	result += signal[z_i][m_i + 1][r_i + 1] * zwhi * mwlo * rwlo;
	result += signal[z_i + 1][m_i][r_i] * zwlo * mwhi * rwhi;
	result += signal[z_i + 1][m_i][r_i + 1] * zwlo * mwhi * rwlo;
	result += signal[z_i + 1][m_i + 1][r_i] * zwlo * mwlo * rwhi;
	result += signal[z_i + 1][m_i + 1][r_i + 1] * zwlo * mwlo * rwlo;

	result /= totweight;

	return result;

}

#endif // end brgastro::tNFW_sig_cache functions

// brgastro::NFW_offset_sig_cache class methods
#if (1)

const int brgastro::tNFW_offset_sig_cache::load( const bool silent )
{
	std::ifstream in_file;
	std::string file_data;
	bool need_to_calc = false;
	int loop_counter = 0;
	double temp_data;
	int i, j, k, l;

	if ( loaded )
		return 0;

	do
	{
		if ( loop_counter >= 2 )
		{
			if ( !silent )
				std::cerr
						<< "ERROR: infinite loop detected in brgastro::NFW_offset_sig_cache.\n";
			return INFINITE_LOOP_ERROR;
		}
		else
		{
			loop_counter++;
		}
		need_to_calc = false;

		if ( open_file( in_file, file_name ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Check that it has the right header
		getline( in_file, file_data );
		if ( file_data.compare( header_string ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Trim out any other commented lines
		if ( trim_comments_all_at_top( in_file ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Load range parameters;
		if ( !( in_file >> halo_z_min >> halo_z_max >> halo_z_step
				>> halo_m_min >> halo_m_max >> halo_m_step >> r_min >> r_max
				>> r_step >> offset_r_min >> offset_r_max >> offset_r_step ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Set up data
		halo_z_res = (int)( ( halo_z_max - halo_z_min ) / halo_z_step ) + 1;
		halo_m_res = (int)( ( halo_m_max - halo_m_min ) / halo_m_step ) + 1;
		r_res = (int)( ( r_max - r_min ) / r_step ) + 1;
		offset_r_res = (int)( ( offset_r_max - offset_r_min ) / offset_r_step )
				+ 1;
		if ( make_array4d( signal, halo_z_res, halo_m_res, r_res,
				offset_r_res ) )
			return 1;

		// Read in data
		i = j = k = l = 0;
		while ( !in_file.eof() )
		{
			in_file >> temp_data;
			signal[i][j][k][l] = temp_data;
			l++;
			if ( l >= offset_r_res )
			{
				l = 0;
				k++;
				if ( k >= r_res )
				{
					k = 0;
					j++;
					if ( j >= halo_m_res )
					{
						j = 0;
						i++;
						if ( i >= halo_z_res )
							break;
					}
				}
			}
		}
		// Check that it was all read properly
		if ( i < halo_z_res )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

	} while ( need_to_calc );

	// Finish up
	in_file.close();
	in_file.clear();
	loaded = true;
	return 0;
}

const int brgastro::tNFW_offset_sig_cache::unload()
{
	loaded = false;
	signal.clear();
	return 0;
}

const int brgastro::tNFW_offset_sig_cache::calc( const bool silent )
{
	double z, lm, r, l_offset_r;
	double mass, offset_r;

	// Test that range is sane
	if ( ( halo_z_max <= halo_z_min ) || ( halo_z_step <= 0 )
			|| ( halo_m_max <= halo_m_min ) || ( halo_m_step <= 0 )
			|| ( r_max <= r_min ) || ( r_step <= 0 )
			|| ( offset_r_max <= offset_r_min ) || ( offset_r_step <= 0 ) )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Bad range passed to tNFW_offset_sig_cache::calc(silent)\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	// Set up data
	halo_z_res = (int)( ( halo_z_max - halo_z_min ) / halo_z_step ) + 1;
	halo_m_res = (int)( ( halo_m_max - halo_m_min ) / halo_z_step ) + 1;
	r_res = (int)( ( r_max - r_min ) / r_step ) + 1;
	offset_r_res = (int)( ( offset_r_max - offset_r_min ) / offset_r_step )
			+ 1;
	if ( make_array4d( signal, halo_z_res, halo_m_res, r_res, offset_r_res ) )
		return 1;

	z = halo_z_min;
	for ( int i = 0; i < halo_z_res; i++ )
	{
		lm = halo_m_min;
		for ( int j = 0; j < halo_m_res; j++ )
		{
			mass = std::pow( 10, lm ) * unitconv::Msuntokg;
			brgastro::tNFW_profile profile(mass,z);
			r = r_min;
			for ( int k = 0; k < r_res; k++ )
			{
				l_offset_r = offset_r_min;
				for ( int l = 0; l < offset_r_res; l++ )
				{
					offset_r = std::pow( 10, l_offset_r ) * unitconv::kpctom;
					signal[i][j][k][l] = profile.offset_WLsig( r, offset_r );
					l_offset_r += offset_r_step;
				}
				r += r_step;
			}
			lm += halo_m_step;
		}
		z += halo_z_step;
	}

	return 0;
}

const int brgastro::tNFW_offset_sig_cache::output( const bool silent )
{
	std::ofstream out_file;
	std::string file_data;

	if ( !loaded )
	{
		if ( calc( silent ) )
			return 1;
	}

	if ( open_file( out_file, file_name ) )
		return 1;

	// Output header
	out_file << header_string << "\n#\n";

	// Output range
	out_file << halo_z_min << "\t" << halo_z_max << "\t" << halo_z_step << "\n"
			<< halo_m_min << "\t" << halo_m_max << "\t" << halo_m_step << "\n"
			<< r_min << "\t" << r_max << "\t" << r_step << "\n" << offset_r_min
			<< "\t" << offset_r_max << "\t" << offset_r_step << "\n";

	// Output data
	for ( int i = 0; i < halo_z_res; i++ )
	{
		for ( int j = 0; j < halo_m_res; j++ )
		{
			for ( int k = 0; k < r_res; k++ )
			{
				for ( int l = 0; l < offset_r_res; l++ )
				{
					if ( !( out_file << signal[i][j][k][l] << "\n" ) )
						return errorNOS();
				}
			}
		}
	}

	out_file.close();
	out_file.clear();

	return 0;

}

const int brgastro::tNFW_offset_sig_cache::set_file_name(
		const std::string new_name )
{
	file_name = new_name;
	if ( loaded )
	{
		if ( unload() )
			return errorNOS();
	}
	return 0;
}

const int brgastro::tNFW_offset_sig_cache::set_range(
		const double new_halo_z_min, const double new_halo_z_max,
		const double new_halo_z_step, const double new_halo_m_min,
		const double new_halo_m_max, const double new_halo_m_step,
		const double new_r_min, const double new_r_max,
		const double new_r_step, const double new_offset_r_min,
		const double new_offset_r_max, const double new_offset_r_step,
		const bool silent )
{
	// First we try to load
	if ( !loaded )
		load( silent );

	// Go through variables, check if any are actually changed. If so, recalculate cache
	if ( ( halo_z_min != new_halo_z_min ) || ( halo_z_max != new_halo_z_max )
			|| ( halo_z_step != new_halo_z_step )
			|| ( halo_m_min != new_halo_m_min )
			|| ( halo_m_max != new_halo_m_max )
			|| ( halo_m_step != new_halo_m_step )
			|| ( halo_z_min != new_r_min ) || ( halo_z_max != new_r_max )
			|| ( halo_z_step != new_r_step )
			|| ( halo_z_min != new_offset_r_min )
			|| ( halo_z_max != new_offset_r_max )
			|| ( halo_z_step != new_offset_r_step ) )
	{
		halo_z_min = new_halo_z_min;
		halo_z_max = new_halo_z_max;
		halo_z_step = new_halo_z_step;
		halo_m_min = new_halo_m_min;
		halo_m_max = new_halo_m_max;
		halo_m_step = new_halo_m_step;
		r_min = new_r_min;
		r_max = new_r_max;
		r_step = new_r_step;
		offset_r_min = new_offset_r_min;
		offset_r_max = new_offset_r_max;
		offset_r_step = new_offset_r_step;

		if ( unload() )
			return errorNOS( silent );
		if ( calc( silent ) )
			return 1;
	}
	return 0;
}

const int brgastro::tNFW_offset_sig_cache::set_precision(
		const int new_precision, const bool silent )
{
	if ( new_precision > 0 )
	{
		sig_digits = min( new_precision, DBL_MAX_PRECISION );
		return 0;
	}
	else
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Precision for tNFW_offset_sig_cache must be > 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
}

const BRG_UNITS brgastro::tNFW_offset_sig_cache::get( const double z,
		const BRG_MASS m, const BRG_DISTANCE r, const BRG_DISTANCE offset_r,
		const bool silent )
{
	double zlo, zhi, mlo, mhi;
	BRG_DISTANCE rlo, rhi, orlo, orhi;
	int z_i, m_i, r_i, or_i; // Lower nearby array points
	double zwlo, zwhi, mwlo, mwhi, rwlo, rwhi, orwlo, orwhi;
	double lm = log10( m * unitconv::kgtoMsun );
	double lor = log10( offset_r * unitconv::mtokpc );
	double result = 0;

	if ( !loaded )
	{
		#pragma omp critical(tNFW_offset_sig_load)
		{
			if ( load( silent ) )
			{
				result = -1;
			}
		}
		if ( result == -1 )
			return result;
	}

	z_i = (int)( ( z - halo_z_min ) / halo_z_step );
	z_i = max( z_i, 0 );
	z_i = min( z_i, halo_z_res - 2 );

	zlo = halo_z_min + halo_z_step * z_i;
	zhi = halo_z_min + halo_z_step * ( z_i + 1 );
	zwlo = z - zlo;
	zwhi = zhi - z;

	m_i = (int)( ( lm - halo_m_min ) / halo_m_step );
	m_i = max( m_i, 0 );
	m_i = min( m_i, halo_m_res - 2 );

	mlo = halo_m_min + halo_m_step * m_i;
	mhi = halo_m_min + halo_m_step * ( m_i + 1 );
	mwlo = lm - mlo;
	mwhi = mhi - lm;

	r_i = (int)( ( r - r_min ) / r_step );
	r_i = max( r_i, 0 );
	r_i = min( r_i, r_res - 2 );

	rlo = r_min + r_step * r_i;
	rhi = r_min + r_step * ( r_i + 1 );
	rwlo = r - rlo;
	rwhi = rhi - r;

	or_i = (int)( ( lor - offset_r_min ) / offset_r_step );
	or_i = max( or_i, 0 );
	or_i = min( or_i, offset_r_res - 2 );

	orlo = offset_r_min + offset_r_step * r_i;
	orhi = offset_r_min + offset_r_step * ( r_i + 1 );
	orwlo = lor - orlo;
	orwhi = orhi - lor;

	double totweight = halo_z_step * halo_m_step * r_step * offset_r_step;

#ifdef _BRG_USE_UNITS_
	result = BRG_UNITS(0,-2,0,1,0,0,0);
#else
	result = 0;
#endif

	result += signal[z_i][m_i][r_i][or_i] * zwhi * mwhi * rwhi * orwhi;
	result += signal[z_i][m_i][r_i + 1][or_i] * zwhi * mwhi * rwlo * orwhi;
	result += signal[z_i][m_i + 1][r_i][or_i] * zwhi * mwlo * rwhi * orwhi;
	result += signal[z_i][m_i + 1][r_i + 1][or_i] * zwhi * mwlo * rwlo * orwhi;
	result += signal[z_i + 1][m_i][r_i][or_i] * zwlo * mwhi * rwhi * orwhi;
	result += signal[z_i + 1][m_i][r_i + 1][or_i] * zwlo * mwhi * rwlo * orwhi;
	result += signal[z_i + 1][m_i + 1][r_i][or_i] * zwlo * mwlo * rwhi * orwhi;
	result += signal[z_i + 1][m_i + 1][r_i + 1][or_i] * zwlo * mwlo * rwlo
			* orwhi;
	result += signal[z_i][m_i][r_i][or_i + 1] * zwhi * mwhi * rwhi * orwlo;
	result += signal[z_i][m_i][r_i + 1][or_i + 1] * zwhi * mwhi * rwlo * orwlo;
	result += signal[z_i][m_i + 1][r_i][or_i + 1] * zwhi * mwlo * rwhi * orwlo;
	result += signal[z_i][m_i + 1][r_i + 1][or_i + 1] * zwhi * mwlo * rwlo
			* orwlo;
	result += signal[z_i + 1][m_i][r_i][or_i + 1] * zwlo * mwhi * rwhi * orwlo;
	result += signal[z_i + 1][m_i][r_i + 1][or_i + 1] * zwlo * mwhi * rwlo
			* orwlo;
	result += signal[z_i + 1][m_i + 1][r_i][or_i + 1] * zwlo * mwlo * rwhi
			* orwlo;
	result += signal[z_i + 1][m_i + 1][r_i + 1][or_i + 1] * zwlo * mwlo * rwlo
			* orwlo;

	result /= totweight;

	return result;

}

#endif // end brgastro::NFW_offset_sig_cache functions

// brgastro::tNFW_group_sig_cache class methods
#if (1)

const int brgastro::tNFW_group_sig_cache::load( const bool silent )
{
	std::ifstream in_file;
	std::string file_data;
	bool need_to_calc = false;
	int loop_counter = 0;
	double temp_data;
	int i, j, k, l;

	if ( loaded )
		return 0;

	do
	{
		if ( loop_counter >= 2 )
		{
			if ( !silent )
				std::cerr
						<< "ERROR: infinite loop detected in brgastro::NFW_offset_sig_cache.\n";
			return INFINITE_LOOP_ERROR;
		}
		else
		{
			loop_counter++;
		}
		need_to_calc = false;

		if ( open_file( in_file, file_name ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Check that it has the right header
		getline( in_file, file_data );
		if ( file_data.compare( header_string ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Trim out any other commented lines
		if ( trim_comments_all_at_top( in_file ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Load range parameters;
		if ( !( in_file >> halo_z_min >> halo_z_max >> halo_z_step
				>> halo_m_min >> halo_m_max >> halo_m_step >> r_min >> r_max
				>> r_step >> group_c_min >> group_c_max >> group_c_step ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Set up data
		halo_z_res = (int)( ( halo_z_max - halo_z_min ) / halo_z_step ) + 1;
		halo_m_res = (int)( ( halo_m_max - halo_m_min ) / halo_m_step ) + 1;
		r_res = (int)( ( r_max - r_min ) / r_step ) + 1;
		group_c_res = (int)( ( group_c_max - group_c_min ) / group_c_step )
				+ 1;
		if ( make_array4d( signal, halo_z_res, halo_m_res, r_res,
				group_c_res ) )
			return 1;

		// Read in data
		i = j = k = l = 0;
		while ( !in_file.eof() )
		{
			in_file >> temp_data;
			signal[i][j][k][l] = temp_data;
			l++;
			if ( l >= group_c_res )
			{
				l = 0;
				k++;
				if ( k >= r_res )
				{
					k = 0;
					j++;
					if ( j >= halo_m_res )
					{
						j = 0;
						i++;
						if ( i >= halo_z_res )
							break;
					}
				}
			}
		}
		// Check that it was all read properly
		if ( i < halo_z_res )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

	} while ( need_to_calc );

	// Finish up
	in_file.close();
	in_file.clear();
	loaded = true;
	return 0;
}

const int brgastro::tNFW_group_sig_cache::unload()
{
	loaded = false;
	signal.clear();
	return 0;
}

const int brgastro::tNFW_group_sig_cache::calc( const bool silent )
{
	double z, lm, r, l_group_c;
	double mass, group_c;

	// Test that range is sane
	if ( ( halo_z_max <= halo_z_min ) || ( halo_z_step <= 0 )
			|| ( halo_m_max <= halo_m_min ) || ( halo_m_step <= 0 )
			|| ( r_max <= r_min ) || ( r_step <= 0 )
			|| ( group_c_max <= group_c_min ) || ( group_c_step <= 0 ) )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Bad range passed to tNFW_group_sig_cache::calc(silent)\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	// Set up data
	halo_z_res = (int)( ( halo_z_max - halo_z_min ) / halo_z_step ) + 1;
	halo_m_res = (int)( ( halo_m_max - halo_m_min ) / halo_z_step ) + 1;
	r_res = (int)( ( r_max - r_min ) / r_step ) + 1;
	group_c_res = (int)( ( group_c_max - group_c_min ) / group_c_step ) + 1;
	if ( make_array4d( signal, halo_z_res, halo_m_res, r_res, group_c_res ) )
		return 1;

	z = halo_z_min;
	for ( int i = 0; i < halo_z_res; i++ )
	{
		lm = halo_m_min;
		for ( int j = 0; j < halo_m_res; j++ )
		{
			mass = std::pow( 10, lm ) * unitconv::Msuntokg;
			r = r_min;
			for ( int k = 0; k < r_res; k++ )
			{
				l_group_c = group_c_min;
				for ( int l = 0; l < group_c_res; l++ )
				{
					group_c = std::pow( 10, l_group_c ) * unitconv::kpctom;
					signal[i][j][k][l] =
							brgastro::tNFW_profile( mass, z ).group_WLsig( r,
									group_c );
					l_group_c += group_c_step;
				}
				r += r_step;
			}
			lm += halo_m_step;
		}
		z += halo_z_step;
	}

	return 0;
}

const int brgastro::tNFW_group_sig_cache::output( const bool silent )
{
	std::ofstream out_file;
	std::string file_data;

	if ( !loaded )
	{
		if ( calc( silent ) )
			return 1;
	}

	if ( open_file( out_file, file_name ) )
		return 1;

	// Output header
	out_file << header_string << "\n#\n";

	// Output range
	out_file << halo_z_min << "\t" << halo_z_max << "\t" << halo_z_step << "\n"
			<< halo_m_min << "\t" << halo_m_max << "\t" << halo_m_step << "\n"
			<< r_min << "\t" << r_max << "\t" << r_step << "\n" << group_c_min
			<< "\t" << group_c_max << "\t" << group_c_step << "\n";

	// Output data
	for ( int i = 0; i < halo_z_res; i++ )
	{
		for ( int j = 0; j < halo_m_res; j++ )
		{
			for ( int k = 0; k < r_res; k++ )
			{
				for ( int l = 0; l < group_c_res; l++ )
				{
					if ( !( out_file << signal[i][j][k][l] << "\n" ) )
						return FILE_ACCESS_ERROR;
				}
			}
		}
	}

	out_file.close();
	out_file.clear();

	return 0;

}

const int brgastro::tNFW_group_sig_cache::set_file_name(
		const std::string new_name )
{
	file_name = new_name;
	if ( loaded )
	{
		if ( unload() )
			return errorNOS();
	}
	return 0;
}

const int brgastro::tNFW_group_sig_cache::set_range(
		const double new_halo_z_min, const double new_halo_z_max,
		const double new_halo_z_step, const double new_halo_m_min,
		const double new_halo_m_max, const double new_halo_m_step,
		const double new_r_min, const double new_r_max,
		const double new_r_step, const double new_group_c_min,
		const double new_group_c_max, const double new_group_c_step,
		const bool silent )
{
	// First we try to load
	if ( !loaded )
		load( silent );

	// Go through variables, check if any are actually changed. If so, recalculate cache
	if ( ( halo_z_min != new_halo_z_min ) || ( halo_z_max != new_halo_z_max )
			|| ( halo_z_step != new_halo_z_step )
			|| ( halo_m_min != new_halo_m_min )
			|| ( halo_m_max != new_halo_m_max )
			|| ( halo_m_step != new_halo_m_step )
			|| ( halo_z_min != new_r_min ) || ( halo_z_max != new_r_max )
			|| ( halo_z_step != new_r_step )
			|| ( halo_z_min != new_group_c_min )
			|| ( halo_z_max != new_group_c_max )
			|| ( halo_z_step != new_group_c_step ) )
	{
		halo_z_min = new_halo_z_min;
		halo_z_max = new_halo_z_max;
		halo_z_step = new_halo_z_step;
		halo_m_min = new_halo_m_min;
		halo_m_max = new_halo_m_max;
		halo_m_step = new_halo_m_step;
		r_min = new_r_min;
		r_max = new_r_max;
		r_step = new_r_step;
		group_c_min = new_group_c_min;
		group_c_max = new_group_c_max;
		group_c_step = new_group_c_step;

		if ( unload() )
			return errorNOS( silent );
		if ( calc( silent ) )
			return 1;
	}
	return 0;
}

const int brgastro::tNFW_group_sig_cache::set_precision(
		const int new_precision, const bool silent )
{
	if ( new_precision > 0 )
	{
		sig_digits = min( new_precision, DBL_MAX_PRECISION );
		return 0;
	}
	else
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Precision for tNFW_group_sig_cache must be > 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
}

const BRG_UNITS brgastro::tNFW_group_sig_cache::get( const double z,
		const BRG_MASS m, const BRG_DISTANCE r, const double group_c,
		const bool silent )
{
	double zlo, zhi, mlo, mhi, rlo, rhi, gclo, gchi;
	int z_i, m_i, r_i, gc_i; // Lower nearby array points
	double zwlo, zwhi, mwlo, mwhi, rwlo, rwhi, gcwlo, gcwhi;
	double lm = log10( m * unitconv::kgtoMsun );
#ifdef _BRG_USE_UNITS_
	BRG_UNITS result(0,-2,0,1,0,0,0);
#else
	double result = 0;
#endif

	if ( !loaded )
	{
		#pragma omp critical(tNFW_group_sig_load)
		{
			if ( load( silent ) )
			{
				result = -1;
			}
		}
		if ( result == -1 )
			return result;
	}

	z_i = (int)( ( z - halo_z_min ) / halo_z_step );
	z_i = max( z_i, 0 );
	z_i = min( z_i, halo_z_res - 2 );

	zlo = halo_z_min + halo_z_step * z_i;
	zhi = halo_z_min + halo_z_step * ( z_i + 1 );
	zwlo = z - zlo;
	zwhi = zhi - z;

	m_i = (int)( ( lm - halo_m_min ) / halo_m_step );
	m_i = max( m_i, 0 );
	m_i = min( m_i, halo_m_res - 2 );

	mlo = halo_m_min + halo_m_step * m_i;
	mhi = halo_m_min + halo_m_step * ( m_i + 1 );
	mwlo = m - mlo;
	mwhi = mhi - m;

	r_i = (int)( ( r - r_min ) / r_step );
	r_i = max( r_i, 0 );
	r_i = min( r_i, r_res - 2 );

	rlo = r_min + r_step * r_i;
	rhi = r_min + r_step * ( r_i + 1 );
	rwlo = r - rlo;
	rwhi = rhi - r;

	gc_i = (int)( ( group_c - group_c_min ) / group_c_step );
	gc_i = max( gc_i, 0 );
	gc_i = min( gc_i, group_c_res - 2 );

	gclo = group_c_min + group_c_step * r_i;
	gchi = group_c_min + group_c_step * ( r_i + 1 );
	gcwlo = group_c - gclo;
	gcwhi = gchi - group_c;

	double totweight = halo_z_step * halo_m_step * r_step * group_c_step;

	result = 0;

	result += signal[z_i][m_i][r_i][gc_i] * zwhi * mwhi * rwhi * gcwhi;
	result += signal[z_i][m_i][r_i + 1][gc_i] * zwhi * mwhi * rwlo * gcwhi;
	result += signal[z_i][m_i + 1][r_i][gc_i] * zwhi * mwlo * rwhi * gcwhi;
	result += signal[z_i][m_i + 1][r_i + 1][gc_i] * zwhi * mwlo * rwlo * gcwhi;
	result += signal[z_i + 1][m_i][r_i][gc_i] * zwlo * mwhi * rwhi * gcwhi;
	result += signal[z_i + 1][m_i][r_i + 1][gc_i] * zwlo * mwhi * rwlo * gcwhi;
	result += signal[z_i + 1][m_i + 1][r_i][gc_i] * zwlo * mwlo * rwhi * gcwhi;
	result += signal[z_i + 1][m_i + 1][r_i + 1][gc_i] * zwlo * mwlo * rwlo
			* gcwhi;
	result += signal[z_i][m_i][r_i][gc_i + 1] * zwhi * mwhi * rwhi * gcwlo;
	result += signal[z_i][m_i][r_i + 1][gc_i + 1] * zwhi * mwhi * rwlo * gcwlo;
	result += signal[z_i][m_i + 1][r_i][gc_i + 1] * zwhi * mwlo * rwhi * gcwlo;
	result += signal[z_i][m_i + 1][r_i + 1][gc_i + 1] * zwhi * mwlo * rwlo
			* gcwlo;
	result += signal[z_i + 1][m_i][r_i][gc_i + 1] * zwlo * mwhi * rwhi * gcwlo;
	result += signal[z_i + 1][m_i][r_i + 1][gc_i + 1] * zwlo * mwhi * rwlo
			* gcwlo;
	result += signal[z_i + 1][m_i + 1][r_i][gc_i + 1] * zwlo * mwlo * rwhi
			* gcwlo;
	result += signal[z_i + 1][m_i + 1][r_i + 1][gc_i + 1] * zwlo * mwlo * rwlo
			* gcwlo;

	result /= totweight;

	return result;

}

#endif // end brgastro::tNFW_group_sig_cache functions


