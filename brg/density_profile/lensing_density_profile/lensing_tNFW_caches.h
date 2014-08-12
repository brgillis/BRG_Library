/**       @file tNFW_caches.h
 *
 *     Project: brg
 *        Path: /brg/density_profile/tNFW_caches.h
 *
 *  Created on: 11 Aug 2014
 *      Author: brg
 */

#ifndef _BRG_TNFW_CACHES_H_
#define _BRG_TNFW_CACHES_H_

#include <string>
#include <vector>

#include "../../brg_global.h"

namespace brgastro {

class tNFW_sig_cache
{
	// Weak lensing from tNFW profile cache
	static bool loaded;
	static std::string file_name;
	static double halo_z_min, halo_z_max, halo_z_step;
	static int halo_z_res;
	static double r_min, r_max, r_step;
	static int r_res;
	static double halo_m_min, halo_m_max, halo_m_step;
	static int halo_m_res;
	static std::vector< std::vector< std::vector< double > > > signal;
	static std::string header_string;
	static int sig_digits;
	const int load( const bool silent = false );
	const int unload();
	const int calc( const bool silent = false );
	const int output( const bool silent = false );
public:
	const int set_file_name( const std::string new_name );
	const int set_range( const double new_z_halo_min,
			const double new_z_halo_max, const double new_z_halo_step,
			const double new_m_halo_min, const double new_m_halo_max,
			const double new_m_halo_step, const double new_r_min,
			const double new_r_max, const double new_r_step,
			const bool silent = false );
	const int set_precision( const int new_precision,
			const bool silent = false );

	const BRG_UNITS get( const double z_halo, const BRG_MASS m_halo,
			const BRG_DISTANCE r_halo, const bool silent = false ); // Takes in standard units

};
// class tNFW_sig_cache

class tNFW_offset_sig_cache
{
	// Offset weak lensing signal from tNFW profile cache
	static bool loaded;
	static std::string file_name;
	static double halo_z_min, halo_z_max, halo_z_step;
	static int halo_z_res;
	static double r_min, r_max, r_step;
	static int r_res;
	static double halo_m_min, halo_m_max, halo_m_step;
	static int halo_m_res;
	static double offset_r_min, offset_r_max, offset_r_step;
	static int offset_r_res;
	static std::vector< std::vector< std::vector< std::vector< double > > > > signal;
	static std::string header_string;
	static int sig_digits;
	const int load( const bool silent = false );
	const int unload();
	const int calc( const bool silent = false );
	const int output( const bool silent = false );
public:
	const int set_file_name( const std::string new_name );
	const int set_range( const double new_z_halo_min,
			const double new_z_halo_max, const double new_z_halo_step,
			const double new_m_halo_min, const double new_m_halo_max,
			const double new_m_halo_step, const double new_r_min,
			const double new_r_max, const double new_r_step,
			const double new_offset_r_min, const double new_offset_r_max,
			const double new_offset_r_step, const bool silent = false );
	const int set_precision( const int new_precision,
			const bool silent = false );

	const BRG_UNITS get( const double z_halo, const BRG_MASS m_halo,
			const BRG_DISTANCE r, const BRG_DISTANCE offset_r,
			const bool silent = false ); // Takes in standard units

};
// class tNFW_offset_sig_cache

class tNFW_group_sig_cache
{
	// Group weak lensing signal from tNFW profile cache
	static bool loaded;
	static std::string file_name;
	static double halo_z_min, halo_z_max, halo_z_step;
	static int halo_z_res;
	static double r_min, r_max, r_step;
	static int r_res;
	static double halo_m_min, halo_m_max, halo_m_step;
	static int halo_m_res;
	static double group_c_min, group_c_max, group_c_step;
	static int group_c_res;
	static std::vector< std::vector< std::vector< std::vector< double > > > > signal;
	static std::string header_string;
	static int sig_digits;
	const int load( const bool silent = false );
	const int unload();
	const int calc( const bool silent = false );
	const int output( const bool silent = false );
public:
	const int set_file_name( const std::string new_name );
	const int set_range( const double new_z_halo_min,
			const double new_z_halo_max, const double new_z_halo_step,
			const double new_m_halo_min, const double new_m_halo_max,
			const double new_m_halo_step, const double new_r_min,
			const double new_r_max, const double new_r_step,
			const double new_group_c_min, const double new_group_c_max,
			const double new_group_c_step, const bool silent = false );
	const int set_precision( const int new_precision,
			const bool silent = false );

	const BRG_UNITS get( const double z_halo, const BRG_MASS m_halo,
			const BRG_DISTANCE r, const double group_c, const bool silent =
					false ); // Takes in standard units

};
// class tNFW_group_sig_cache

} // end namespace brgastro

#endif /* _BRG_TNFW_CACHES_H_ */
