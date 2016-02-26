/**********************************************************************\
 @file IceBRGpy.i
 ------------------

 TODO <Insert file description here>

 **********************************************************************

 Copyright (C) 2015 brg

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

// SWIG includes
%include "typemaps.i"
%include "std_string.i"
%include "std_vector.i"

%module IceBRGpy

%{
	 
	/* Include the headers in the wrapper code */
	#include "IceBRG_main/units/unit_conversions.hpp"

	#include "IceBRG_physics/abundance_matching.hpp"
	#include "IceBRG_physics/cluster_visibility.hpp"
	#include "IceBRG_physics/constants.hpp"
	#include "IceBRG_physics/cosmology.hpp"
	#include "IceBRG_physics/distance_measures.hpp"
	#include "IceBRG_physics/galaxy_visibility.hpp"
	#include "IceBRG_physics/luminosity.hpp"
	#include "IceBRG_physics/mass_function.hpp"

	#include "IceBRG_physics/detail/redshift_obj.hpp"
	#include "IceBRG_physics/density_profile/detail/density_profile.hpp"
	#include "IceBRG_physics/density_profile/point_mass_profile.hpp"
	#include "IceBRG_physics/density_profile/tNFW_profile.hpp"
	#include "IceBRG_physics/sky_obj/detail/sky_obj.hpp"
	#include "IceBRG_physics/sky_obj/galaxy.hpp"

	#include "IceBRG_lensing/detail/lensing_profile_extension.hpp"
	#include "IceBRG_lensing/detail/pair_bin_summary.hpp"
	#include "IceBRG_lensing/detail/pair_bin.hpp"
	#include "IceBRG_lensing/detail/pair_bins_summary.hpp"
	#include "IceBRG_lensing/detail/source_obj.hpp"
	#include "IceBRG_lensing/lens_source_pair.hpp"
	#include "IceBRG_lensing/lensing_tNFW_profile.hpp"
	#include "IceBRG_lensing/pair_binner.hpp"
	#include "IceBRG_lensing/source_galaxy.hpp"
	
	using namespace IceBRG;
	 
%}
 
// Parse the header files to generate wrappers
%include "/disk2/brg/include/IceBRG_main/units/unit_conversions.hpp"

%include "/disk2/brg/include/IceBRG_physics/abundance_matching.hpp"
%include "/disk2/brg/include/IceBRG_physics/cluster_visibility.hpp"
%include "/disk2/brg/include/IceBRG_physics/constants.hpp"
%include "/disk2/brg/include/IceBRG_physics/cosmology.hpp"
%include "/disk2/brg/include/IceBRG_physics/distance_measures.hpp"
%include "/disk2/brg/include/IceBRG_physics/galaxy_visibility.hpp"
%include "/disk2/brg/include/IceBRG_physics/luminosity.hpp"
%include "/disk2/brg/include/IceBRG_physics/mass_function.hpp"

%include "/disk2/brg/include/IceBRG_physics/detail/redshift_obj.hpp"
%include "/disk2/brg/include/IceBRG_physics/density_profile/detail/density_profile.hpp"
%include "/disk2/brg/include/IceBRG_physics/density_profile/point_mass_profile.hpp"
%include "/disk2/brg/include/IceBRG_physics/density_profile/tNFW_profile.hpp"
%include "/disk2/brg/include/IceBRG_physics/sky_obj/detail/sky_obj.hpp"
%include "/disk2/brg/include/IceBRG_physics/sky_obj/galaxy.hpp"

%include "/disk2/brg/include/IceBRG_lensing/detail/lensing_profile_extension.hpp"
%include "/disk2/brg/include/IceBRG_lensing/detail/pair_bin_summary.hpp"
%include "/disk2/brg/include/IceBRG_lensing/detail/pair_bin.hpp"
%include "/disk2/brg/include/IceBRG_lensing/detail/pair_bins_summary.hpp"
%include "/disk2/brg/include/IceBRG_lensing/detail/source_obj.hpp"
%include "/disk2/brg/include/IceBRG_lensing/lens_source_pair.hpp"
%template(lensing_tNFW_profile_extension) IceBRG::lensing_profile_extension< IceBRG::lensing_tNFW_profile >;
%include "/disk2/brg/include/IceBRG_lensing/lensing_tNFW_profile.hpp"
%include "/disk2/brg/include/IceBRG_lensing/pair_binner.hpp"
%include "/disk2/brg/include/IceBRG_lensing/source_galaxy.hpp"

// Tell Swig about typedefs in use
typedef double flt_t;

typedef flt_t dimensionless_type;

typedef flt_t distance_type;
typedef flt_t area_type;
typedef flt_t volume_type;
typedef flt_t inverse_distance_type;
typedef flt_t inverse_area_type;
typedef flt_t inverse_volume_type;

typedef flt_t time_type;
typedef flt_t inverse_time_type;

typedef flt_t mass_type;

typedef flt_t angle_type;
typedef flt_t square_angle_type;
typedef flt_t inverse_angle_type;
typedef flt_t inverse_square_angle_type;

typedef flt_t temperature_type;

typedef flt_t velocity_type;
typedef flt_t acceleration_type;

typedef flt_t density_type;
typedef flt_t inverse_density_type;
typedef flt_t surface_density_type;
typedef flt_t inverse_surface_density_type;

typedef flt_t inverse_volume_inverse_mass_type;

typedef flt_t any_units_type;