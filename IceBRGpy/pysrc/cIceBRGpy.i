/**********************************************************************\
 @file cIceBRGpy.i
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
%include "std_pair.i"
%include "std_string.i"
%include "std_vector.i"

%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
import_array();
%}
// some additional backward compatibility declarations for supporting numpy < 1.7.0
%{
#if NPY_API_VERSION < 0x00000007
#define NPY_ARRAY_C_CONTIGUOUS NPY_C_CONTIGUOUS
#define NPY_ARRAY_F_CONTIGUOUS NPY_F_CONTIGUOUS
#define NPY_ARRAY_ALIGNED  NPY_ALIGNED
#endif
%}

%module cIceBRGpy

%{
	 
	/* Include the headers in the wrapper code */
	
	#include <utility>

	#include "SHE_GST_IceBRG_main/globals.hpp"
	#include "SHE_GST_IceBRG_main/units/unit_conversions.hpp"
	#include "SHE_GST_IceBRG_main/vector/rebin.hpp"

	#include "SHE_GST_IceBRG_physics/abundance_matching.hpp"
	#include "SHE_GST_IceBRG_physics/cluster_visibility.hpp"
	#include "SHE_GST_IceBRG_physics/constants.hpp"
	#include "SHE_GST_IceBRG_physics/cosmology.hpp"
	#include "SHE_GST_IceBRG_physics/distance_measures.hpp"
	#include "SHE_GST_IceBRG_physics/galaxy_visibility.hpp"
	#include "SHE_GST_IceBRG_physics/luminosity.hpp"
	#include "SHE_GST_IceBRG_physics/mass_function.hpp"

	#include "SHE_GST_IceBRG_physics/detail/redshift_obj.hpp"
	#include "SHE_GST_IceBRG_physics/density_profile/detail/density_profile.hpp"
	#include "SHE_GST_IceBRG_physics/density_profile/point_mass_profile.hpp"
	#include "SHE_GST_IceBRG_physics/density_profile/tNFW_profile.hpp"
	#include "SHE_GST_IceBRG_physics/sky_obj/detail/sky_obj.hpp"
	#include "SHE_GST_IceBRG_physics/sky_obj/galaxy.hpp"

	#include "SHE_GST_IceBRG_lensing/detail/lensing_profile_extension.hpp"
	#include "SHE_GST_IceBRG_lensing/detail/pair_bin_summary.hpp"
	#include "SHE_GST_IceBRG_lensing/detail/pair_bin.hpp"
	#include "SHE_GST_IceBRG_lensing/detail/pair_bins_summary.hpp"
	#include "SHE_GST_IceBRG_lensing/detail/source_obj.hpp"
	#include "SHE_GST_IceBRG_lensing/lens_source_pair.hpp"
	#include "SHE_GST_IceBRG_lensing/lensing_tNFW_profile.hpp"
	#include "SHE_GST_IceBRG_lensing/pair_binner.hpp"
	#include "SHE_GST_IceBRG_lensing/source_galaxy.hpp"
	
	using namespace IceBRG;

	template< typename T >
	std::pair<int,int> rebin_wrap( T * p_image,
			int ss_nx,
			int ss_ny,
			int x_offset=0,
			int y_offset=0,
			int subsampling_factor=5)
	{
		auto out_array = IceBRG::rebin(p_image,
							ss_nx,
							ss_ny,
							x_offset,
							y_offset,
							subsampling_factor);
		
		for(int i=0; i<ss_nx*ss_ny; ++i)
		{
			p_image[i] = 0;
		}
		
		int ncols = out_array.cols();
		int nrows = out_array.rows();

		for(int i=0; i<ncols; ++i)
		{
			for(int j=0; j<nrows; ++j)
			{
				p_image[j+i*nrows] = out_array(j,i);
			}
		}
		
		return std::make_pair(nrows,ncols);
	}
	 
%}
 
// Parse the header files to generate wrappers
%include "/home/brg/include/SHE_GST_IceBRG_main/units/unit_conversions.hpp"

%apply (int* INPLACE_ARRAY2, int DIM1, int DIM2)
	{( int * p_image,
			int ss_nx,
			int ss_ny)}
%apply (long* INPLACE_ARRAY2, int DIM1, int DIM2)
	{( long* p_image,
			int ss_nx,
			int ss_ny)}
%apply (unsigned int* INPLACE_ARRAY2, int DIM1, int DIM2)
	{( unsigned int * p_image,
			int ss_nx,
			int ss_ny)}
%apply (unsigned long* INPLACE_ARRAY2, int DIM1, int DIM2)
	{( unsigned long* p_image,
			int ss_nx,
			int ss_ny)}
%apply (float* INPLACE_ARRAY2, int DIM1, int DIM2)
	{( float * p_image,
			int ss_nx,
			int ss_ny)}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2)
	{( double * p_image,
			int ss_nx,
			int ss_ny)}

%template() std::pair<int,int>;

%include "/home/brg/include/SHE_GST_IceBRG_main/vector/rebin.hpp"


template< typename T >
std::pair<int,int>
rebin_wrap( T * p_image,
		int ss_nx,
		int ss_ny,
		int x_offset=0,
		int y_offset=0,
		int subsampling_factor=5 );

%template(rebin_int) rebin_wrap<int>;
%template(rebin_long) rebin_wrap<long>;
%template(rebin_uint) rebin_wrap<unsigned int>;
%template(rebin_ulong) rebin_wrap<unsigned long>;
%template(rebin_float) rebin_wrap<float>;
%template(rebin_double) rebin_wrap<double>;

%inline {
	IceBRG::str_t get_workdir() { return IceBRG::globals::workdir; }
	void set_workdir( IceBRG::str_t const & workdir ) { IceBRG::globals::workdir = workdir; }
	
	error_behavior_type get_error_behavior() { return IceBRG::globals::error_behavior; }
	void set_error_behavior( error_behavior_type const & error_behavior )
	{ IceBRG::globals::error_behavior = error_behavior; }
}

%include "/home/brg/include/SHE_GST_IceBRG_physics/abundance_matching.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_physics/cluster_visibility.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_physics/constants.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_physics/cosmology.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_physics/distance_measures.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_physics/galaxy_visibility.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_physics/luminosity.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_physics/mass_function.hpp"

%include "/home/brg/include/SHE_GST_IceBRG_physics/detail/redshift_obj.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_physics/density_profile/detail/density_profile.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_physics/density_profile/point_mass_profile.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_physics/density_profile/tNFW_profile.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_physics/sky_obj/detail/sky_obj.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_physics/sky_obj/galaxy.hpp"

%include "/home/brg/include/SHE_GST_IceBRG_lensing/detail/lensing_profile_extension.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_lensing/detail/pair_bin_summary.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_lensing/detail/pair_bin.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_lensing/detail/pair_bins_summary.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_lensing/detail/source_obj.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_lensing/lens_source_pair.hpp"

%template(lensing_tNFW_profile_extension)
	IceBRG::lensing_profile_extension< IceBRG::lensing_tNFW_profile >;
%include "/home/brg/include/SHE_GST_IceBRG_lensing/lensing_tNFW_profile.hpp"

%include "/home/brg/include/SHE_GST_IceBRG_lensing/pair_binner.hpp"
%include "/home/brg/include/SHE_GST_IceBRG_lensing/source_galaxy.hpp"

// Tell Swig about typedefs in use
typedef double flt_t;
typedef std::string str_t;

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
