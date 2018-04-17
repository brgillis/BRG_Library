/**********************************************************************\
 @file manipulations.hpp
 ------------------

 Functions to manipulate vectors in various ways.

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

#ifndef _BRG_MANIPULATIONS_HPP_INCLUDED_
#define _BRG_MANIPULATIONS_HPP_INCLUDED_

#include <vector>

#include "IceBRG_main/common.hpp"

#include "IceBRG_main/math/misc_math.hpp"
#include "IceBRG_main/utility.hpp"
#include "IceBRG_main/vector/make_vector.hpp"

namespace IceBRG {

template<typename T1>
void concatenate(T1 & v1, const T1 & v2)
{
	v1.reserve(v1.size()+v2.size());
	v1.insert(v1.end(),v2.begin(),v2.end());
}

template<typename T1, typename T2>
void concatenate(T1 & v1, const T2 & v2)
{
	v1.reserve(v1.size()+v2.size());
	for(auto it=v2.begin(); it!=v2.end(); ++it)
	{
		v1.push_back(*it);
	}
}

template<typename T1>
void concatenate_map(T1 & v1, const T1 & v2)
{
	v1.insert(v2.begin(),v2.end());
}

template<typename T>
std::vector< std::vector<T> > pad(const std::vector< std::vector<T> > & v, const T & default_value=T())
{
	size_t n_cols = v.size();
	size_t n_rows = 0;
	for(size_t i=0; i<n_cols; ++i)
		if(v[i].size() > n_rows) n_rows = v[i].size();

	std::vector< std::vector<T> > result = v;

	for(size_t i=0; i<n_cols; ++i)
		if(result[i].size() < n_rows) result[i].resize(n_rows,default_value);

	return result;
}

template<typename T>
std::vector< std::vector<T> > transpose(const std::vector< std::vector<T> > & v, const T & default_value=T())
{
	size_t n_cols = v.size();
	size_t n_rows = 0;

	size_t num_col_lengths=0;

	for(size_t i=0; i<v.size(); ++i)
	{
		if(v[i].size() != n_rows)
		{
			++num_col_lengths;
			if(v[i].size() > n_rows)
				n_rows = v[i].size();
		}
	}

	std::vector< std::vector<T> > cv;
	const std::vector< std::vector<T> > *pv = &v;

	if(num_col_lengths>1)
	{
		cv = pad(v,default_value);
		pv = &cv;
	}

	std::vector< std::vector<T> > result;

	auto result_function = [&] (size_t i, size_t j)
		{
			return (*pv).at(j).at(i);
		};

	make_vector_function(result,result_function,n_rows,n_cols);

	return result;
}

template<typename T>
std::vector< std::vector<T> > reverse_vertical(const std::vector< std::vector<T> > & v)
{
	size_t n_cols = v.size();
	size_t n_rows = 0;
	for(size_t i=0; i<v.size(); ++i)
		if(v[i].size() > n_rows) n_rows = v[i].size();

	std::vector< std::vector<T> > result;

	auto result_function = [&] (size_t i, size_t j)
		{
			return v[i][n_rows-j-1];
		};

	make_vector_function(result,result_function,n_cols,n_rows);

	return result;
}

template<typename T>
std::vector<T> remove_outliers( std::vector<T> const & v_in, int const & min_remaining_members = 2 )
{

    assert(min_remaining_members>=2);

    int num_tot = v_in.size();
    int num_masked = 0;
    bool outliers_found = true;

    std::vector<bool> mask(v_in.size(),false);


    while(outliers_found and ((num_tot-num_masked)>min_remaining_members))
    {

        outliers_found = false;

        // Find any outliers and remove them

        // Get the mean and sigma of the unmasked elements
        double mean = 0;
        for(int i=0; i<num_tot; ++i)
        {
            if(!mask[i])
                mean += v_in[i];
        }
        mean /= (num_tot-num_masked);

        double sigma = 0;
        for(int i=0; i<num_tot; ++i)
        {
            if(!mask[i])
                sigma += square(v_in[i]-mean);
        }
        sigma /= (num_tot-num_masked-1);
        sigma = std::sqrt(sigma);

        for(int i=0; i<num_tot; ++i)
        {
            if(mask[i])
                continue;

            double prob = 1-0.5 * std::erfc(-(std::abs(v_in[i] - mean) / sigma) * M_SQRT1_2);

            if(prob*num_tot<0.5)
            {
                mask[i] = true;
                outliers_found = true;
                ++num_masked;
            }
        }
    }

    // Did we remove any at all?
    if(num_masked==0)
        return v_in;

    // Did we remove too many?
    if(num_tot-num_masked<min_remaining_members)
    {
        throw std::runtime_error("Too few elements left after removing outliers.");
    }

    std::vector<T> v_out;
    v_out.reserve(num_tot);
    for(int i=0; i<num_tot; ++i)
    {
        if(!mask[i])
        {
            v_out.push_back(v_in[i]);
        }
    }

    return v_out;
}

} // namespace IceBRG

#endif // _BRG_MANIPULATIONS_HPP_INCLUDED_
