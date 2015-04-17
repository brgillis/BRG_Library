/*
 * common.h
 *
 *  Created on: 17 Apr 2015
 *      Author: brg
 */

#ifndef BRG_LENSING_COMMON_H_
#define BRG_LENSING_COMMON_H_

#include <vector>

#include "Eigen/Core"

typedef double flt_type;
typedef long double long_flt_type;

typedef int int_type;
typedef long int long_int_type;

typedef Eigen::Array<flt_type,Eigen::Dynamic,1> flt_array_type;
typedef std::vector<flt_type> flt_vector_type;

typedef Eigen::Array<int_type,Eigen::Dynamic,1> int_array_type;
typedef std::vector<int_type> int_vector_type;

#endif /* BRG_LENSING_COMMON_H_ */
