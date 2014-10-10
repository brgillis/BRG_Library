/**********************************************************************\
 @file sgsmooth.h
 ------------------

 TODO <Insert file description here>

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


// body file: sgsmooth.cpp


#ifndef _BRG_EXTERNAL_SGSMOOTH_H_INCLUDED_
#define _BRG_EXTERNAL_SGSMOOTH_H_INCLUDED_

#include <vector>

namespace brgastro{ namespace sgsmooth{

/* common definitions for the sgsmooth plugin */

std::vector<double> sg_smooth(const std::vector<double> & input,
                             const size_t window, const size_t order);

std::vector<double> sg_derivative(const std::vector<double> & input,
                             const size_t window, const size_t order,
                             const double delta=1.0);

void sgs_error(const char *errmsg);

} // end namespace brgastro::sgsmooth

using sgsmooth::sg_smooth;
using sgsmooth::sg_derivative;

}// end all namespaces

#endif // _BRG_EXTERNAL_SGSMOOTH_H_INCLUDED_
