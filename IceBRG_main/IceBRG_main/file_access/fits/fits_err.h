/*
 * Copyright (C) 2012-2020 Euclid Science Ground Segment
 *
 * This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either version 3.0 of the License, or (at your option)
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/**
 * @file fits_err.h
 *
 * @date Oct 22, 2015
 * @author user
 */


#ifndef ICEBRG_MAIN_FILE_ACCESS_FITS_FITS_ERR_H_
#define ICEBRG_MAIN_FILE_ACCESS_FITS_FITS_ERR_H_

#include "IceBRG_main/common.h"

namespace IceBRG {

void fitsio_error(const int_type & status, const bool no_except=false);

} // namespace IceBRG

#endif /* ifndef ICEBRG_MAIN_FILE_ACCESS_FITS_FITS_ERR_H_ */
