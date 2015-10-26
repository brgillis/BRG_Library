/*
 * Copyright (C) 2012-2020 Bryan R. Gillis
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
 * @file fitsfile_deleter.hpp
 *
 * @date Oct 22, 2015
 * @author brg
 */

#ifndef ICEBRG_MAIN_FILE_ACCESS_FITS_FITSFILE_DELETER_H_
#define ICEBRG_MAIN_FILE_ACCESS_FITS_FITSFILE_DELETER_H_

namespace IceBRG {

struct fitsfile_deleter
{
  void operator()(void *p);
};

} /* namespace IceBRG */
#endif /* ifndef ICEBRG_MAIN_FILE_ACCESS_FITS_FITSFILE_DELETER_H_ */
