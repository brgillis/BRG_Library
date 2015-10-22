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
 * @file fits_file.cpp
 *
 * @date Oct 22, 2015
 * @author brg
 */

namespace fitsio {
extern "C" {
#include <cfitsio/fitsio.h>
}
} // namespace fitsio

#include "IceBRG_main/common.h"
#include "IceBRG_main/file_access/fits/fits_err.h"

namespace IceBRG {

// Open/close a file
void fits_file::open( const str_type & filename, const int & mode = 1 )
{
  int status;
  p_fitsfile_type new_file(new fitsio::fitsfile);

  fitsio::fits_open_file(&(new_file.get()),
                            filename.c_str(),
                            mode,
                            &status);

  if(status) fitsio_error(status);

  return;
}
void fits_file::open_readonly( const str_type & filename )
{
  open(filename,READONLY);
}
void fits_file::open_readwrite( const str_type & filename )
{
  open(filename,READWRITE);
}

void fits_file::close()
{
  _p_fits_file.release();
}

} // namespace IceBRG



