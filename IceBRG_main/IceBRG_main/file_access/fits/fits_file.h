/*  
 * Copyright (C) 2015 Bryan R. Gillis
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
 * @file fits_file.h
 *
 * @date Oct 22, 2015
 * @author user
 */

#ifndef ICEBRG_MAIN_FILE_ACCESS_FITS_FITS_FILE_HPP_
#define ICEBRG_MAIN_FILE_ACCESS_FITS_FITS_FILE_HPP

#include <memory>

#include "IceBRG_main/common.h"
#include "IceBRG_main/file_access/fits/fitsfile_deleter.h"

// Forward declare the fitsfile struct
namespace fitsio {
struct fitsfile;
}

namespace IceBRG {


class fits_file {

  typedef std::unique_ptr<fitsio::fitsfile,fitsfile_deleter>
    p_fitsfile_type;

  p_fitsfile_type _p_fits_file;

public:

  /* Default Constructor */
  fits_file ()
  : _p_fits_file(nullptr)
  { }

  // Construct from filename
  fits_file (const str_type & filename )
  : fits_file()
  { open_readonly(filename); }

  /* Deconstructor */
  virtual ~fits_file () {}

  // Open/close a file
  void open( const str_type & filename, const int & mode = 1 );
  void open_readonly( const str_type & filename );
  void open_readwrite( const str_type & filename );
  void close();

private:

};

} /* namespace IceBRG */
#endif /* ifndef ICEBRG_MAIN_FILE_ACCESS_FITS_FITS_FILE_HPP_ */
