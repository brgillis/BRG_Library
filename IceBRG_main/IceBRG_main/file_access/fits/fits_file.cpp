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

extern "C" {
#include <cfitsio/fitsio.h>
}

#include "IceBRG_main/common.h"
#include "IceBRG_main/file_access/fits/fits_err.h"
#include "IceBRG_main/file_access/fits/fits_file.h"

namespace IceBRG {

// Impl struct
struct fits_file::fits_file_impl {

	typedef std::unique_ptr<fitsfile,fitsfile_deleter> p_fitsfile_type;

	p_fitsfile_type p_fitsfile;

	// Constructor
	fits_file_impl(fitsfile *p_fitsfile = nullptr)
	: p_fitsfile(p_fitsfile)
	{
	}

};

/* Default Constructor */
fits_file::fits_file ()
{
}

// Construct from filename
fits_file::fits_file (const str_type & filename )
: fits_file()
{ open_readonly(filename); }

/* Default Destructor */
fits_file::~fits_file ()
{
}

// Open/close a file
void fits_file::open( const str_type & filename, const int & mode )
{
  int status;
  fits_file_impl::p_fitsfile_type new_file(new fitsfile);

  auto p_new_file = new_file.get();

  fits_open_file(&p_new_file,
                 filename.c_str(),
                 mode,
                 &status);

  if(status) fitsio_error(status);

  std::swap(_p_impl->p_fitsfile,new_file);

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
	_p_impl->p_fitsfile.release();
}

} // namespace IceBRG



