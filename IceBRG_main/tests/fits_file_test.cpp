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
 * @file fits_file_test.cpp
 *
 * @date Oct 22, 2015
 * @author user
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "IceBRG_main/common.h"
#include "IceBRG_main/file_access/fits/fits_file.h"

namespace IceBRG {

struct fits_file_fixture {

  const str_type aux_path = "auxdir/tests/fits";
  const str_type good_filename = "empty.fits";
  const str_type bad_filename = "invalid.fits";
  const str_type no_filename = "nonexistant.fits";

};

BOOST_AUTO_TEST_SUITE (FITS_FILE_TEST)

BOOST_FIXTURE_TEST_CASE(test_fits_file, fits_file_fixture) {

  BOOST_CHECK_NO_THROW(fits_file(good_filename));

  BOOST_CHECK_THROW(fits_file(bad_filename),std::runtime_error);

  BOOST_CHECK_THROW(fits_file(no_filename),std::runtime_error);

}

BOOST_AUTO_TEST_SUITE_END ()

} /* namespace IceBRG */
