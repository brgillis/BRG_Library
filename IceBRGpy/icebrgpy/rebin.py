""" @file /disk2/brg/git/brg_library/IceBRGpy/icebrgpy/rebin.py

    Created 4 Mar 2016

    @TODO: File docstring

    ---------------------------------------------------------------------

    Copyright (C) 2016 brg

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
"""

from copy import deepcopy
import numpy as np

try:
    import IceBRGpy
except ImportError as _e:
    from Release import IceBRGpy


def rebin(a, x_offset=0, y_offset=0, subsampling_factor=5, conserve=False):
    """ Rebins an array with a given offset and subsampling factor. Note that
        unless 'conserve' is set to True, the input array will be overwritten.
    """
    
    # If we want to conserve, do so by operating on a copy of the array
    if(conserve):
        a = deepcopy(a)
    
    # Use the proper function for the data type
    if a.dtype == 'float32':
        new_shape = IceBRGpy.rebin_float(a,x_offset,y_offset,subsampling_factor)
    elif a.dtype == 'float64':
        new_shape = IceBRGpy.rebin_double(a,x_offset,y_offset,subsampling_factor)
    elif a.dtype == 'int32':
        new_shape = IceBRGpy.rebin_int(a,x_offset,y_offset,subsampling_factor)
    elif a.dtype == 'int64':
        new_shape = IceBRGpy.rebin_long(a,x_offset,y_offset,subsampling_factor)
    elif a.dtype == 'uint32':
        new_shape = IceBRGpy.rebin_uint(a,x_offset,y_offset,subsampling_factor)
    elif a.dtype == 'uint64':
        new_shape = IceBRGpy.rebin_ulong(a,x_offset,y_offset,subsampling_factor)
    else:
        raise Exception("Unsupported data type for rebinning: " + str(a.dtype))
    
    # Resort the new array into the proper shape
    new_size = np.product(new_shape)
    
    rebinned_array = np.reshape(np.ravel(a)[0:new_size],new_shape)
    
    return rebinned_array