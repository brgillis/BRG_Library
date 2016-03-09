#!/usr/bin/env python

""" @file /disk2/brg/git/brg_library/IceBRGpy/plot_n_of_m.py

    Created 23 Feb 2016

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

import numpy as np
import sys
from matplotlib import pyplot

import IceBRG

def main(argv):
    """ @TODO main docstring
    """
    
    l10_masses = np.linspace(8, 16, 100)
    
    z1 = 0.
    z2 = 0.5
    z3 = 1.0
    z4 = 1.5
    z5 = 2.0
    
    n_z1 = np.zeros_like(l10_masses)
    n_z2 = np.zeros_like(l10_masses)
    n_z3 = np.zeros_like(l10_masses)
    n_z4 = np.zeros_like(l10_masses)
    n_z5 = np.zeros_like(l10_masses)
    
    for i in range(len(l10_masses)):
        n_z1[i] = np.log10(IceBRG.log10_mass_function(l10_masses[i],z1)*np.power(IceBRG.Mpctom,3))
        n_z2[i] = np.log10(IceBRG.log10_mass_function(l10_masses[i],z2)*np.power(IceBRG.Mpctom,3))
        n_z3[i] = np.log10(IceBRG.log10_mass_function(l10_masses[i],z3)*np.power(IceBRG.Mpctom,3))
        n_z4[i] = np.log10(IceBRG.log10_mass_function(l10_masses[i],z4)*np.power(IceBRG.Mpctom,3))
        n_z5[i] = np.log10(IceBRG.log10_mass_function(l10_masses[i],z5)*np.power(IceBRG.Mpctom,3))
        
    fig = pyplot.figure()
    
    ax = fig.add_subplot(1,1,1)
    
    ax.plot(l10_masses,n_z1,label="z="+str(z1))
    ax.plot(l10_masses,n_z2,label="z="+str(z2))
    ax.plot(l10_masses,n_z3,label="z="+str(z3))
    ax.plot(l10_masses,n_z4,label="z="+str(z4))
    ax.plot(l10_masses,n_z5,label="z="+str(z5))
    
    ax.set_ylim(-10,2)
    
    ax.legend(loc='lower left')
    
    pyplot.show()

if __name__ == "__main__":
    main(sys.argv)
