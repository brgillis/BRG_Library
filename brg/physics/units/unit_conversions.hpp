/**********************************************************************\
  @file unit_conversions.hpp

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

#ifndef _BRG_UNIT_CONVERSIONS_HPP_INCLUDED_
#define _BRG_UNIT_CONVERSIONS_HPP_INCLUDED_

namespace brgastro{ namespace unitconv {
// All unit conversions are exact unless noted

// Distance
// Default unit: meter (m)
constexpr double mtom = 1;
constexpr double mtomm = 1e3;
constexpr double mmtom = 1 / mtomm;
constexpr double mtocm = 1e2;
constexpr double cmtom = 1 / mtocm;
constexpr double mtoum = 1e6;
constexpr double umtom = 1 / mtoum;
constexpr double mtonm = 1e9;
constexpr double nmtom = 1 / mtonm;
constexpr double mtoangstrom = 1e10;
constexpr double angstromtom = 1 / mtoangstrom;
constexpr double mtokm = 1e-3;
constexpr double kmtom = 1 / mtokm;
constexpr double ltyrtom = 9460730472580800;
constexpr double mtoltyr = 1 / ltyrtom;
constexpr double AUtom = 149597870700;
constexpr double mtoAU = 1 / AUtom;
constexpr double pctom = AUtom * 648000 / pi;
constexpr double mtopc = 1 / pctom;
constexpr double kpctom = 1000 * pctom;
constexpr double mtokpc = 1 / kpctom;
constexpr double Mpctom = 1000000 * pctom;
constexpr double mtoMpc = 1 / Mpctom;
constexpr double mitom = 1609.344;
constexpr double mtomi = 1 / mitom;
constexpr double Mmitom = 1e6 * mitom;
constexpr double mtoMmi = 1 / Mmitom;
constexpr double fttom = 0.3048;
constexpr double mtoft = 1 / fttom;
constexpr double intom = .0254;
constexpr double mtoin = 1 / intom;
constexpr double ydtom = 0.9144;
constexpr double mtoyd = 1 / ydtom;

// Time
// Default unit: second (s)
constexpr double stos = 1;
constexpr double stocs = 1e2;
constexpr double cstos = 1 / stocs;
constexpr double stoms = 1e3;
constexpr double mstos = 1 / stoms;
constexpr double stous = 1e6;
constexpr double ustos = 1 / stous;
constexpr double stons = 1e9;
constexpr double nstos = 1 / stons;
constexpr double mintos = 60;
constexpr double stomin = 1 / mintos;
constexpr double hrtos = mintos * 60;
constexpr double stohr = 1 / hrtos;
constexpr double daytos = hrtos * 24;
constexpr double stoday = 1 / daytos;
constexpr double weektos = daytos * 7;
constexpr double stoweek = 1 / weektos;
constexpr double yrtos = daytos * 365.24219; // Approximate tropical year
constexpr double stoyr = 1 / yrtos; // Approximate
constexpr double monthtos = yrtos / 12; // Mean month length for tropical year
constexpr double stomonth = 1 / monthtos; // Approximate
constexpr double kyrtos = yrtos * 1e3; // Approximate
constexpr double stokyr = 1 / kyrtos; // Approximate
constexpr double Myrtos = yrtos * 1e6; // Approximate
constexpr double stoMyr = 1 / Myrtos; // Approximate
constexpr double Gyrtos = yrtos * 1e9; // Approximate
constexpr double stoGyr = 1 / Gyrtos; // Approximate

// Velocity
// Default units: meters per second (mps)
constexpr double mpstomps = 1;
constexpr double mpstokmps = 1e-3;
constexpr double kmpstomps = 1 / mpstokmps;
constexpr double ctomps = 299792458;
constexpr double mpstoc = 1 / ctomps;
constexpr double mpstomiphr = mtomi / stohr;
constexpr double miphr = 1 / mpstomiphr;

// Mass
// Default unit: kilogram (kg)
constexpr double kgtokg = 1;
constexpr double kgtogm = 1e3;
constexpr double gmtokg = 1 / kgtogm;
constexpr double Mearthtokg = 5.9736e24; // Approximate
constexpr double kgtoMearth = 1 / Mearthtokg; // Approximate
constexpr double Msuntokg = 1.9891e30; // Approximate
constexpr double kgtoMsun = 1 / Msuntokg; // Approximate
constexpr double kgtottMsun = kgtoMsun * 1e-10; // Approximate
constexpr double ttMsuntokg = 1 / kgtottMsun; // Approximate

// Temperature
// Default unit: Kelvin (K)
constexpr double KtoK = 1;
constexpr double KtodegF = 1.8;
constexpr double degCtoK = KtoK;
constexpr double KtodegC = 1 / degCtoK;
constexpr double degFtoK = 1 / KtodegF;
constexpr double degCtodegF = KtodegF;
constexpr double degFtodegC = degFtoK;
constexpr double KtodegR = KtodegF;
constexpr double degRtoK = degFtoK;
constexpr double degCtodegR = KtodegF;
constexpr double degRtodegC = degFtoK;

// Angle
// Default unit: radian (rad)
constexpr double radtorad = 1;
constexpr double degtorad = pi / 180;
constexpr double radtodeg = 1 / degtorad;
constexpr double degtoamin = 60;
constexpr double amintodeg = 1 / degtoamin;
constexpr double amintoasec = 60;
constexpr double asectoamin = 1 / amintoasec;
constexpr double asectodeg = asectoamin * amintodeg;
constexpr double degtoasec = 1 / asectodeg;
constexpr double amintorad = amintodeg * degtorad;
constexpr double radtoamin = 1 / amintorad;
constexpr double asectorad = asectodeg * degtorad;
constexpr double radtoasec = 1 / asectorad;

// Charge
// Default unit: Coulomb (C)
constexpr double CtoC = 1;
constexpr double Ctoesu = 6.241509324e18; // Approximate
constexpr double esutoC = 1 / Ctoesu;

} } // namespace brgastro::unitconv



#endif /* _BRG_UNIT_CONVERSIONS_HPP_INCLUDED_ */
