/*
 * Units.h
 *
 *  Created on: 6. 10. 2015
 *      Author: Z. Drasal (CERN)
 */

#ifndef INCLUDE_UNITS_H_
#define INCLUDE_UNITS_H_

#include <math.h>

namespace Units {

  static const float MeV = 1.0;
  static const float GeV = 1E3;
  static const float TeV = 1E6;

  static const float T   = 1.0;

  static const float mm  = 1.0;
  static const float cm  = 1E1;
  static const float m   = 1E3;
  static const float um  = 1E-3;

  static const float g   = 1.0;
  static const float kg  = 1E3;

  static const float mm2 = mm*mm;
  static const float cm2 = cm*cm;
  static const float m2  = m*m;

  static const float mm3 = mm*mm*mm;
  static const float cm3 = cm*cm*cm;
  static const float m3  = m*m*m;

  static const float k   = 1E3;
  static const float M   = 1E6;
  static const float G   = 1E9;

  static const float s   = 1.0;

  static const float   b   = 1.0;
  static const float   B   = 8.0;
  static const float   kb  = 1.0 * 1024;
  static const double  Mb  = 1.0 * 1024 * 1024;
  static const double  Gb  = 1.0 * 1024 * 1024 * 1024;
  static const double  Tb  = 1.0 * 1024 * 1024 * 1024 * 1024;
  static const float   kB  = 8.0 * 1024;
  static const float   MB  = 8.0 * 1024 * 1024;
  static const double  GB  = 8.0 * 1024 * 1024 * 1024;
  static const double  TB  = 8.0 * 1024 * 1024 * 1024 * 1024;

  static const double percent = 1E-2;
  static const double degrees = M_PI/180.;
  static const double radians = 1;
} // Namespace

#endif /* INCLUDE_UNITS_H_ */
