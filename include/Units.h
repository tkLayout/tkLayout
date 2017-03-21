/*
 * Units.h
 *
 *  Created on: 6. 10. 2015
 *      Author: Z. Drasal (CERN)
 */

#ifndef INCLUDE_UNITS_H_
#define INCLUDE_UNITS_H_

namespace Units {

  static const float MeV = 1.0;
  static const float GeV = 1.E3;
  static const float TeV = 1.E6;

  static const float T   = 1.0;

  static const float mm  = 1.0;
  static const float cm  = 1.E1;
  static const float m   = 1.E3;
  static const float um  = 1.E-3;

  static const float g   = 1.0;
  static const float kg  = 1.E3;

  static const float mm2 = mm*mm;
  static const float cm2 = cm*cm;
  static const float m2  = m*m;

  static const float mm3 = mm*mm*mm;
  static const float cm3 = cm*cm*cm;
  static const float m3  = m*m*m;

  static const float s   = 1.0;

  static const float mW   = 1.E-3;
  static const float kW   = 1.E3;

  static const float k   = 1.E3;
  static const float M   = 1.E6;
  static const float G   = 1.E9;

  static const double deg = M_PI/180.; // one degree in radians
} // Namespace

#endif /* INCLUDE_UNITS_H_ */
