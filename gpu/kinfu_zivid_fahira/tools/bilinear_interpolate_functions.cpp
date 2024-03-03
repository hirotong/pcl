#ifndef BILINEAR_INTERPOLATE_FUNCTIONS_CPP
#define BILINEAR_INTERPOLATE_FUNCTIONS_CPP

#include <iostream>
#include <vector>
#include <cmath>
//#include "noise_data.h"

// wikipedia bilinear interpolate

double b_interpolate(const noise_data &p1, const noise_data &p2,
                     const noise_data &p3, const noise_data &p4, double x,
                     double y) {
  double x1 = p1.dist, q11 = p1.noise;
  double x2 = p2.dist, q21 = p2.noise;
  double y1 = p1.angle, q12 = p3.noise;
  double y2 = p3.angle, q22 = p4.noise;

  std::cout << "picked valus" << x1 << " " << y1 << " " << q11 << " " << x2
            << " " << y1 << " " << q21 << " " << x1 << " " << y2 << " " << q12
            << " " << x2 << " " << y2 << " " << q22 << std::endl;

  std::cout << "picked points" << p1.dist << " " << p1.angle << " p2 "
            << p2.dist << " " << p2.angle << " p3 " << p3.dist << " "
            << p3.angle << " p4 " << p4.dist << " " << p4.angle << " "
            << std::endl;

  double x20 = x2 - x;
  double x21 = x2 - x1;
  double x01 = x - x1;

  double f_x_y1 = (x20 / x21) * q11 + (x01 / x21) * q21;
  double f_x_y2 = (x20 / x21) * q12 + (x01 / x21) * q22;

  double y20 = y2 - y;
  double y21 = y2 - y1;
  double y01 = y - y1;

  double f_x_y = (y20 / y21) * f_x_y1 + (y01 / y21) * f_x_y2;

  return f_x_y;
}
// linear interpolate just angles given distance already exists in the data
double l_interpolate_angle(const noise_data &p1, const noise_data &p2, double x,
                           double y) {
  double y1 = p1.angle, q11 = p1.noise;
  double y2 = p2.angle, q21 = p2.noise;

  double y20 = y2 - y;
  double y21 = y2 - y1;
  double y01 = y - y1;

  double f_x_y = (y20 / y21) * q11 + (y01 / y21) * q21;
  return f_x_y;
}
// linear interpolate just distance given angle already exists in the data
double l_interpolate_dist(const noise_data &p1, const noise_data &p2, double x,
                          double y) {
  double x1 = p1.dist, q11 = p1.noise;
  double x2 = p2.dist, q21 = p2.noise;

  double x20 = x2 - x;
  double x21 = x2 - x1;
  double x01 = x - x1;

  double f_x_y = (x20 / x21) * q11 + (x01 / x21) * q21;
  return f_x_y;
}

#endif /* BILINEAR_INTERPOLATE_FUNCTIONS_CPP */
