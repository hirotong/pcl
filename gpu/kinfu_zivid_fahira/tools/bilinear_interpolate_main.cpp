#ifndef BILINEAR_INTERPOLATE_MAIN_CPP
#define BILINEAR_INTERPOLATE_MAIN_CPP
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/program_options.hpp>
#include "noise_data.h"
#include "bilinear_interpolate_functions.cpp"

int main(int argc, char *argv[]) {

  namespace po = boost::program_options;

  po::options_description desc("Options");
  desc.add_options()("help", "Show the help message")(
      "dist_value", po::value<float>(),
      " float/int value for a distance between 370 and 1170")(
      "angle_value", po::value<float>(),
      " float/int value for a angle between 0 and 70");

  // Configure positional arguments
  po::positional_options_description pos_args;
  pos_args.add("dist_value", 1);
  pos_args.add("angle_value", 1);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv)
                .options(desc)
                .positional(pos_args)
                .run(),
            vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("dist_value") || !vm.count("angle_value")) {
    std::cout << desc << std::endl;
    return 1;
  }

  // Target coordinates for interpolation
  float target_x = vm["dist_value"].as<float>();  // 375;
  float target_y = vm["angle_value"].as<float>(); // 15;

  noise_data p1, p2, p3, p4;

  // check if the given value matches exactly with one in the data
  for (size_t i = 0; i < data.size(); ++i) {
    if ((data[i].dist == target_x) && (data[i].angle == target_y)) {

      std::cout << "this is exact data point, we do not need interpolation"
                << std::endl;

      std::cout << "the exact value at (" << target_x << ", " << target_y
                << ") is " << data[i].noise << std::endl;
      return 0;
    }
  }
  //=================================================
  // check if only distance is the same only. we reach here if above block
  // returns false only
  bool got_points = 0;
  for (size_t i = 0; i < data.size(); ++i) {
    if ((data[i].dist == target_x) && (data[i].angle < target_y) &&
        ((std::abs(data[i].angle - target_y)) <= 10)) {

      p1 = data[i];
      std::cout << "the selected point1 to interpolate angle is " << p1.dist
                << " " << p1.angle << " " << p1.noise << std::endl;
    }
  }

  for (size_t i = 0; i < data.size(); ++i) {
    if ((data[i].dist == target_x) && (data[i].angle > target_y) &&
        ((std::abs(data[i].angle - target_y)) <= 10)) {

      p2 = data[i];
      std::cout << "the selected point2 to interpolate angle is " << p2.dist
                << " " << p2.angle << " " << p2.noise << std::endl;
      // return 0;
      got_points = 1;
    }
  }
  if (got_points) {
    double interpolated_value1 =
        l_interpolate_angle(p1, p2, target_x, target_y);
    std::cout << "Interpolated value at (" << target_x << ", " << target_y
              << ") is " << interpolated_value1 << std::endl;
    return 0;
  }

  //==============================================================
  // check if we have only angle matches but need to interpolate distances
  got_points = 0;
  for (size_t i = 0; i < data.size(); ++i) {
    if ((data[i].angle == target_y) && (data[i].dist < target_x) &&
        (data[i].dist >= target_x - 25)) {
      p1 = data[i];
      std::cout << "the selected point1 is to interpolate dist is " << p1.dist
                << " " << p1.angle << " " << p1.noise << std::endl;
    }
  }

  for (size_t i = 0; i < data.size(); ++i) {
    if ((data[i].angle == target_y) && (data[i].dist >= target_x) &&
        (data[i].dist <= target_x + 25)) {
      p2 = data[i];
      std::cout << "the selected point2 is to interpolate dist is " << p2.dist
                << " " << p2.angle << " " << p2.noise << std::endl;
      got_points = 1;
    }
  }
  if (got_points) {
    double interpolated_value2 = l_interpolate_dist(p1, p2, target_x, target_y);
    std::cout << "Interpolated value at (" << target_x << ", " << target_y
              << ") is " << interpolated_value2 << std::endl;
    return 0;
  }

  //===================================================
  // Find the four nearest points surrounding the target point

  for (size_t i = 0; i < data.size(); ++i) {
    if ((data[i].dist < target_x) && (data[i].dist >= target_x - 25) &&
        (data[i].angle < target_y) &&
        ((std::abs(data[i].angle - target_y)) <= 10)) {
      p1 = data[i];

      break;
    }
  }
  for (size_t i = 0; i < data.size(); ++i) {
    if ((data[i].dist >= target_x) && (data[i].dist <= target_x + 25) &&
        (data[i].angle < target_y) &&
        ((std::abs(data[i].angle - target_y)) <= 10)) {
      p2 = data[i];

      break;
    }
  }

  for (size_t i = 0; i < data.size(); ++i) {
    if ((data[i].dist < target_x) && (data[i].dist >= target_x - 25) &&
        (data[i].angle > target_y) &&
        ((std::abs(data[i].angle - target_y)) <= 10)) {
      p3 = data[i];

      break;
    }
  }

  for (size_t i = 0; i < data.size(); ++i) {
    if ((data[i].dist >= target_x) && (data[i].dist <= target_x + 25) &&
        (data[i].angle > target_y) &&
        ((std::abs(data[i].angle - target_y)) <= 10)) {
      p4 = data[i];
      break;
    }
  }

  // Perform bilinear interpolation

  double interpolated_value = b_interpolate(p1, p2, p3, p4, target_x, target_y);

  std::cout << "Interpolated value at (" << target_x << ", " << target_y
            << ") is " << interpolated_value << std::endl;

  return 0;
}

#endif /* BILINEAR_INTERPOLATE_MAIN_CPP */
