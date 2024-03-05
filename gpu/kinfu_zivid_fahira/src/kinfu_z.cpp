#ifndef KINFU_CPP
#define KINFU_CPP
// valinila kinfu
/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2011, Willow Garage, Inc.
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Willow Garage, Inc. nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <pcl/common/time.h>
#include <pcl/gpu/kinfu_zivid/kinfu_z.h>
#include <pcl/gpu/kinfu_zivid/volume_related.h>
#include <boost/asio/detail/atomic_count.hpp>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/SVD>

#include "internal.h"

#include <algorithm>
#include <iostream>

#ifdef HAVE_OPENCV
#include <opencv2/gpu/gpu.hpp>
#include <opencv2/opencv.hpp>
#endif

using namespace pcl::device;
using namespace pcl::gpu;

using Eigen::AngleAxisf;
using Eigen::Array3f;
using Eigen::Vector3f;
using Eigen::Vector3i;

std::string sfile = write_pose_file_name;
// std::ifstream pose_file("pose_icp_noise0changed.txt");// for 900 volume size
// std::ifstream pose_file2("noise0_500vol.txt"); // for 500 volume size
std::ifstream pose_file(POSE_FILE_NAME);
// settings
// but obtained with volume of 450
// other icp drifts towrads the
// end

std::vector<int>
generate_numbers()
{
  std::vector<int> numbers;

  for (int i = -1; i >= -179; --i) {
    numbers.push_back(i);
  }

  numbers.push_back(180);
  for (int i = 179; i >= 0; --i) {
    numbers.push_back(i);
  }
  numbers.erase(numbers.begin() + 35);
  return numbers;
}
// auto deg_angles = generate_numbers();
// int count_angles = 0;

Eigen::Matrix4f
params_to_matrix(std::vector<double> params)
{
  Eigen::Affine3f transform(Eigen::Affine3d::Identity());
  Eigen::Matrix3f rot;
  rot = Eigen::AngleAxisf(params[5], Eigen::Vector3f::UnitZ()) *
        Eigen::AngleAxisf(params[4], Eigen::Vector3f::UnitY()) *
        Eigen::AngleAxisf(params[3], Eigen::Vector3f::UnitX());

  transform.translate(Eigen::Vector3f(params[0], params[1], params[2]));
  transform.rotate(rot);

  return transform.matrix();
}

Eigen::Matrix4f
read_transformation_from_file(int check_line)
{

  if (check_line == 1) {
    if (!pose_file.is_open()) {
      std::cerr << "Error: Unable to open the file." << std::endl;
    }
    else {
      std::cout << "file is fine" << std::endl;
    }
  }
  pose_file.clear();
  pose_file.seekg(0, std::ios::beg);
  std::string line;
  for (int line_no = 0; line_no < check_line + 2; line_no++) {
    getline(pose_file, line);
    // std::cout << line_no << " " << check_line << " " << line << std::endl;
    if (line_no == (check_line))

    {
      // std::cout << line_no << " " << check_line << " " << line << std::endl;
      std::vector<double> parameters; // vector to store line
      std::istringstream line_ss(line);
      std::string double_the_entry;

      while (line_ss >> double_the_entry) {
        // std::cout << double_the_entry << std::endl;
        parameters.push_back(stod(double_the_entry));
      }
      // parameters[3] = deg_angles[count_angles] * PI / 180.f;
      // count_angles++;
      // std::cout << "done" << std::endl;
      std::cout << parameters.back() << std::endl;
      return params_to_matrix(parameters);
    }
  }
}

void
write_transform_to_text(
    Eigen::Matrix<float, 3, 3, Eigen::RowMajor> transformation_matrix,
    Eigen::Vector3f ts,
    int frame_count)
{
  double x = ts[0];
  double y = ts[1];
  double z = ts[2];

  double roll = std::atan2(transformation_matrix(2, 1), transformation_matrix(2, 2));
  double pitch = std::asin(-transformation_matrix(2, 0));
  double yaw = std::atan2(transformation_matrix(1, 0), transformation_matrix(0, 0));
  std::ofstream ofile(sfile, std::ios::out | std::ios::app);
  ofile << x << "  " << y << "  " << z << "  " << roll << "  " << pitch << "  " << yaw
        << " " << roll * 180.f / PI << "  " << pitch * 180.f / PI << "  "
        << yaw * 180.f / PI << " " << frame_count << std::endl;
  // ofile << std::endl;
  ofile.close();
}

namespace pcl {
namespace gpu {
Eigen::Vector3f
rodrigues2(const Eigen::Matrix3f& matrix);
}
} // namespace pcl

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
pcl::gpu::KinfuTracker::KinfuTracker(int rows, int cols, int noise_components)
: rows_(rows)
, cols_(cols)
, global_time_(0)
, noise_components_(noise_components)
, max_icp_distance_(0)
, integration_metric_threshold_(0.f)
, disable_icp_(false)
{
  const Vector3f volume_size = Vector3f::Constant(VOLUME_SIZE);
  const Vector3i volume_resolution(VOLUME_X, VOLUME_Y, VOLUME_Z);

  tsdf_volume_ = TsdfVolume::Ptr(new TsdfVolume(volume_resolution));
  tsdf_volume_->setSize(volume_size);

  setDepthIntrinsics(KINFU_DEFAULT_DEPTH_FOCAL_X,
                     KINFU_DEFAULT_DEPTH_FOCAL_Y); // default values, can be overwritten

  init_Rcam_ = Eigen::Matrix3f::Identity(); // * AngleAxisf(-30.f/180*3.1415926,
  if (!shrink_volume_) {                    // Vector3f::UnitX());
    init_tcam_ =
        volume_size * 0.5f -
        Vector3f(0, 0, volume_size(2) / 2 * 1.2f); // with volume size from 900 and on
    std::cout << "init_tcam_ with standard volume" << init_tcam_ << std::endl;
  }
  // if you change init_tcam_ here, please change "t" with the same values in
  // kinfu_app_zivid.cpp aroud
  // line 744 which becomes part of "pose"
  // in vector 3f first part adjusts the height,2nd adjusts the width and 3rd
  // adjusts the depth. idealy first and second should be roughly half of the
  // volume size
  if (shrink_volume_) {
    // float v = 0.15f;
    // init_tcam_ = Vector3f(0.4, 0.45, v) -
    //            Vector3f(0, 0, v * 3.f); // with volume_size 800,700
    // float vv = 0.2f;
    // init_tcam_ =
    //   Vector3f(0.3, 0.3, vv) - Vector3f(0, 0, vv * 4.f); // with volume 500.f

    // float vvv = 0.2f;
    // init_tcam_ = Vector3f(0.29, 0.24, vvv) -
    //            Vector3f(0, 0, vvv * 4.f); // with volume 450.f

    // /*
    // float vvvv = 0.2f;
    // init_tcam_ = Vector3f(0.25, 0.23, -0.7f); // with 300 volume shiva

    // init_tcam_ = Vector3f(0.12, 0.12, 0.2f); // with 300 volume ganesh0
    init_tcam_ = Vector3f(shrunk_volume_x, shrunk_volume_y, shrunk_volume_z);
    std::cout << "init_tcam_ " << init_tcam_ << std::endl;
    //*/ //-
    // bigger height, more top will be included
    // Vector3f(0, 0, vvvv * 4.4f); // with volume 300.f

    // std::cout << "volume size " << volume_size << " " << volume_size(2)
    //<< std::endl;
  }

  const int iters[] = {10, 5, 4};
  std::copy(iters, iters + LEVELS, icp_iterations_);

  const float default_distThres = default_threshold_dist_; // 0.10f; //meters
  const float default_angleThres = sin(20.f * 3.14159254f / 180.f);
  const float default_tranc_dist = trunc_dist_m; // 0.03f; //meters

  setIcpCorespFilteringParams(default_distThres, default_angleThres);
  tsdf_volume_->setTsdfTruncDist(default_tranc_dist);

  allocateBufffers(rows, cols);

  rmats_.reserve(30000);
  tvecs_.reserve(30000);

  reset();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::gpu::KinfuTracker::setDepthIntrinsics(float fx, float fy, float cx, float cy)
{
  fx_ = fx;
  fy_ = fy;
  cx_ = (cx == -1) ? cols_ / 2 - 0.5f : cx;
  cy_ = (cy == -1) ? rows_ / 2 - 0.5f : cy;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::gpu::KinfuTracker::getDepthIntrinsics(float& fx,
                                           float& fy,
                                           float& cx,
                                           float& cy) const
{
  fx = fx_;
  fy = fy_;
  cx = cx_;
  cy = cy_;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::gpu::KinfuTracker::setInitalCameraPose(const Eigen::Affine3f& pose)
{
  init_Rcam_ = pose.rotation();
  init_tcam_ = pose.translation();
  reset();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::gpu::KinfuTracker::setDepthTruncationForICP(float max_icp_distance)
{
  max_icp_distance_ = max_icp_distance;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::gpu::KinfuTracker::setCameraMovementThreshold(float threshold)
{
  integration_metric_threshold_ = threshold;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::gpu::KinfuTracker::setIcpCorespFilteringParams(float distThreshold,
                                                    float sineOfAngle)
{
  distThres_ = distThreshold; // mm
  angleThres_ = sineOfAngle;
}
void
pcl::gpu::KinfuTracker::setNoiseComponents(int noise_components)
{
  noise_components_ = noise_components;
}

void
pcl::gpu::KinfuTracker::setFrameCounter(int frame_count)
{
  frame_count_ = frame_count;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int
pcl::gpu::KinfuTracker::cols()
{
  return (cols_);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int
pcl::gpu::KinfuTracker::rows()
{
  return (rows_);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::gpu::KinfuTracker::reset(int num)
{
  if (num == 0) {
    if (global_time_)
      std::cout << "global time Reset" << std::endl;

    global_time_ = 0;
    rmats_.clear();
    tvecs_.clear();

    rmats_.push_back(init_Rcam_);
    tvecs_.push_back(init_tcam_);

    tsdf_volume_->reset();

    if (color_volume_) // color integration mode is enabled
      color_volume_->reset();
    std::cout << std::endl;
  }
  else {
    frame_count_++;
    std::cout << std::endl;
    std::cout << "global time Reset at " << num << std::endl;
  }
  resets_count++;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::gpu::KinfuTracker::allocateBufffers(int rows, int cols)
{
  depths_curr_.resize(LEVELS);
  vmaps_g_curr_.resize(LEVELS);
  nmaps_g_curr_.resize(LEVELS);

  vmaps_g_prev_.resize(LEVELS);
  nmaps_g_prev_.resize(LEVELS);

  vmaps_curr_.resize(LEVELS);
  nmaps_curr_.resize(LEVELS);

  coresps_.resize(LEVELS);

  for (int i = 0; i < LEVELS; ++i) {
    int pyr_rows = rows >> i;
    int pyr_cols = cols >> i;

    depths_curr_[i].create(pyr_rows, pyr_cols);

    vmaps_g_curr_[i].create(pyr_rows * 3, pyr_cols);
    nmaps_g_curr_[i].create(pyr_rows * 3, pyr_cols);

    vmaps_g_prev_[i].create(pyr_rows * 3, pyr_cols);
    nmaps_g_prev_[i].create(pyr_rows * 3, pyr_cols);

    vmaps_curr_[i].create(pyr_rows * 3, pyr_cols);
    nmaps_curr_[i].create(pyr_rows * 3, pyr_cols);

    coresps_[i].create(pyr_rows, pyr_cols);
  }
  depthRawScaled_.create(rows, cols);
  // see estimate transform for the magic numbers
  gbuf_.create(27, 20 * 60);
  sumbuf_.create(27);
}

// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// bool pcl::gpu::KinfuTracker::operator()(const DepthMap &depth_raw,
//                                         Eigen::Affine3f *hint) {
//   device::Intr intr(fx_, fy_, cx_, cy_);

//   if (!disable_icp_) {
//     {
//       // ScopeTime time(">>> Bilateral, pyr-down-all, create-maps-all");
//       //  depth_raw.copyTo(depths_curr_[0]);
//       device::bilateralFilter(depth_raw, depths_curr_[0]);

//       if (max_icp_distance_ > 0)
//         device::truncateDepth(depths_curr_[0], max_icp_distance_);

//       for (int i = 1; i < LEVELS; ++i)
//         device::pyrDown(depths_curr_[i - 1], depths_curr_[i]);

//       for (int i = 0; i < LEVELS; ++i) {
//         device::createVMap(intr(i), depths_curr_[i], vmaps_curr_[i]);
//         // device::createNMap(vmaps_curr_[i], nmaps_curr_[i]);
//         computeNormalsEigen(vmaps_curr_[i], nmaps_curr_[i]);
//       }
//       pcl::device::sync();
//     }

//     // can't perform more on first frame
//     if (global_time_ == 0) {
//       Matrix3frm init_Rcam = rmats_[0]; //  [Ri|ti] - pos of camera, i.e.
//       Vector3f init_tcam = tvecs_[0];   //  transform from camera to global coo
//                                         //  space for (i-1)th camera pose

//       Mat33 &device_Rcam = device_cast<Mat33>(init_Rcam);
//       float3 &device_tcam = device_cast<float3>(init_tcam);

//       Matrix3frm init_Rcam_inv = init_Rcam.inverse();
//       Mat33 &device_Rcam_inv = device_cast<Mat33>(init_Rcam_inv);
//       float3 device_volume_size =
//           device_cast<const float3>(tsdf_volume_->getSize());

//       // integrateTsdfVolume(depth_raw, intr, device_volume_size,
//       // device_Rcam_inv, device_tcam, tranc_dist, volume_);
//       device::integrateTsdfVolume(
//           depth_raw, nmaps_curr_[0], intr, device_volume_size, device_Rcam_inv,
//           device_tcam, tsdf_volume_->getTsdfTruncDist(), tsdf_volume_->data(),
//           depthRawScaled_, noise_components_);

//       for (int i = 0; i < LEVELS; ++i)
//         device::tranformMaps(vmaps_curr_[i], nmaps_curr_[i], device_Rcam,
//                              device_tcam, vmaps_g_prev_[i], nmaps_g_prev_[i]);

//       ++global_time_;
//       ++frame_count_;
//       // std::cout << POSE_FILE_NAME << std::endl;
//       std::cout << "reading file " << POSE_FILE_NAME << std::endl;
//       return (false);
//     }

//     ///////////////////////////////////////////////////////////////////////////////////////////
//     // Iterative Closest Point
//     Matrix3frm Rprev =
//         rmats_[global_time_ - 1];              //  [Ri|ti] - pos of camera, i.e.
//     Vector3f tprev = tvecs_[global_time_ - 1]; //  tranfrom from camera to
//                                                //  global coo space for (i-1)th
//                                                //  camera pose
//     Matrix3frm Rprev_inv = Rprev.inverse();    // Rprev.t();

//     // Mat33&  device_Rprev     = device_cast<Mat33> (Rprev);
//     Mat33 &device_Rprev_inv = device_cast<Mat33>(Rprev_inv);
//     float3 &device_tprev = device_cast<float3>(tprev);
//     Matrix3frm Rcurr;
//     Vector3f tcurr;
//     if (hint) {
//       Rcurr = hint->rotation().matrix();
//       tcurr = hint->translation().matrix();
//     } else {
//       Rcurr = Rprev; // transform to global coo for ith camera pose
//       tcurr = tprev;
//     }
//     {
//       // ScopeTime time("icp-all");
//       for (int level_index = LEVELS - 1; level_index >= 0; --level_index) {
//         int iter_num = icp_iterations_[level_index];

//         MapArr &vmap_curr = vmaps_curr_[level_index];
//         MapArr &nmap_curr = nmaps_curr_[level_index];

//         // MapArr& vmap_g_curr = vmaps_g_curr_[level_index];
//         // MapArr& nmap_g_curr = nmaps_g_curr_[level_index];

//         MapArr &vmap_g_prev = vmaps_g_prev_[level_index];
//         MapArr &nmap_g_prev = nmaps_g_prev_[level_index];

//         // CorespMap& coresp = coresps_[level_index];

//         for (int iter = 0; iter < iter_num; ++iter) {
//           Mat33 &device_Rcurr = device_cast<Mat33>(Rcurr);
//           float3 &device_tcurr = device_cast<float3>(tcurr);

//           Eigen::Matrix<double, 6, 6, Eigen::RowMajor> A;
//           Eigen::Matrix<double, 6, 1> b;
// #if 0
//             device::tranformMaps(vmap_curr, nmap_curr, device_Rcurr, device_tcurr,
//             vmap_g_curr, nmap_g_curr); findCoresp(vmap_g_curr, nmap_g_curr,
//             device_Rprev_inv, device_tprev, intr(level_index), vmap_g_prev,
//             nmap_g_prev, distThres_, angleThres_, coresp);
//             device::estimateTransform(vmap_g_prev, nmap_g_prev, vmap_g_curr, coresp,
//             gbuf_, sumbuf_, A.data(), b.data());

//             //cv::gpu::GpuMat ma(coresp.rows(), coresp.cols(), CV_32S, coresp.ptr(),
//             coresp.step());
//             //cv::Mat cpu;
//             //ma.download(cpu);
//             //cv::imshow(names[level_index] + string(" --- coresp white == -1"), cpu
//             == -1);
// #else
//           estimateCombined(device_Rcurr, device_tcurr, vmap_curr, nmap_curr,
//                            device_Rprev_inv, device_tprev, intr(level_index),
//                            vmap_g_prev, nmap_g_prev, distThres_, angleThres_,
//                            gbuf_, sumbuf_, A.data(), b.data());
// #endif
//           // checking nullspace
//           double det = A.determinant();

//           if (std::abs(det) < 1e-15 || std::isnan(det)) {
//             if (std::isnan(det))
//               std::cout << "qnan" << std::endl;

//             reset(global_time_);
//             return (false);
//           }
//           // float maxc = A.maxCoeff();

//           Eigen::Matrix<float, 6, 1> result = A.llt().solve(b).cast<float>();
//           // Eigen::Matrix<float, 6, 1> result = A.jacobiSvd(ComputeThinU |
//           // ComputeThinV).solve(b);

//           float alpha = result(0);
//           float beta = result(1);
//           float gamma = result(2);

//           Eigen::Matrix3f Rinc =
//               (Eigen::Matrix3f)AngleAxisf(gamma, Vector3f::UnitZ()) *
//               AngleAxisf(beta, Vector3f::UnitY()) *
//               AngleAxisf(alpha, Vector3f::UnitX());
//           Vector3f tinc = result.tail<3>();

//           // compose
//           tcurr = Rinc * tcurr + tinc;
//           Rcurr = Rinc * Rcurr;
//         }
//       }
//     }
//     if (read_saved_pose_) {
//       if (frame_count_ % 50 == 0) {
//         std::cout
//             << "currently set to run for 360 frames. after this it will do "
//                "seg fault"
//             << std::endl;
//       }
//       if ((frame_count_ > 0) && (frame_count_ < 361)) {
//         // int increment = 0;
//         // if (global_time_ > 35) { // at 36, it will look for pose at 37th
//         // pose increment = 1;
//         //}

//         // auto transform_mat = read_transformation_from_file(
//         //    global_time_ - 1); //  + increment); no need to add increment
//         //  in the file. at
//         //  36th frame we have 37
//         // degrees already

//         auto transform_mat = read_transformation_from_file(frame_count_ - 1);
//         Rcurr = transform_mat.block<3, 3>(0, 0);
//         tcurr = transform_mat.block<3, 1>(0, 3);
//         // std::cout << "transl " << tcurr.transpose() << " global_time_ "
//         //        << global_time_ << std::endl;
//         if (global_time_ == 358)
//           std::cout << "resets_count for global time " << resets_count
//                     << std::endl;
//       }
//     }

//     // std::cout << "glob time " << global_time_ << std::endl;

//     // save transform
//     rmats_.push_back(Rcurr);
//     tvecs_.push_back(tcurr);
//   } else /* if (disable_icp_) */
//   {
//     if (global_time_ == 0) {
//       ++global_time_;
//       ++frame_count_;
//     }

//     /*Matrix3frm Rcurr = rmats_[global_time_ - 1];
//     Vector3f tcurr = tvecs_[global_time_ - 1];
//      rmats_.push_back(Rcurr);
//     tvecs_.push_back(tcurr);
//     */
//     if (read_saved_pose_) {
//       if (frame_count_ % 50 == 0) {
//         std::cout
//             << "currently set to run for 360 frames. after this it will do "
//                "seg fault"
//             << std::endl;
//       }
//       if ((frame_count_ > 0) && (frame_count_ < 361)) {

//         auto transform_mat = read_transformation_from_file(frame_count_ - 1);
//         Matrix3frm Rcurr = transform_mat.block<3, 3>(0, 0);
//         Vector3f tcurr = transform_mat.block<3, 1>(0, 3);
//         // std::cout << "transl " << tcurr.transpose() << " global_time_ "
//         //        << global_time_ << std::endl;
//         if (global_time_ == 358)
//           std::cout << "resets_count for global time " << resets_count
//                     << std::endl;
//         rmats_.push_back(Rcurr);
//         tvecs_.push_back(tcurr);
//       }
//     }

//     // std::cout << "same? " << tcurr.transpose() << std::endl;
//   }

//   Matrix3frm Rprev = rmats_[global_time_ - 1];
//   Vector3f tprev = tvecs_[global_time_ - 1];

//   Matrix3frm Rcurr = rmats_.back();
//   Vector3f tcurr = tvecs_.back();

//   /*if (write_pose_to_file_)
//     write_transform_to_text(Rcurr, tcurr,frame_count);*/

//   ///////////////////////////////////////////////////////////////////////////////////////////
//   // Integration check - We do not integrate volume if camera does not move.
//   float rnorm = rodrigues2(Rcurr.inverse() * Rprev).norm();
//   float tnorm = (tcurr - tprev).norm();
//   const float alpha = 1.f;
//   bool integrate = (rnorm + alpha * tnorm) / 2 >= integration_metric_threshold_;

//   if (disable_icp_)
//     integrate = true;

//   ///////////////////////////////////////////////////////////////////////////////////////////
//   // Volume integration
//   float3 device_volume_size =
//       device_cast<const float3>(tsdf_volume_->getSize());

//   Matrix3frm Rcurr_inv = Rcurr.inverse();
//   Mat33 &device_Rcurr_inv = device_cast<Mat33>(Rcurr_inv);
//   float3 &device_tcurr = device_cast<float3>(tcurr);
//   if (integrate) {
//     // ScopeTime time("tsdf");
//     // integrateTsdfVolume(depth_raw, intr, device_volume_size,
//     // device_Rcurr_inv, device_tcurr, tranc_dist, volume_);
//     //      integrateWeightedTsdfVolume (depth_raw, nmaps_curr_[0], intr,
//     //      device_volume_size, device_Rcurr_inv, device_tcurr,
//     //      tsdf_volume_->getTsdfTruncDist(), tsdf_volume_->data(),
//     //      depthRawScaled_, noise_components_);
//     if (reconst_every_n_frame) {
//       if (global_time_ % n_th_frame == 0) {
//         integrateTsdfVolume(
//             depth_raw, nmaps_curr_[0], intr, device_volume_size,
//             device_Rcurr_inv, device_tcurr, tsdf_volume_->getTsdfTruncDist(),
//             tsdf_volume_->data(), depthRawScaled_, noise_components_);
//         // Ray casting
//         Mat33 &device_Rcurr = device_cast<Mat33>(Rcurr);
//         {
//           // ScopeTime time("ray-cast-all");
//           raycast(intr, device_Rcurr, device_tcurr,
//                   tsdf_volume_->getTsdfTruncDist(), device_volume_size,
//                   tsdf_volume_->data(), vmaps_g_prev_[0], nmaps_g_prev_[0]);
//           for (int i = 1; i < LEVELS; ++i) {
//             resizeVMap(vmaps_g_prev_[i - 1], vmaps_g_prev_[i]);
//             resizeNMap(nmaps_g_prev_[i - 1], nmaps_g_prev_[i]);
//           }
//           pcl::device::sync();
//         }
//         if (write_pose_to_file_)
//           write_transform_to_text(Rcurr, tcurr, frame_count_);
//       }
//     } else {
//       integrateTsdfVolume(
//           depth_raw, nmaps_curr_[0], intr, device_volume_size, device_Rcurr_inv,
//           device_tcurr, tsdf_volume_->getTsdfTruncDist(), tsdf_volume_->data(),
//           depthRawScaled_, noise_components_);
//       // Ray casting
//       Mat33 &device_Rcurr = device_cast<Mat33>(Rcurr);
//       {
//         // ScopeTime time("ray-cast-all");
//         raycast(intr, device_Rcurr, device_tcurr,
//                 tsdf_volume_->getTsdfTruncDist(), device_volume_size,
//                 tsdf_volume_->data(), vmaps_g_prev_[0], nmaps_g_prev_[0]);
//         for (int i = 1; i < LEVELS; ++i) {
//           resizeVMap(vmaps_g_prev_[i - 1], vmaps_g_prev_[i]);
//           resizeNMap(nmaps_g_prev_[i - 1], nmaps_g_prev_[i]);
//         }
//         pcl::device::sync();
//       }
//       if (write_pose_to_file_)
//         write_transform_to_text(Rcurr, tcurr, frame_count_);
//     }
//   }

//   ///////////////////////////////////////////////////////////////////////////////////////////

//   ++global_time_;
//   ++frame_count_;
//   return (true);
// }

bool
pcl::gpu::KinfuTracker::operator()(const DepthMap& depth_raw, const Eigen::Affine3f* hint)
{
  device::Intr intr(fx_, fy_, cx_, cy_);

  if (hint != nullptr) {


    // still do thesee to keep the same as the original kinfu
    {
      // ScopeTime time(">>> Bilateral, pyr-down-all, create-maps-all");
      // depth_raw.copyTo(depths_curr[0]);
      device::bilateralFilter(depth_raw, depths_curr_[0]);

      if (max_icp_distance_ > 0)
        device::truncateDepth(depths_curr_[0], max_icp_distance_);

      for (int i = 1; i < LEVELS; ++i)
        device::pyrDown(depths_curr_[i - 1], depths_curr_[i]);

      for (int i = 0; i < LEVELS; ++i) {
        device::createVMap(intr(i), depths_curr_[i], vmaps_curr_[i]);
        // device::createNMap(vmaps_curr_[i], nmaps_curr_[i]);
        computeNormalsEigen(vmaps_curr_[i], nmaps_curr_[i]);
      }
      pcl::device::sync();
    }
    disable_icp_ = true;
    bool integrate = false;
    Matrix3frm Rcurr;
    Vector3f tcurr;
    if (global_time_ == 0) {
      // Rcurr = init_Rcam_;
      // tcurr = init_tcam_;
      // Eigen::Affine3f init_pose = Eigen::Affine3f::Identity();
      // init_pose.rotate(Rcurr);
      // init_pose.translate(tcurr);
      // delta_pose_ = init_pose * (*hint).inverse();
      integrate = true;
      auto init_pose = (*hint);
      Rcurr = init_pose.rotation();
      tcurr = init_pose.translation();
      rmats_[0] = Rcurr;
      tvecs_[0] = tcurr;

    }
    else {
      // pose_prev.translate(Tprev);
      // Eigen::Affine3f pose_curr = delta_pose_ * (*hint);
      auto pose_curr = (*hint);
      Rcurr = pose_curr.rotation();
      tcurr = pose_curr.translation();

      rmats_.push_back(Rcurr);
      tvecs_.push_back(tcurr);
      // pose_prev.rotate(Rprev);
      Matrix3frm Rprev = rmats_[global_time_ - 1]; //  [Ri|ti] - pos of camera, i.e.
      Vector3f Tprev = tvecs_[global_time_ - 1];   //  (i-1)th camera pose

      float rnorm = rodrigues2(Rcurr.inverse() * Rprev).norm();
      float tnorm = (tcurr - Tprev).norm();
      const float alpha = 1.f;
      integrate = (rnorm + alpha * tnorm) / 2 >= integration_metric_threshold_;
    }
    // Volume integration
    float3 device_volume_size = device_cast<const float3>(tsdf_volume_->getSize());

    Matrix3frm Rcurr_inv = Rcurr.inverse();
    Mat33& device_Rcurr_inv = device_cast<Mat33>(Rcurr_inv);
    float3& device_tcurr = device_cast<float3>(tcurr);

    if (integrate) {
      integrateTsdfVolume(depth_raw,
                                  nmaps_curr_[0],
                                  intr,
                                  device_volume_size,
                                  device_Rcurr_inv,
                                  device_tcurr,
                                  tsdf_volume_->getTsdfTruncDist(),
                                  tsdf_volume_->data(),
                                  depthRawScaled_,
                                  noise_components_);
    }

    ///////////////////////////////////////////////////////////////////////////////////
    // Ray casting
    Mat33& device_Rcurr = device_cast<Mat33>(Rcurr);
    {
      // ScopeTime time("ray-cast-all");
      raycast(intr,
              device_Rcurr,
              device_tcurr,
              tsdf_volume_->getTsdfTruncDist(),
              device_volume_size,
              tsdf_volume_->data(),
              vmaps_g_prev_[0],
              nmaps_g_prev_[0]);
      for (int i = 1; i < LEVELS; ++i) {
        resizeVMap(vmaps_g_prev_[i - 1], vmaps_g_prev_[i]);
        resizeNMap(nmaps_g_prev_[i - 1], nmaps_g_prev_[i]);
      }
      pcl::device::sync();
    }

    ++global_time_;
    return (true);
  }

  if (!disable_icp_) {
    {
      // ScopeTime time(">>> Bilateral, pyr-down-all, create-maps-all");
      // depth_raw.copyTo(depths_curr[0]);
      device::bilateralFilter(depth_raw, depths_curr_[0]);

      if (max_icp_distance_ > 0)
        device::truncateDepth(depths_curr_[0], max_icp_distance_);

      for (int i = 1; i < LEVELS; ++i)
        device::pyrDown(depths_curr_[i - 1], depths_curr_[i]);

      for (int i = 0; i < LEVELS; ++i) {
        device::createVMap(intr(i), depths_curr_[i], vmaps_curr_[i]);
        // device::createNMap(vmaps_curr_[i], nmaps_curr_[i]);
        computeNormalsEigen(vmaps_curr_[i], nmaps_curr_[i]);
      }
      pcl::device::sync();
    }

    // can't perform more on first frame
    if (global_time_ == 0) {
      Matrix3frm init_Rcam = rmats_[0]; //  [Ri|ti] - pos of camera, i.e.
      Vector3f init_tcam = tvecs_[0]; //  transform from camera to global coo space for
                                      //  (i-1)th camera pose

      Mat33& device_Rcam = device_cast<Mat33>(init_Rcam);
      float3& device_tcam = device_cast<float3>(init_tcam);

      Matrix3frm init_Rcam_inv = init_Rcam.inverse();
      Mat33& device_Rcam_inv = device_cast<Mat33>(init_Rcam_inv);
      float3 device_volume_size = device_cast<const float3>(tsdf_volume_->getSize());

      // integrateTsdfVolume(depth_raw, intr, device_volume_size, device_Rcam_inv,
      // device_tcam, tranc_dist, volume_);
      // device::integrateTsdfVolume(depth_raw,
      //                             intr,
      //                             device_volume_size,
      //                             device_Rcam_inv,
      //                             device_tcam,
      //                             tsdf_volume_->getTsdfTruncDist(),
      //                             tsdf_volume_->data(),
      //                             depthRawScaled_);
      device::integrateTsdfVolume(depth_raw,
                                          nmaps_curr_[0],
                                          intr,
                                          device_volume_size,
                                          device_Rcam_inv,
                                          device_tcam,
                                          tsdf_volume_->getTsdfTruncDist(),
                                          tsdf_volume_->data(),
                                          depthRawScaled_,
                                          noise_components_);

      for (int i = 0; i < LEVELS; ++i)
        device::tranformMaps(vmaps_curr_[i],
                             nmaps_curr_[i],
                             device_Rcam,
                             device_tcam,
                             vmaps_g_prev_[i],
                             nmaps_g_prev_[i]);

      ++global_time_;
      return (false);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    // Iterative Closest Point
    Matrix3frm Rprev = rmats_[global_time_ - 1]; //  [Ri|ti] - pos of camera, i.e.
    //  tranfrom from camera to global coo space for (i-1)th camera pose
    Vector3f tprev = tvecs_[global_time_ - 1];
    Matrix3frm Rprev_inv = Rprev.inverse(); // Rprev.t();

    // Mat33&  device_Rprev     = device_cast<Mat33> (Rprev);
    Mat33& device_Rprev_inv = device_cast<Mat33>(Rprev_inv);
    float3& device_tprev = device_cast<float3>(tprev);
    Matrix3frm Rcurr =
        Rprev; // tranfrom from camera to global coo space for ith camera pose
    Vector3f tcurr = tprev;
    {
      // ScopeTime time("icp-all");
      for (int level_index = LEVELS - 1; level_index >= 0; --level_index) {
        int iter_num = icp_iterations_[level_index];

        MapArr& vmap_curr = vmaps_curr_[level_index];
        MapArr& nmap_curr = nmaps_curr_[level_index];

        // MapArr& vmap_g_curr = vmaps_g_curr_[level_index];
        // MapArr& nmap_g_curr = nmaps_g_curr_[level_index];

        MapArr& vmap_g_prev = vmaps_g_prev_[level_index];
        MapArr& nmap_g_prev = nmaps_g_prev_[level_index];

        // CorespMap& coresp = coresps_[level_index];

        for (int iter = 0; iter < iter_num; ++iter) {
          Mat33& device_Rcurr = device_cast<Mat33>(Rcurr);
          float3& device_tcurr = device_cast<float3>(tcurr);

          Eigen::Matrix<double, 6, 6, Eigen::RowMajor> A;
          Eigen::Matrix<double, 6, 1> b;
#if 0
            device::tranformMaps(vmap_curr, nmap_curr, device_Rcurr, device_tcurr, vmap_g_curr, nmap_g_curr);
            findCoresp(vmap_g_curr, nmap_g_curr, device_Rprev_inv, device_tprev, intr(level_index), vmap_g_prev, nmap_g_prev, distThres_, angleThres_, coresp);
            device::estimateTransform(vmap_g_prev, nmap_g_prev, vmap_g_curr, coresp, gbuf_, sumbuf_, A.data(), b.data());

            //cv::gpu::GpuMat ma(coresp.rows(), coresp.cols(), CV_32S, coresp.ptr(), coresp.step());
            //cv::Mat cpu;
            //ma.download(cpu);
            //cv::imshow(names[level_index] + string(" --- coresp white == -1"), cpu == -1);
#else
          estimateCombined(device_Rcurr,
                           device_tcurr,
                           vmap_curr,
                           nmap_curr,
                           device_Rprev_inv,
                           device_tprev,
                           intr(level_index),
                           vmap_g_prev,
                           nmap_g_prev,
                           distThres_,
                           angleThres_,
                           gbuf_,
                           sumbuf_,
                           A.data(),
                           b.data());
#endif
          // checking nullspace
          double det = A.determinant();

          if (std::abs(det) < 1e-15 || std::isnan(det)) {
            if (std::isnan(det))
              std::cout << "qnan" << std::endl;

            reset();
            return (false);
          }
          // float maxc = A.maxCoeff();

          Eigen::Matrix<float, 6, 1> result = A.llt().solve(b).cast<float>();
          // Eigen::Matrix<float, 6, 1> result = A.jacobiSvd(ComputeThinU |
          // ComputeThinV).solve(b);

          float alpha = result(0);
          float beta = result(1);
          float gamma = result(2);

          Eigen::Matrix3f Rinc = (Eigen::Matrix3f)AngleAxisf(gamma, Vector3f::UnitZ()) *
                                 AngleAxisf(beta, Vector3f::UnitY()) *
                                 AngleAxisf(alpha, Vector3f::UnitX());
          Vector3f tinc = result.tail<3>();

          // compose
          tcurr = Rinc * tcurr + tinc;
          Rcurr = Rinc * Rcurr;
        }
      }
    }
    // save transform
    rmats_.push_back(Rcurr);
    tvecs_.push_back(tcurr);
  }
  else /* if (disable_icp_) */
  {
    if (global_time_ == 0)
      ++global_time_;

    Matrix3frm Rcurr = rmats_[global_time_ - 1];
    Vector3f tcurr = tvecs_[global_time_ - 1];

    rmats_.push_back(Rcurr);
    tvecs_.push_back(tcurr);
  }

  Matrix3frm Rprev = rmats_[global_time_ - 1];
  Vector3f tprev = tvecs_[global_time_ - 1];

  Matrix3frm Rcurr = rmats_.back();
  Vector3f tcurr = tvecs_.back();

  ///////////////////////////////////////////////////////////////////////////////////////////
  // Integration check - We do not integrate volume if camera does not move.
  float rnorm = rodrigues2(Rcurr.inverse() * Rprev).norm();
  float tnorm = (tcurr - tprev).norm();
  const float alpha = 1.f;
  bool integrate = (rnorm + alpha * tnorm) / 2 >= integration_metric_threshold_;

  if (disable_icp_)
    integrate = true;

  ///////////////////////////////////////////////////////////////////////////////////////////
  // Volume integration
  float3 device_volume_size = device_cast<const float3>(tsdf_volume_->getSize());

  Matrix3frm Rcurr_inv = Rcurr.inverse();
  Mat33& device_Rcurr_inv = device_cast<Mat33>(Rcurr_inv);
  float3& device_tcurr = device_cast<float3>(tcurr);
  if (integrate) {
    // ScopeTime time("tsdf");
    // integrateTsdfVolume(depth_raw, intr, device_volume_size, device_Rcurr_inv,
    // device_tcurr, tranc_dist, volume_);
    // integrateTsdfVolume(depth_raw,
    //                     intr,
    //                     device_volume_size,
    //                     device_Rcurr_inv,
    //                     device_tcurr,
    //                     tsdf_volume_->getTsdfTruncDist(),
    //                     tsdf_volume_->data(),
    //                     depthRawScaled_);
    integrateTsdfVolume(depth_raw,
                                nmaps_curr_[0],
                                intr,
                                device_volume_size,
                                device_Rcurr_inv,
                                device_tcurr,
                                tsdf_volume_->getTsdfTruncDist(),
                                tsdf_volume_->data(),
                                depthRawScaled_,
                                noise_components_);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////
  // Ray casting
  Mat33& device_Rcurr = device_cast<Mat33>(Rcurr);
  {
    // ScopeTime time("ray-cast-all");
    raycast(intr,
            device_Rcurr,
            device_tcurr,
            tsdf_volume_->getTsdfTruncDist(),
            device_volume_size,
            tsdf_volume_->data(),
            vmaps_g_prev_[0],
            nmaps_g_prev_[0]);
    for (int i = 1; i < LEVELS; ++i) {
      resizeVMap(vmaps_g_prev_[i - 1], vmaps_g_prev_[i]);
      resizeNMap(nmaps_g_prev_[i - 1], nmaps_g_prev_[i]);
    }
    pcl::device::sync();
  }

  ++global_time_;
  return (true);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Affine3f
pcl::gpu::KinfuTracker::getCameraPose(int time) const
{
  if (time > (int)rmats_.size() || time < 0)
    time = rmats_.size() - 1;

  Eigen::Affine3f aff;
  aff.linear() = rmats_[time];
  aff.translation() = tvecs_[time];
  return (aff);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

size_t
pcl::gpu::KinfuTracker::getNumberOfPoses() const
{
  return rmats_.size();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const TsdfVolume&
pcl::gpu::KinfuTracker::volume() const
{
  return *tsdf_volume_;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TsdfVolume&
pcl::gpu::KinfuTracker::volume()
{
  return *tsdf_volume_;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const ColorVolume&
pcl::gpu::KinfuTracker::colorVolume() const
{
  return *color_volume_;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ColorVolume&
pcl::gpu::KinfuTracker::colorVolume()
{
  return *color_volume_;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::gpu::KinfuTracker::getImage(View& view) const
{
  // Eigen::Vector3f light_source_pose = tsdf_volume_->getSize() * (-3.f);
  Eigen::Vector3f light_source_pose = tvecs_[tvecs_.size() - 1];

  device::LightSource light;
  light.number = 1;
  light.pos[0] = device_cast<const float3>(light_source_pose);

  view.create(rows_, cols_);
  generateImage(vmaps_g_prev_[0], nmaps_g_prev_[0], light, view);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::gpu::KinfuTracker::getLastFrameCloud(DeviceArray2D<PointType>& cloud) const
{
  cloud.create(rows_, cols_);
  DeviceArray2D<float4>& c = (DeviceArray2D<float4>&)cloud;
  device::convert(vmaps_g_prev_[0], c);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::gpu::KinfuTracker::getLastFrameNormals(DeviceArray2D<NormalType>& normals) const
{
  normals.create(rows_, cols_);
  DeviceArray2D<float8>& n = (DeviceArray2D<float8>&)normals;
  device::convert(nmaps_g_prev_[0], n);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::gpu::KinfuTracker::disableIcp()
{
  disable_icp_ = true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
pcl::gpu::KinfuTracker::initColorIntegration(int max_weight)
{
  color_volume_ =
      pcl::gpu::ColorVolume::Ptr(new ColorVolume(*tsdf_volume_, max_weight));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
pcl::gpu::KinfuTracker::operator()(const DepthMap& depth, const View& colors, const Eigen::Affine3f* hint)
{
  bool res = (*this)(depth, hint);

  if (res && color_volume_) {
    const float3 device_volume_size =
        device_cast<const float3>(tsdf_volume_->getSize());
    device::Intr intr(fx_, fy_, cx_, cy_);

    Matrix3frm R_inv = rmats_.back().inverse();
    Vector3f t = tvecs_.back();

    Mat33& device_Rcurr_inv = device_cast<Mat33>(R_inv);
    float3& device_tcurr = device_cast<float3>(t);

    device::updateColorVolume(intr,
                              tsdf_volume_->getTsdfTruncDist(),
                              device_Rcurr_inv,
                              device_tcurr,
                              vmaps_g_prev_[0],
                              colors,
                              device_volume_size,
                              color_volume_->data(),
                              color_volume_->getMaxWeight());
  }

  return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace pcl {
namespace gpu {
PCL_EXPORTS void
paint3DView(const KinfuTracker::View& rgb24,
            KinfuTracker::View& view,
            float colors_weight = 0.5f)
{
  device::paint3DView(rgb24, view, colors_weight);
}

PCL_EXPORTS void
mergePointNormal(const DeviceArray<PointXYZ>& cloud,
                 const DeviceArray<Normal>& normals,
                 DeviceArray<PointNormal>& output)
{
  const std::size_t size = std::min(cloud.size(), normals.size());
  output.create(size);

  const DeviceArray<float4>& c = (const DeviceArray<float4>&)cloud;
  const DeviceArray<float8>& n = (const DeviceArray<float8>&)normals;
  const DeviceArray<float12>& o = (const DeviceArray<float12>&)output;
  device::mergePointNormal(c, n, o);
}

Eigen::Vector3f
rodrigues2(const Eigen::Matrix3f& matrix)
{
  Eigen::JacobiSVD<Eigen::Matrix3f> svd(matrix,
                                        Eigen::ComputeFullV | Eigen::ComputeFullU);
  Eigen::Matrix3f R = svd.matrixU() * svd.matrixV().transpose();

  double rx = R(2, 1) - R(1, 2);
  double ry = R(0, 2) - R(2, 0);
  double rz = R(1, 0) - R(0, 1);

  double s = sqrt((rx * rx + ry * ry + rz * rz) * 0.25);
  double c = (R.trace() - 1) * 0.5;
  c = c > 1. ? 1. : c < -1. ? -1. : c;

  double theta = std::acos(c);

  if (s < 1e-5) {
    if (c > 0)
      rx = ry = rz = 0;
    else {
      double t;
      t = (R(0, 0) + 1) * 0.5;
      rx = sqrt(std::max(t, 0.0));
      t = (R(1, 1) + 1) * 0.5;
      ry = sqrt(std::max(t, 0.0)) * (R(0, 1) < 0 ? -1.0 : 1.0);
      t = (R(2, 2) + 1) * 0.5;
      rz = sqrt(std::max(t, 0.0)) * (R(0, 2) < 0 ? -1.0 : 1.0);

      if (std::abs(rx) < std::abs(ry) && std::abs(rx) < std::abs(rz) &&
          (R(1, 2) > 0) != (ry * rz > 0))
        rz = -rz;
      theta /= sqrt(rx * rx + ry * ry + rz * rz);
      rx *= theta;
      ry *= theta;
      rz *= theta;
    }
  }
  else {
    double vth = 1 / (2 * s);
    vth *= theta;
    rx *= vth;
    ry *= vth;
    rz *= vth;
  }
  return Eigen::Vector3d(rx, ry, rz).cast<float>();
}
} // namespace gpu
} // namespace pcl

#endif /* KINFU_CPP */
