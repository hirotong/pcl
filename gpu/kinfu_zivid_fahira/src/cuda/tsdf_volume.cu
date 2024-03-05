#ifndef TSDF_VOLUME_CU
#define TSDF_VOLUME_CU
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

#include "device.hpp"
// #include "bilinear.cu" //performs bilateral interpolation for noise
#include <thrust/binary_search.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h> //to bring std vector on cpu to device vector

#include "noise_data.h"

bool apply_zivid_noise = apply_zivid_noise_; // apply zivid noise or kinect noise model
__device__ float sig_z_default_value = sig_z_default_value_;
__device__ float factor_ = noise_scaling_factor;
//===========================================

//==========================================
int data_size = data.size();
thrust::device_vector<noise_data> d_data_ = data;
thrust::device_vector<float> d_distances = distances;
thrust::device_vector<float> d_angles = angles;
int dis_size = distances.size(), ang_size = angles.size();
noise_data* d_data_ptr =
    thrust::raw_pointer_cast(d_data_.data()); // Get raw pointer to device data
float* d_dis_ptr = thrust::raw_pointer_cast(d_distances.data());
float* d_ang_ptr = thrust::raw_pointer_cast(d_angles.data());
//===================================
#ifndef BILINEAR_CU
#define BILINEAR_CU

#include <pcl/gpu/kinfu_zivid/kinfu_z.h>
#include <pcl/gpu/kinfu_zivid/volume_related.h>

// performs bilateral interpolations on text data

// wikipedia bilinear interpolate
__device__ double
b_interpolate2(const noise_data& p1,
               const noise_data& p2,
               const noise_data& p3,
               const noise_data& p4,
               double x,
               double y,
               float factor_)
{

  // printf("factor=%f\n", factor_);
  /*ex point for 380,15
  p1 (370,10), p2(395,10), p3(370,20),p4(395,20)*/
  double x1 = p1.dist / 1000.f, q11 = factor_ * p1.noise / 1000.f;
  double x2 = p2.dist / 1000.f, q21 = factor_ * p2.noise / 1000.f;
  double y1 = p1.angle * Pi / 180.f, q12 = factor_ * p3.noise / 1000.f;
  double y2 = p3.angle * Pi / 180.f, q22 = factor_ * p4.noise / 1000.f;

  double x20 = x2 - x;
  double x21 = x2 - x1;
  double x01 = x - x1;

  double f_x_y1 = (x20 / x21) * q11 + (x01 / x21) * q21;
  double f_x_y2 = (x20 / x21) * q12 + (x01 / x21) * q22;

  double y20 = y2 - y;
  double y21 = y2 - y1;
  double y01 = y - y1;

  double f_x_y = (y20 / y21) * f_x_y1 + (y01 / y21) * f_x_y2;
  // printf("bilinear sigmaz =%f,p1 angle=%f, p3 angle=%f,p1 dist=%f, p2 dist=%f,
  // Dp_scaled=%f, theta=%f\n", f_x_y, p1.angle, p3.angle, p1.dist, p2.dist, x, y *
  // 180.f / PI);

  return f_x_y;
}
// linear interpolate just angles given distance already exists in the data
__device__ double
l_interpolate_angle2(
    const noise_data& p1, const noise_data& p2, double x, double y, float factor_)
{
  double y1 = p1.angle * Pi / 180.f, q11 = factor_ * p1.noise / 1000;
  double y2 = p2.angle * Pi / 180.f, q21 = factor_ * p2.noise / 1000;

  double y20 = y2 - y;
  double y21 = y2 - y1;
  double y01 = y - y1;

  double f_x_y = (y20 / y21) * q11 + (y01 / y21) * q21;
  printf("sigmaz %f,p1%f, p2=%f, Dp_scaled=%f, theta=%f\n",
         f_x_y,
         p1.angle,
         p2.angle,
         x,
         y * 180.f / Pi);

  return f_x_y;
}
// linear interpolate just distance given angle already exists in the data
__device__ double
l_interpolate_dist2(
    const noise_data& p1, const noise_data& p2, double x, double y, float factor_)
{
  // double factor_ = 10.f;
  double x1 = p1.dist / 1000.f, q11 = factor_ * p1.noise / 1000.f;
  double x2 = p2.dist / 1000.f, q21 = factor_ * p2.noise / 1000.f;

  double x20 = x2 - x;
  double x21 = x2 - x1;
  double x01 = x - x1;

  double f_x_y = (x20 / x21) * q11 + (x01 / x21) * q21;

  // printf("sigmaz_l_inter %f,p1_ang=%f, p2_ang=%f,p1.dist=%f,p2.dist=%f Dp_scaled=%f,
  // theta=%f\n", f_x_y, p1.angle, p2.angle, p1.dist, p2.dist, x, y * 180.f / PI);

  return f_x_y;
}
//------------------------------points for angle interpolation
__device__ noise_data
get_p1_for_angle_interpolate2(float Dp_scaled,
                              const float theta,
                              const noise_data* d_data,
                              int data_size)
{
  noise_data p1 = {0, 0, 0};

  for (size_t index = 0; index < data_size; ++index) {
    const noise_data& element = d_data[index];
    if ((element.dist / 1000.f == Dp_scaled) && (element.angle < theta * 180.f / PI) &&
        ((std::abs(element.angle - theta * 180.f / PI)) <= 10)) {
      p1 = element;
      printf("p1 %f=%f\n", p1.angle);
      break;
    }
  }
  return p1;
}
__device__ noise_data
get_p2_for_angle_interpolate2(float Dp_scaled,
                              const float theta,
                              const noise_data* d_data,
                              int data_size)
{
  noise_data p2 = {0, 0, 0};
  for (size_t index = 0; index < data_size; ++index) {
    const noise_data& element = d_data[index];
    if ((element.dist / 1000.f == Dp_scaled) && (element.angle > theta * 180.f / PI) &&
        ((std::abs(element.angle - theta * 180.f / PI)) <= 10)) {
      p2 = element;
      printf("p2 %f=%f\n", p2.angle);
      // sigmaz = l_interpolate_angle(p1, p2, Dp_scaled, theta);
      // cond_satisfied = 1;
      break;
    }
  }
  return p2;
}
//=======================//points for dist interpolation
__device__ void
get_neighbours_for_dist_interpolation2()
{}

//------------------------------
__device__ float
compute_sigma2_(float Dp_scaled,
                const float theta,
                float sig_z_default_value,
                float sigmaz,
                const noise_data* d_data,
                int data_size,
                float factor_)
{
  bool check_next_condition = 0;
  bool check_next_condition1 = 0;
  bool check_next_condition2 = 0;

  // when we have exact data points, we do not neeed interpolation
  for (size_t index = 0; index < data_size; ++index) {
    const noise_data& element = d_data[index];
    if ((element.dist / 1000.f == Dp_scaled) &&
        (element.angle == theta * 180.f / PI)) // apparently this will never happen
    {
      sigmaz = factor_ * element.noise;
      printf("Element %d: dist = %f, angle = %f, noise = %f\n",
             index,
             element.dist,
             element.angle,
             element.noise);
      break;
    }
  }

  // if (sigmaz == -100.f)
  if (sigmaz == sig_z_default_value) {
    check_next_condition = 1;
  }
  if (sigmaz != sig_z_default_value) {
    return sigmaz;
  }

  //------------------------------------------angle interpolate------------
  if (check_next_condition) {

    noise_data p1, p2;
    bool cond_satisfied = 0;
    //  index = blockIdx.x * blockDim.x + threadIdx.x;
    //  if (index < data_size)
    for (size_t index = 0; index < data_size; ++index) {
      const noise_data& element = d_data[index];
      if ((element.dist / 1000.f == Dp_scaled) &&
          (element.angle < theta * 180.f / PI) &&
          ((std::abs(element.angle - theta * 180.f / PI)) <= 10)) {
        p1 = element;
        printf("p1 %f=%f\n", p1.angle);
        break;
      }
    }
    for (size_t index = 0; index < data_size; ++index) {
      const noise_data& element = d_data[index];
      if ((element.dist / 1000.f == Dp_scaled) &&
          (element.angle > theta * 180.f / PI) &&
          ((std::abs(element.angle - theta * 180.f / PI)) <= 10)) {
        p2 = element;
        printf("p2 %f=%f\n", p2.angle);
        // sigmaz = l_interpolate_angle(p1, p2, Dp_scaled, theta);
        cond_satisfied = 1;
        break;
      }
    }
    if (cond_satisfied) {
      sigmaz = l_interpolate_angle2(p1, p2, Dp_scaled, theta, factor_);
    }

    // if (sigmaz != -100.f)
    if (sigmaz != sig_z_default_value) {
      printf("angle_interpolatep1.angle=%f,p2.angle=%f,p1.dist=%f,p2.dist=%f,Dp_scaled="
             "%f,theta=%f,sigmaz=%f\n",
             p1.angle,
             p2.angle,
             p1.dist,
             p2.dist,
             Dp_scaled,
             theta,
             sigmaz);
      return sigmaz;
    }
    // if (sigmaz == -100.f)
    if (sigmaz == sig_z_default_value) {
      check_next_condition1 = 1;
    }
  }
  //==================================dist interpolate====================
  if ((check_next_condition1) | (theta == 0.0f)) {

    noise_data p1, p2;
    bool cond_satisfied1 = 0;
    // index = blockIdx.x * blockDim.x + threadIdx.x;
    //  if (index < data_size)
    // printf("distance interepolation\n");
    for (size_t index = 0; index < data_size; ++index) {
      const noise_data& element = d_data[index];
      if ((element.angle == theta * 180.f / PI) &&
          (element.dist / 1000.f < Dp_scaled) &&
          (element.dist / 1000.f >= Dp_scaled - 0.025)) {
        p1 = element;
        break;
      }
    }
    for (size_t index = 0; index < data_size; ++index) {
      const noise_data& element = d_data[index];
      if ((element.angle == theta * 180.f / PI) &&
          (element.dist / 1000.f >= Dp_scaled) &&
          (element.dist / 1000.f <= Dp_scaled + 0.025)) {
        p2 = element;
        cond_satisfied1 = 1;
        /* p5 = p1; // doing just for printing
         p6 = p2;
         p7 = p1;
         p8 = p2;
        */
        break;
      }
    }
    if (cond_satisfied1) {
      sigmaz = l_interpolate_dist2(p1, p2, Dp_scaled, theta, factor_);
    }

    // if (sigmaz == -100.f)
    if (sigmaz == sig_z_default_value) {
      check_next_condition2 = 1;
    }
    if (sigmaz != sig_z_default_value) {
      return sigmaz;
    }
  }
  //------------------bilinear------------------------------------=====================================
  if (check_next_condition2) {

    noise_data p1, p2, p3, p4;
    for (size_t index = 0; index < data_size; ++index) {
      const noise_data& element = d_data[index];

      if ((element.dist / 1000.f < Dp_scaled) &&
          (element.dist / 1000.f >= Dp_scaled - 0.025) &&
          (element.angle < theta * 180.f / PI) &&
          ((std::abs(element.angle - theta * 180.f / PI)) <= 10)) {
        p1 = element;
        break;
      }
    }
    for (size_t index = 0; index < data_size; ++index) {
      const noise_data& element = d_data[index];

      if ((element.dist / 1000.f >= Dp_scaled) &&
          (element.dist / 1000.f <= Dp_scaled + 0.025) &&
          (element.angle < theta * 180.f / PI) &&
          ((std::abs(element.angle - theta * 180.f / PI)) <= 10)) {
        p2 = element;
        break;
      }
    }
    for (size_t index = 0; index < data_size; ++index) {
      const noise_data& element = d_data[index];

      if ((element.dist / 1000.f < Dp_scaled) &&
          (element.dist / 1000.f >= Dp_scaled - 0.025) &&
          (element.angle > theta * 180.f / PI) &&
          ((std::abs(element.angle - theta * 180.f / PI)) <= 10)) {
        p3 = element;
        break;
      }
    }
    for (size_t index = 0; index < data_size; ++index) {
      const noise_data& element = d_data[index];

      if ((element.dist / 1000.f >= Dp_scaled) &&
          (element.dist / 1000.f <= Dp_scaled + 0.025) &&
          (element.angle > theta * 180.f / PI) &&
          ((std::abs(element.angle - theta * 180.f / PI)) <= 10)) {
        p4 = element;
        break;
      }
    }
    /* p5 = p1;
     p6 = p2;
     p7 = p3;
     p8 = p4;
     */

    sigmaz = b_interpolate2(p1, p2, p3, p4, Dp_scaled, theta, factor_);
    return sigmaz;
  }
}

__device__ int
find_insert_position(const float* sorted_arry, size_t array_size, float new_value)
{
  for (int i = 0; i <  array_size; ++i) {
    if (new_value < sorted_arry[i]) {
      return i;
    }
  }
  return array_size;
}

__device__ float
compute_sigma2(float Dp_scaled,
               const float theta,
               float sig_z_default_value,
               float sigmaz,
               const noise_data* d_data,
               int data_size,
               float factor_,
               const float* d_dis_ptr,
               const float* d_ang_ptr,
               int dis_size,
               int ang_size)
{

  auto idx_dis = find_insert_position(d_dis_ptr, dis_size, Dp_scaled * 1000);
  auto idx_ang = find_insert_position(d_ang_ptr, ang_size, theta * 180.f / PI);

  if (idx_dis == dis_size) {
    idx_dis = dis_size - 1;
  }
  if (idx_ang == ang_size) {
    idx_ang = ang_size - 1;
  }

  int d_dis = idx_dis == 0 ? 0 : (idx_dis == dis_size - 1) ? 0 : -1;
  int d_ang = idx_ang == 0 ? 0 : (idx_ang == ang_size) ? 0 : dis_size;

  int idx_tr = dis_size - idx_dis - 1 + idx_ang * dis_size;
  int idx_tl = idx_tr - d_ang;
  int idx_br = idx_tr - d_dis;
  int idx_bl = idx_tr - d_dis - d_ang;
  auto p1 = d_data[idx_bl];
  auto p2 = d_data[idx_tl];
  auto p3 = d_data[idx_br];
  auto p4 = d_data[idx_tr];

  sigmaz = b_interpolate2(p1, p2, p3, p4, Dp_scaled, theta, factor_);
  return sigmaz;
}

#endif /* BILINEAR_CU */

//===================================

using namespace pcl::device;

namespace pcl {
namespace device {
template <typename T>
__global__ void
initializeVolume(PtrStep<T> volume)
{
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;

  if (x < VOLUME_X && y < VOLUME_Y) {
    T* pos = volume.ptr(y) + x;
    int z_step = VOLUME_Y * volume.step / sizeof(*pos);

#pragma unroll
    for (int z = 0; z < VOLUME_Z; ++z, pos += z_step)
      pack_tsdf(0.f, 0, *pos);
  }
}
} // namespace device
} // namespace pcl

void
pcl::device::initVolume(PtrStep<short2> volume)
{
  dim3 block(16, 16);
  dim3 grid(1, 1, 1);
  grid.x = divUp(VOLUME_X, block.x);
  grid.y = divUp(VOLUME_Y, block.y);

  initializeVolume<<<grid, block>>>(volume);
  cudaSafeCall(cudaGetLastError());
  cudaSafeCall(cudaDeviceSynchronize());
}

namespace pcl {
namespace device {
struct Tsdf {
  static constexpr int CTA_SIZE_X = 32;
  static constexpr int CTA_SIZE_Y = 8;
  static constexpr int MAX_WEIGHT = 1 << 7;

  mutable PtrStep<short2> volume;
  float3 cell_size;

  Intr intr;

  Mat33 Rcurr_inv;
  float3 tcurr;

  PtrStepSz<ushort> depth_raw; // depth in mm

  float tranc_dist_mm;

  __device__ __forceinline__ float3
  getVoxelGCoo(int x, int y, int z) const
  {
    float3 coo = make_float3(x, y, z);
    coo += 0.5f; // shift to cell center;

    coo.x *= cell_size.x;
    coo.y *= cell_size.y;
    coo.z *= cell_size.z;

    return coo;
  }

  __device__ __forceinline__ void
  operator()() const
  {
    int x = threadIdx.x + blockIdx.x * CTA_SIZE_X;
    int y = threadIdx.y + blockIdx.y * CTA_SIZE_Y;

    if (x >= VOLUME_X || y >= VOLUME_Y)
      return;

    short2* pos = volume.ptr(y) + x;
    int elem_step = volume.step * VOLUME_Y / sizeof(*pos);

    for (int z = 0; z < VOLUME_Z; ++z, pos += elem_step) {
      float3 v_g = getVoxelGCoo(x, y, z); // 3 // p

      // transform to curr cam coo space
      float3 v = Rcurr_inv * (v_g - tcurr); // 4

      int2 coo; // project to current cam
      coo.x = __float2int_rn(v.x * intr.fx / v.z + intr.cx);
      coo.y = __float2int_rn(v.y * intr.fy / v.z + intr.cy);

      if (v.z > 0 && coo.x >= 0 && coo.y >= 0 && coo.x < depth_raw.cols &&
          coo.y < depth_raw.rows) // 6
      {
        int Dp = depth_raw.ptr(coo.y)[coo.x];

        if (Dp != 0) {
          float xl = (coo.x - intr.cx) / intr.fx;
          float yl = (coo.y - intr.cy) / intr.fy;
          float lambda_inv = rsqrtf(xl * xl + yl * yl + 1);

          float sdf = 1000 * norm(tcurr - v_g) * lambda_inv - Dp; // mm

          sdf *= (-1);

          if (sdf >= -tranc_dist_mm) {
            float tsdf = fmin(1.f, sdf / tranc_dist_mm);

            float weight_prev;
            float tsdf_prev;

            // read and unpack
            unpack_tsdf(*pos, tsdf_prev, weight_prev);

            const int Wrk = 1;

            float tsdf_new =
                (tsdf_prev * weight_prev + Wrk * tsdf) / (weight_prev + Wrk);
            float weight_new = min(weight_prev + Wrk, float(MAX_WEIGHT));

            pack_tsdf(tsdf_new, weight_new, *pos);
          }
        }
      }
    }
  }
};

__global__ void
integrateTsdfKernel(const Tsdf tsdf)
{
  tsdf();
}

__global__ void
tsdf2(PtrStep<short2> volume,
      const float tranc_dist_mm,
      const Mat33 Rcurr_inv,
      float3 tcurr,
      const Intr intr,
      const PtrStepSz<ushort> depth_raw,
      const float3 cell_size)
{
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;

  if (x >= VOLUME_X || y >= VOLUME_Y)
    return;

  short2* pos = volume.ptr(y) + x;
  int elem_step = volume.step * VOLUME_Y / sizeof(short2);

  float v_g_x = (x + 0.5f) * cell_size.x - tcurr.x;
  float v_g_y = (y + 0.5f) * cell_size.y - tcurr.y;
  float v_g_z = (0 + 0.5f) * cell_size.z - tcurr.z;

  float v_x = Rcurr_inv.data[0].x * v_g_x + Rcurr_inv.data[0].y * v_g_y +
              Rcurr_inv.data[0].z * v_g_z;
  float v_y = Rcurr_inv.data[1].x * v_g_x + Rcurr_inv.data[1].y * v_g_y +
              Rcurr_inv.data[1].z * v_g_z;
  float v_z = Rcurr_inv.data[2].x * v_g_x + Rcurr_inv.data[2].y * v_g_y +
              Rcurr_inv.data[2].z * v_g_z;

  // #pragma unroll
  for (int z = 0; z < VOLUME_Z; ++z) {
    float3 vr;
    vr.x = v_g_x;
    vr.y = v_g_y;
    vr.z = (v_g_z + z * cell_size.z);

    float3 v;
    v.x = v_x + Rcurr_inv.data[0].z * z * cell_size.z;
    v.y = v_y + Rcurr_inv.data[1].z * z * cell_size.z;
    v.z = v_z + Rcurr_inv.data[2].z * z * cell_size.z;

    int2 coo; // project to current cam
    coo.x = __float2int_rn(v.x * intr.fx / v.z + intr.cx);
    coo.y = __float2int_rn(v.y * intr.fy / v.z + intr.cy);

    if (v.z > 0 && coo.x >= 0 && coo.y >= 0 && coo.x < depth_raw.cols &&
        coo.y < depth_raw.rows) // 6
    {
      int Dp = depth_raw.ptr(coo.y)[coo.x]; // mm

      if (Dp != 0) {
        float xl = (coo.x - intr.cx) / intr.fx;
        float yl = (coo.y - intr.cy) / intr.fy;
        float lambda_inv = rsqrtf(xl * xl + yl * yl + 1);

        float sdf = Dp - norm(vr) * lambda_inv * 1000; // mm

        if (sdf >= -tranc_dist_mm) {
          float tsdf = fmin(1.f, sdf / tranc_dist_mm);

          float weight_prev;
          float tsdf_prev;

          // read and unpack
          unpack_tsdf(*pos, tsdf_prev, weight_prev);

          const int Wrk = 1;

          float tsdf_new = (tsdf_prev * weight_prev + Wrk * tsdf) / (weight_prev + Wrk);
          float weight_new = min(weight_prev + Wrk, float(Tsdf::MAX_WEIGHT));

          pack_tsdf(tsdf_new, weight_new, *pos);
        }
      }
    }
    pos += elem_step;
  } /* for(int z = 0; z < VOLUME_Z; ++z) */
} /* __global__ */
} // namespace device
} // namespace pcl

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::device::integrateTsdfVolume(const PtrStepSz<ushort>& depth_raw,
                                 const Intr& intr,
                                 const float3& volume_size,
                                 const Mat33& Rcurr_inv,
                                 const float3& tcurr,
                                 float tranc_dist,
                                 PtrStep<short2> volume)
{
  Tsdf tsdf;

  tsdf.volume = volume;
  tsdf.cell_size.x = volume_size.x / VOLUME_X;
  tsdf.cell_size.y = volume_size.y / VOLUME_Y;
  tsdf.cell_size.z = volume_size.z / VOLUME_Z;

  tsdf.intr = intr;

  tsdf.Rcurr_inv = Rcurr_inv;
  tsdf.tcurr = tcurr;
  tsdf.depth_raw = depth_raw;

  tsdf.tranc_dist_mm = tranc_dist * 1000; // mm

  dim3 block(Tsdf::CTA_SIZE_X, Tsdf::CTA_SIZE_Y);
  dim3 grid(divUp(VOLUME_X, block.x), divUp(VOLUME_Y, block.y));

#if 0
   //tsdf2<<<grid, block>>>(volume, tranc_dist, Rcurr_inv, tcurr, intr, depth_raw, tsdf.cell_size);
   integrateTsdfKernel<<<grid, block>>>(tsdf);
#endif
  cudaSafeCall(cudaGetLastError());
  cudaSafeCall(cudaDeviceSynchronize());
}

namespace pcl {
namespace device {
__global__ void
scaleDepth(const PtrStepSz<ushort> depth, PtrStep<float> scaled, const Intr intr)
{
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;

  if (x >= depth.cols || y >= depth.rows)
    return;

  int Dp = depth.ptr(y)[x];

  float xl = (x - intr.cx) / intr.fx;
  float yl = (y - intr.cy) / intr.fy;
  float lambda = sqrtf(xl * xl + yl * yl + 1);

  scaled.ptr(y)[x] = Dp * lambda / 1000.f; // meters
}

__global__ void
tsdf23(const PtrStepSz<float> depthScaled,
       PtrStep<short2> volume,
       const float tranc_dist,
       const Mat33 Rcurr_inv,
       const float3 tcurr,
       const Intr intr,
       const float3 cell_size)
{
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;

  if (x >= VOLUME_X || y >= VOLUME_Y)
    return;

  float v_g_x = (x + 0.5f) * cell_size.x - tcurr.x;
  float v_g_y = (y + 0.5f) * cell_size.y - tcurr.y;
  float v_g_z = (0 + 0.5f) * cell_size.z - tcurr.z;

  float v_g_part_norm = v_g_x * v_g_x + v_g_y * v_g_y;

  float v_x = (Rcurr_inv.data[0].x * v_g_x + Rcurr_inv.data[0].y * v_g_y +
               Rcurr_inv.data[0].z * v_g_z) *
              intr.fx;
  float v_y = (Rcurr_inv.data[1].x * v_g_x + Rcurr_inv.data[1].y * v_g_y +
               Rcurr_inv.data[1].z * v_g_z) *
              intr.fy;
  float v_z = (Rcurr_inv.data[2].x * v_g_x + Rcurr_inv.data[2].y * v_g_y +
               Rcurr_inv.data[2].z * v_g_z);

  float z_scaled = 0;

  float Rcurr_inv_0_z_scaled = Rcurr_inv.data[0].z * cell_size.z * intr.fx;
  float Rcurr_inv_1_z_scaled = Rcurr_inv.data[1].z * cell_size.z * intr.fy;

  float tranc_dist_inv = 1.0f / tranc_dist;

  short2* pos = volume.ptr(y) + x;
  int elem_step = volume.step * VOLUME_Y / sizeof(short2);

  // #pragma unroll
  for (int z = 0; z < VOLUME_Z; ++z,
           v_g_z += cell_size.z,
           z_scaled += cell_size.z,
           v_x += Rcurr_inv_0_z_scaled,
           v_y += Rcurr_inv_1_z_scaled,
           pos += elem_step) {
    float inv_z = 1.0f / (v_z + Rcurr_inv.data[2].z * z_scaled);
    if (inv_z < 0)
      continue;

    // project to current cam
    int2 coo = {__float2int_rn(v_x * inv_z + intr.cx),
                __float2int_rn(v_y * inv_z + intr.cy)};

    if (coo.x >= 0 && coo.y >= 0 && coo.x < depthScaled.cols &&
        coo.y < depthScaled.rows) // 6
    {
      float Dp_scaled = depthScaled.ptr(coo.y)[coo.x]; // meters

      float sdf = Dp_scaled - sqrtf(v_g_z * v_g_z + v_g_part_norm);

      if (Dp_scaled != 0 && sdf >= -tranc_dist) // meters
      {
        float tsdf = fmin(1.0f, sdf * tranc_dist_inv);

        // read and unpack
        float tsdf_prev;
        float weight_prev;
        unpack_tsdf(*pos, tsdf_prev, weight_prev);

        const int Wrk = 1;

        float tsdf_new = (tsdf_prev * weight_prev + Wrk * tsdf) / (weight_prev + Wrk);
        float weight_new = min(weight_prev + Wrk, float(Tsdf::MAX_WEIGHT));

        pack_tsdf(tsdf_new, weight_new, *pos);
      }
    }
  } // for(int z = 0; z < VOLUME_Z; ++z)
} // __global__

//==========================================kinect noise
// fucntions=================================
//================================================kinect axial noise====================
__global__ void
    tsdf23_axial_kinect /*(const PtrStepSz<float> depthScaled, PtrStep<short2> volume,
                        const float tranc_dist, const Mat33 Rcurr_inv, const float3
                        tcurr, const Intr intr, const float3 cell_size)*/

    (const PtrStepSz<float> depthScaled,
     const PtrStep<float> nmap,
     int rows,
     int cols, /*added*/
     PtrStep<short2> volume,
     const Mat33 Rcurr_inv,
     const float3 tcurr,
     const Intr intr,
     const float3 cell_size)
{
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;

  if (x >= VOLUME_X || y >= VOLUME_Y)
    return;

  float v_g_x = (x + 0.5f) * cell_size.x - tcurr.x;
  float v_g_y = (y + 0.5f) * cell_size.y - tcurr.y;
  float v_g_z = (0 + 0.5f) * cell_size.z - tcurr.z;

  float v_g_part_norm = v_g_x * v_g_x + v_g_y * v_g_y;

  float v_x = (Rcurr_inv.data[0].x * v_g_x + Rcurr_inv.data[0].y * v_g_y +
               Rcurr_inv.data[0].z * v_g_z) *
              intr.fx;
  float v_y = (Rcurr_inv.data[1].x * v_g_x + Rcurr_inv.data[1].y * v_g_y +
               Rcurr_inv.data[1].z * v_g_z) *
              intr.fy;
  float v_z = (Rcurr_inv.data[2].x * v_g_x + Rcurr_inv.data[2].y * v_g_y +
               Rcurr_inv.data[2].z * v_g_z);

  float z_scaled = 0;

  float Rcurr_inv_0_z_scaled = Rcurr_inv.data[0].z * cell_size.z * intr.fx;
  float Rcurr_inv_1_z_scaled = Rcurr_inv.data[1].z * cell_size.z * intr.fy;

  short2* pos = volume.ptr(y) + x;
  int elem_step = volume.step * VOLUME_Y / sizeof(short2);

  // #pragma unroll
  for (int z = 0; z < VOLUME_Z; ++z,
           v_g_z += cell_size.z,
           z_scaled += cell_size.z,
           v_x += Rcurr_inv_0_z_scaled,
           v_y += Rcurr_inv_1_z_scaled,
           pos += elem_step) {
    float inv_z = 1.0f / (v_z + Rcurr_inv.data[2].z * z_scaled);
    if (inv_z < 0)
      continue;

    // project to current cam
    int2 coo = {__float2int_rn(v_x * inv_z + intr.cx),
                __float2int_rn(v_y * inv_z + intr.cy)};

    if (coo.x >= 0 && coo.y >= 0 && coo.x < depthScaled.cols &&
        coo.y < depthScaled.rows) // 6
    {
      float Dp_scaled = depthScaled.ptr(coo.y)[coo.x]; // meters

      // calculate sigma_L and sigma_z
      float3 nvect;
      nvect.x = nmap.ptr(coo.y)[coo.x];
      nvect.y = nmap.ptr(coo.y + rows)[coo.x];
      nvect.z = nmap.ptr(coo.y + 2 * rows)[coo.x];
      const float theta = acosf(fabsf(nvect.z));
      if (theta > PI / 2.5f)
        continue;
      const float sigmaz = 0.0012f + 0.0019f * (Dp_scaled - 0.4f) * (Dp_scaled - 0.4f) +
                           0.0001f / sqrtf(Dp_scaled) * theta * theta /
                               (PI / 2.f - theta) / (PI / 2.f - theta);
      float sdf = Dp_scaled - sqrtf(v_g_z * v_g_z + v_g_part_norm);

      if (Dp_scaled != 0 && sdf >= -6.f * sigmaz) // meters
      {
        float tsdf = sqrtf(1.0f - __expf(-2.f / PI * sdf * sdf / (sigmaz * sigmaz)));
        if (sdf < 0)
          tsdf = -tsdf;

        // read and unpack
        float tsdf_prev, weight_prev;
        unpack_tsdf(*pos, tsdf_prev, weight_prev);

        const float Wrk =
            0.0012f / sigmaz * (0.4f / Dp_scaled) *
            (0.4f / Dp_scaled); //* __expf(-du * du / (2.f * sigmaL * sigmaL) - dz * dz
                                /// (2.f * sigmaz * sigmaz)); // range from ~0 to 1.0

        float tsdf_new = (tsdf_prev * weight_prev + Wrk * tsdf) / (weight_prev + Wrk);
        float weight_new = min(weight_prev + Wrk, float(2 * Tsdf::MAX_WEIGHT));

        pack_tsdf(tsdf_new, weight_new, *pos);
      }
    }
  } // for(int z = 0; z < VOLUME_Z; ++z)
} // __global__
  //================================tsdf kinect both
  // noises=======================================================
__global__ void
    tsdf23_both_noise_kinect /*(const PtrStepSz<float> depthScaled, PtrStep<short2>
                              volume, const float tranc_dist, const Mat33 Rcurr_inv,
                              const float3 tcurr, const Intr intr, const float3
                              cell_size)*/
    (const PtrStepSz<float> depthScaled,
     const PtrStep<float> nmap,
     int rows,
     int cols, /*added*/
     PtrStep<short2> volume,
     const Mat33 Rcurr_inv,
     const float3 tcurr,
     const Intr intr,
     const float3 cell_size)
{
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;

  if (x >= VOLUME_X || y >= VOLUME_Y)
    return;

  float v_g_x = (x + 0.5f) * cell_size.x - tcurr.x;
  float v_g_y = (y + 0.5f) * cell_size.y - tcurr.y;
  float v_g_z = (0 + 0.5f) * cell_size.z - tcurr.z;

  float v_g_part_norm = v_g_x * v_g_x + v_g_y * v_g_y;

  float v_x = (Rcurr_inv.data[0].x * v_g_x + Rcurr_inv.data[0].y * v_g_y +
               Rcurr_inv.data[0].z * v_g_z) *
              intr.fx;
  float v_y = (Rcurr_inv.data[1].x * v_g_x + Rcurr_inv.data[1].y * v_g_y +
               Rcurr_inv.data[1].z * v_g_z) *
              intr.fy;
  float v_z = (Rcurr_inv.data[2].x * v_g_x + Rcurr_inv.data[2].y * v_g_y +
               Rcurr_inv.data[2].z * v_g_z);

  float z_scaled = 0;

  float Rcurr_inv_0_z_scaled = Rcurr_inv.data[0].z * cell_size.z * intr.fx;
  float Rcurr_inv_1_z_scaled = Rcurr_inv.data[1].z * cell_size.z * intr.fy;

  // float tranc_dist_inv = 1.0f / tranc_dist;

  short2* pos = volume.ptr(y) + x;
  int elem_step = volume.step * VOLUME_Y / sizeof(short2);

  // #pragma unroll
  for (int z = 0; z < VOLUME_Z; ++z,
           v_g_z += cell_size.z,
           z_scaled += cell_size.z,
           v_x += Rcurr_inv_0_z_scaled,
           v_y += Rcurr_inv_1_z_scaled,
           pos += elem_step) {
    float inv_z = 1.0f / (v_z + Rcurr_inv.data[2].z * z_scaled);
    if (inv_z < 0)
      continue;

    // project to current cam
    int2 coo = {__float2int_rn(v_x * inv_z + intr.cx),
                __float2int_rn(v_y * inv_z + intr.cy)};
    // and subpixel
    float2 dcoo = {(coo.x - (v_x * inv_z + intr.cx)),
                   (coo.y - (v_y * inv_z + intr.cy))};

    if (coo.x >= 1 && coo.y >= 1 && coo.x < depthScaled.cols - 1 &&
        coo.y < depthScaled.rows - 1) {
      float Dp_scaled_mid = depthScaled.ptr(coo.y)[coo.x]; // meters
      if (Dp_scaled_mid == 0)
        continue;

      // read and unpack
      float tsdf_prev, weight_prev;
      unpack_tsdf(*pos, tsdf_prev, weight_prev);
      // loop within 3x3 pixel area around the current position

#pragma unroll
      for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
          const int ii = coo.x + i;
          const int jj = coo.y + j;
          const float Dp_scaled = depthScaled.ptr(jj)[ii]; // meters
          const float dz = fabs(Dp_scaled - Dp_scaled_mid);
          const float du =
              sqrtf((i + dcoo.x) * (i + dcoo.x) + (j + dcoo.y) * (j + dcoo.y));

          // calculate sigma_L and sigma_z
          float3 nvect;
          nvect.x = nmap.ptr(jj)[ii];
          nvect.y = nmap.ptr(jj + rows)[ii];
          nvect.z = nmap.ptr(jj + 2 * rows)[ii];
          const float theta = acosf(fabsf(nvect.z));
          if (theta > PI / 2.5f) // ignore normal angle larger than 72 degrees
            continue;
          const float sigmaL = 0.8f + 0.035f * theta / (PI / 2.0f - theta);
          const float sigmaz = 0.0012f +
                               0.0019f * (Dp_scaled - 0.4f) * (Dp_scaled - 0.4f) +
                               0.0001f / sqrtf(Dp_scaled) * theta * theta /
                                   (PI / 2.f - theta) / (PI / 2.f - theta);
          const float sdf = Dp_scaled - sqrtf(v_g_z * v_g_z + v_g_part_norm);

          if ((Dp_scaled != 0) && (dz < 3.f * sigmaz) &&
              (sdf > -6.f * sigmaz)) // meters
          {
            // new form of truncated signed distance function
            float tsdf = sqrtf(1.f - __expf(-2.f / PI * sdf * sdf / (sigmaz * sigmaz)));
            if (sdf < 0)
              tsdf = -tsdf;

            const float Wrk =
                0.0012f / sigmaz * (0.4f / Dp_scaled) * (0.4f / Dp_scaled) *
                __expf(-du * du / (2.f * sigmaL * sigmaL) -
                       dz * dz / (2.f * sigmaz * sigmaz)); // range from ~0 to 1.0

            tsdf_prev = (tsdf_prev * weight_prev + Wrk * tsdf) / (weight_prev + Wrk);
            weight_prev = min(weight_prev + Wrk, float(2 * Tsdf::MAX_WEIGHT));
          }
        }
      }
      pack_tsdf(tsdf_prev, weight_prev, *pos);
    }
  } // for(int z = 0; z < VOLUME_Z; ++z)
} // __global__

//================================================zivid noise
// fucntions===================
__global__ void
    tsdf23_axial_zivid /*(const PtrStepSz<float> depthScaled, PtrStep<short2> volume,
                        const float tranc_dist, const Mat33 Rcurr_inv, const float3
                        tcurr, const Intr intr, const float3 cell_size)*/

    (const PtrStepSz<float> depthScaled,
     const PtrStep<float> nmap,
     int rows,
     int cols, /*added*/
     PtrStep<short2> volume,
     const Mat33 Rcurr_inv,
     const float3 tcurr,
     const Intr intr,
     const float3 cell_size,
     const noise_data* d_data,
     const float* d_dis_ptr,
     const float* d_ang_ptr,
     int data_size,
     int dis_size,
     int ang_size)
{
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;

  if (x >= VOLUME_X || y >= VOLUME_Y)
    return;

  float v_g_x = (x + 0.5f) * cell_size.x - tcurr.x;
  float v_g_y = (y + 0.5f) * cell_size.y - tcurr.y;
  float v_g_z = (0 + 0.5f) * cell_size.z - tcurr.z;

  float v_g_part_norm = v_g_x * v_g_x + v_g_y * v_g_y;

  float v_x = (Rcurr_inv.data[0].x * v_g_x + Rcurr_inv.data[0].y * v_g_y +
               Rcurr_inv.data[0].z * v_g_z) *
              intr.fx;
  float v_y = (Rcurr_inv.data[1].x * v_g_x + Rcurr_inv.data[1].y * v_g_y +
               Rcurr_inv.data[1].z * v_g_z) *
              intr.fy;
  float v_z = (Rcurr_inv.data[2].x * v_g_x + Rcurr_inv.data[2].y * v_g_y +
               Rcurr_inv.data[2].z * v_g_z);

  float z_scaled = 0;

  float Rcurr_inv_0_z_scaled = Rcurr_inv.data[0].z * cell_size.z * intr.fx;
  float Rcurr_inv_1_z_scaled = Rcurr_inv.data[1].z * cell_size.z * intr.fy;

  short2* pos = volume.ptr(y) + x;
  int elem_step = volume.step * VOLUME_Y / sizeof(short2);

  // #pragma unroll
  for (int z = 0; z < VOLUME_Z; ++z,
           v_g_z += cell_size.z,
           z_scaled += cell_size.z,
           v_x += Rcurr_inv_0_z_scaled,
           v_y += Rcurr_inv_1_z_scaled,
           pos += elem_step) {
    float inv_z = 1.0f / (v_z + Rcurr_inv.data[2].z * z_scaled);
    if (inv_z < 0)
      continue;

    // project to current cam
    int2 coo = {__float2int_rn(v_x * inv_z + intr.cx),
                __float2int_rn(v_y * inv_z + intr.cy)};

    if (coo.x >= 0 && coo.y >= 0 && coo.x < depthScaled.cols &&
        coo.y < depthScaled.rows) // 6
    {
      float Dp_scaled = depthScaled.ptr(coo.y)[coo.x]; // meters

      // calculate sigma_L and sigma_z
      float3 nvect;
      nvect.x = nmap.ptr(coo.y)[coo.x];
      nvect.y = nmap.ptr(coo.y + rows)[coo.x];
      nvect.z = nmap.ptr(coo.y + 2 * rows)[coo.x];
      const float theta = acosf(fabsf(nvect.z));
      if ((theta * 180.f / PI > 70) || (Dp_scaled < 0.037) || (Dp_scaled > 1.070))
        continue; // this requires extrapolation
      float sigmaz = sig_z_default_value;
      sigmaz = compute_sigma2(Dp_scaled,
                              theta,
                              sig_z_default_value,
                              sigmaz,
                              d_data,
                              data_size,
                              factor_,
                              d_dis_ptr,
                              d_ang_ptr,
                              dis_size,
                              ang_size);
      // float sigmaz_ = compute_sigma2_(
      //     Dp_scaled, theta, sig_z_default_value, sigmaz, d_data, data_size, factor_);

      /* if (sigmaz >2.0592945567003293 / 1000.f)
       {
           sigmaz = 2.0592945567003293 / 1000.f;
       }
       */
      // printf("factor=%f,sig_z_default_value=%f\n", factor_, sig_z_default_value);
      if (sigmaz != sig_z_default_value) // if sigma value has  been updated from
                                         // initial bad value
      {
        float sdf = Dp_scaled - sqrtf(v_g_z * v_g_z + v_g_part_norm);

        if (Dp_scaled != 0 && sdf >= -6.f * sigmaz) // meters
        {
          float tsdf = sqrtf(1.0f - __expf(-2.f / PI * sdf * sdf / (100 * sigmaz * sigmaz)));
          if (sdf < 0)
            tsdf = -tsdf;
          // read and unpack
          float tsdf_prev, weight_prev;
          unpack_tsdf(*pos, tsdf_prev, weight_prev);

          // const float Wrk = 0.0012f / sigmaz * (0.4f / Dp_scaled) * (0.4f /
          // Dp_scaled); // range from ~0 to 1.0 const float Wrk = (1 - ((sigmaz -
          // (0.03586010878107277f / 1000.f) * factor_) / ((1.8498489058744052 / 1000.f)
          // * factor_ - (0.03586010878107277f / 1000.f) * factor_))) +
          // __expf(-Dp_scaled * Dp_scaled / (2.f * sigmaz * sigmaz));

          //  const float Wrk = (1 - ((sigmaz - (0.03586010878107277f / 1000.f) *
          //  factor_) / ((2.1527348301713722 / 1000.f) * factor_ -
          //  (0.03586010878107277f / 1000.f) * factor_)));

          //
          //  const float Wrk = 1 / sigmaz;
          //  const float Wrk = 1.f;
          const float Wrk =
              0.00003545f / sigmaz * (0.4f / Dp_scaled) * (0.4f / Dp_scaled);
          ; // range from ~0 to 1.0
          //  weight saved is with kine                        _ _cexpf(-dz * dz / (2.f
          //  * sigmaz * sigmaz))t noise tailored for zivid.it was worst const float Wrk
          //  = 1 - ((sigmaz - (0.03586f / 1000.f) * factor_) / ((1.5 / 1000.f) *
          //  factor_ - (0.03586f / 1000.f) * factor_)); weight2 saved is with 2.15 max.
          //  it is worse than when max is 1.84 weight3 saved is with 1.5 max with
          //  condition that if sigmaz>1.5/1000, set it to 1.5.1000.(bad) otherwise for
          //  weight3 setting we will have flying points appearing in frame around 45
          //  and more weight4, same as above but with float(Tsdf::MAX_WEIGHT)) instead
          //  of float(2*Tsdf::MAX_WEIGHT)) weight5 is with max 1.8498 but with
          //  float(Tsdf::MAX_WEIGHT)) instead of float(2*Tsdf::MAX_WEIGHT))
          float tsdf_new = (tsdf_prev * weight_prev + Wrk * tsdf) / (weight_prev + Wrk);
          float weight_new = min(weight_prev + Wrk, float(2 * Tsdf::MAX_WEIGHT));

          pack_tsdf(tsdf_new, weight_new, *pos);
        }
      }
    }
  } // for(int z = 0; z < VOLUME_Z; ++z)
} // __global__

//=======================================================both axial and lateral
__global__ void
    tsdf23_both_noise_zivid /*(const PtrStepSz<float> depthScaled, PtrStep<short2>
                              volume, const float tranc_dist, const Mat33 Rcurr_inv,
                              const float3 tcurr, const Intr intr, const float3
                              cell_size)*/
    (const PtrStepSz<float> depthScaled,
     const PtrStep<float> nmap,
     int rows,
     int cols, /*added*/
     PtrStep<short2> volume,
     const Mat33 Rcurr_inv,
     const float3 tcurr,
     const Intr intr,
     const float3 cell_size,
     const noise_data* d_data,
     const float* d_dis_ptr,
     const float* d_ang_ptr,
     int data_size,
     int dis_size,
     int ang_size)
{
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;

  if (x >= VOLUME_X || y >= VOLUME_Y)
    return;

  float v_g_x = (x + 0.5f) * cell_size.x - tcurr.x;
  float v_g_y = (y + 0.5f) * cell_size.y - tcurr.y;
  float v_g_z = (0 + 0.5f) * cell_size.z - tcurr.z;

  float v_g_part_norm = v_g_x * v_g_x + v_g_y * v_g_y;

  float v_x = (Rcurr_inv.data[0].x * v_g_x + Rcurr_inv.data[0].y * v_g_y +
               Rcurr_inv.data[0].z * v_g_z) *
              intr.fx;
  float v_y = (Rcurr_inv.data[1].x * v_g_x + Rcurr_inv.data[1].y * v_g_y +
               Rcurr_inv.data[1].z * v_g_z) *
              intr.fy;
  float v_z = (Rcurr_inv.data[2].x * v_g_x + Rcurr_inv.data[2].y * v_g_y +
               Rcurr_inv.data[2].z * v_g_z);

  float z_scaled = 0;

  float Rcurr_inv_0_z_scaled = Rcurr_inv.data[0].z * cell_size.z * intr.fx;
  float Rcurr_inv_1_z_scaled = Rcurr_inv.data[1].z * cell_size.z * intr.fy;

  // float tranc_dist_inv = 1.0f / tranc_dist;

  short2* pos = volume.ptr(y) + x;
  int elem_step = volume.step * VOLUME_Y / sizeof(short2);

  // #pragma unroll
  for (int z = 0; z < VOLUME_Z; ++z,
           v_g_z += cell_size.z,
           z_scaled += cell_size.z,
           v_x += Rcurr_inv_0_z_scaled,
           v_y += Rcurr_inv_1_z_scaled,
           pos += elem_step) {
    float inv_z = 1.0f / (v_z + Rcurr_inv.data[2].z * z_scaled);
    if (inv_z < 0)
      continue;

    // project to current cam
    int2 coo = {__float2int_rn(v_x * inv_z + intr.cx),
                __float2int_rn(v_y * inv_z + intr.cy)};
    // and subpixel
    float2 dcoo = {(coo.x - (v_x * inv_z + intr.cx)),
                   (coo.y - (v_y * inv_z + intr.cy))};

    if (coo.x >= 1 && coo.y >= 1 && coo.x < depthScaled.cols - 1 &&
        coo.y < depthScaled.rows - 1) {
      float Dp_scaled_mid = depthScaled.ptr(coo.y)[coo.x]; // meters
      if (Dp_scaled_mid == 0)
        continue;

      // read and unpack
      float tsdf_prev, weight_prev;
      unpack_tsdf(*pos, tsdf_prev, weight_prev);
      // loop within 3x3 pixel area around the current position

#pragma unroll
      for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
          const int ii = coo.x + i;
          const int jj = coo.y + j;
          const float Dp_scaled = depthScaled.ptr(jj)[ii]; // meters
          const float dz = fabs(Dp_scaled - Dp_scaled_mid);
          const float du =
              sqrtf((i + dcoo.x) * (i + dcoo.x) + (j + dcoo.y) * (j + dcoo.y));

          // calculate sigma_L and sigma_z
          float3 nvect;
          nvect.x = nmap.ptr(jj)[ii];
          nvect.y = nmap.ptr(jj + rows)[ii];
          nvect.z = nmap.ptr(jj + 2 * rows)[ii];
          const float theta = acosf(fabsf(nvect.z));
          // if (theta > PI / 2.5f) // ignore normal angle larger than 72 degrees
          if ((theta < 0) || (theta * 180.f / PI > 70) || (Dp_scaled < 0.037) ||
              (Dp_scaled > 1.070))
            continue;
          // const float sigmaL = 0.8f + 0.035f * theta / (PI / 2.0f - theta);
          // const float sigmaz = 0.0012f + 0.0019f * (Dp_scaled - 0.4f) * (Dp_scaled -
          // 0.4f) + 0.0001f / sqrtf(Dp_scaled) * theta * theta / (PI / 2.f - theta) /
          // (PI / 2.f - theta);
          float sigmaz = sig_z_default_value;
          sigmaz = compute_sigma2(Dp_scaled,
                                  theta,
                                  sig_z_default_value,
                                  sigmaz,
                                  d_data,
                                  data_size,
                                  factor_,
                                  d_dis_ptr,
                                  d_ang_ptr,
                                  dis_size,
                                  ang_size);
          // const float sigmaL = 0.001431 * Dp_scaled;
          const float sigmaL =
              0.00041 * (KINFU_DEFAULT_DEPTH_FOCAL_X + KINFU_DEFAULT_DEPTH_FOCAL_Y) / 2;
          const float sdf = Dp_scaled - sqrtf(v_g_z * v_g_z + v_g_part_norm);

          if ((Dp_scaled != 0) && (dz < 3.f * sigmaz) &&
              (sdf > -6.f * sigmaz)) // meters
          {
            // new form of truncated signed distance function
            float tsdf = sqrtf(1.f - __expf(-2.f / PI * sdf * sdf / (sigmaz * sigmaz)));
            if (sdf < 0)
              tsdf = -tsdf;

            // range from ~0 to 1.0
            // const float Wrk = 0.0012f / sigmaz * (0.4f / Dp_scaled) *
            //                   (0.4f / Dp_scaled) *
            //                   __expf(-du * du / (2.f * sigmaL * sigmaL) -
            //                          dz * dz / (2.f * sigmaz * sigmaz));
            const float Wrk = 0.00003545f / sigmaz * (0.4f / Dp_scaled) *
                              (0.4f / Dp_scaled) *
                              __expf(-du * du / (2.f * sigmaL * sigmaL) -
                                     dz * dz / (2.f * sigmaz * sigmaz));

            // const float Wrk =
            //     1 - (((sigmaz - (0.03586f / 1000.f) * factor_) /
            //           ((1.8498 / 1000.f) * factor_ - (0.03586f / 1000.f) * factor_))
            //           +
            //          ((sigmaL - 0.001431 * (0.03586f / 1000.f)) /
            //           ((0.001431 * 1.070) - (0.03586f / 1000.f))));

            tsdf_prev = (tsdf_prev * weight_prev + Wrk * tsdf) / (weight_prev + Wrk);
            weight_prev = min(weight_prev + Wrk, float(2 * Tsdf::MAX_WEIGHT));
          }
        }
      }
      pack_tsdf(tsdf_prev, weight_prev, *pos);
    }
  } // for(int z = 0; z < VOLUME_Z; ++z)
} // __global__

//===============================================================================================
__global__ void
tsdf23normal_hack(const PtrStepSz<float> depthScaled,
                  PtrStep<short2> volume,
                  const float tranc_dist,
                  const Mat33 Rcurr_inv,
                  const float3 tcurr,
                  const Intr intr,
                  const float3 cell_size)
{
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;

  if (x >= VOLUME_X || y >= VOLUME_Y)
    return;

  const float v_g_x = (x + 0.5f) * cell_size.x - tcurr.x;
  const float v_g_y = (y + 0.5f) * cell_size.y - tcurr.y;
  float v_g_z = (0 + 0.5f) * cell_size.z - tcurr.z;

  float v_g_part_norm = v_g_x * v_g_x + v_g_y * v_g_y;

  float v_x = (Rcurr_inv.data[0].x * v_g_x + Rcurr_inv.data[0].y * v_g_y +
               Rcurr_inv.data[0].z * v_g_z) *
              intr.fx;
  float v_y = (Rcurr_inv.data[1].x * v_g_x + Rcurr_inv.data[1].y * v_g_y +
               Rcurr_inv.data[1].z * v_g_z) *
              intr.fy;
  float v_z = (Rcurr_inv.data[2].x * v_g_x + Rcurr_inv.data[2].y * v_g_y +
               Rcurr_inv.data[2].z * v_g_z);

  float z_scaled = 0;

  float Rcurr_inv_0_z_scaled = Rcurr_inv.data[0].z * cell_size.z * intr.fx;
  float Rcurr_inv_1_z_scaled = Rcurr_inv.data[1].z * cell_size.z * intr.fy;

  float tranc_dist_inv = 1.0f / tranc_dist;

  short2* pos = volume.ptr(y) + x;
  int elem_step = volume.step * VOLUME_Y / sizeof(short2);

  // #pragma unroll
  for (int z = 0; z < VOLUME_Z; ++z,
           v_g_z += cell_size.z,
           z_scaled += cell_size.z,
           v_x += Rcurr_inv_0_z_scaled,
           v_y += Rcurr_inv_1_z_scaled,
           pos += elem_step) {
    float inv_z = 1.0f / (v_z + Rcurr_inv.data[2].z * z_scaled);
    if (inv_z < 0)
      continue;

    // project to current cam
    int2 coo = {__float2int_rn(v_x * inv_z + intr.cx),
                __float2int_rn(v_y * inv_z + intr.cy)};

    if (coo.x >= 0 && coo.y >= 0 && coo.x < depthScaled.cols &&
        coo.y < depthScaled.rows) // 6
    {
      float Dp_scaled = depthScaled.ptr(coo.y)[coo.x]; // meters

      float sdf = Dp_scaled - sqrtf(v_g_z * v_g_z + v_g_part_norm);

      if (Dp_scaled != 0 && sdf >= -tranc_dist) // meters
      {
        float tsdf = fmin(1.0f, sdf * tranc_dist_inv);

        bool integrate = true;
        if ((x > 0 && x < VOLUME_X - 2) && (y > 0 && y < VOLUME_Y - 2) &&
            (z > 0 && z < VOLUME_Z - 2)) {
          const float qnan = std::numeric_limits<float>::quiet_NaN();
          float3 normal = make_float3(qnan, qnan, qnan);

          float Fn, Fp;
          float Wn = 0, Wp = 0;
          unpack_tsdf(*(pos + elem_step), Fn, Wn);
          unpack_tsdf(*(pos - elem_step), Fp, Wp);

          if (Wn > 16 && Wp > 16)
            normal.z = (Fn - Fp) / cell_size.z;

          unpack_tsdf(*(pos + volume.step / sizeof(short2)), Fn, Wn);
          unpack_tsdf(*(pos - volume.step / sizeof(short2)), Fp, Wp);

          if (Wn > 16 && Wp > 16)
            normal.y = (Fn - Fp) / cell_size.y;

          unpack_tsdf(*(pos + 1), Fn, Wn);
          unpack_tsdf(*(pos - 1), Fp, Wp);

          if (Wn > 16 && Wp > 16)
            normal.x = (Fn - Fp) / cell_size.x;

          if (normal.x != qnan && normal.y != qnan && normal.z != qnan) {
            float norm2 = dot(normal, normal);
            if (norm2 >= 1e-10) {
              normal *= rsqrt(norm2);

              float nt = v_g_x * normal.x + v_g_y * normal.y + v_g_z * normal.z;
              float cosine = nt * rsqrt(v_g_x * v_g_x + v_g_y * v_g_y + v_g_z * v_g_z);

              if (cosine < 0.5)
                integrate = false;
            }
          }
        }

        if (integrate) {
          // read and unpack
          float tsdf_prev;
          float weight_prev;
          unpack_tsdf(*pos, tsdf_prev, weight_prev);

          const int Wrk = 1;

          float tsdf_new = (tsdf_prev * weight_prev + Wrk * tsdf) / (weight_prev + Wrk);
          float weight_new = min(weight_prev + Wrk, float(Tsdf::MAX_WEIGHT));

          pack_tsdf(tsdf_new, weight_new, *pos);
        }
      }
    }
  } // for(int z = 0; z < VOLUME_Z; ++z)
} // __global__
} // namespace device
} // namespace pcl
//=============================================================
void
pcl::device::integrateTsdfVolume(const PtrStepSz<ushort>& depth,
                                 const MapArr& nmap,
                                 const Intr& intr,
                                 const float3& volume_size,
                                 const Mat33& Rcurr_inv,
                                 const float3& tcurr,
                                 float tranc_dist,
                                 PtrStep<short2> volume,
                                 DeviceArray2D<float>& depthScaled,
                                 int noise_components)
{
  depthScaled.create(depth.rows, depth.cols);

  dim3 block_scale(32, 8);
  dim3 grid_scale(divUp(depth.cols, block_scale.x), divUp(depth.rows, block_scale.y));

  // scales depth along ray and converts mm -> meters.
  scaleDepth<<<grid_scale, block_scale>>>(depth, depthScaled, intr);
  cudaSafeCall(cudaGetLastError());

  float3 cell_size;
  cell_size.x = volume_size.x / VOLUME_X;
  cell_size.y = volume_size.y / VOLUME_Y;
  cell_size.z = volume_size.z / VOLUME_Z;

  // dim3 block(Tsdf::CTA_SIZE_X, Tsdf::CTA_SIZE_Y);
  dim3 block(16, 16);
  dim3 grid(divUp(VOLUME_X, block.x), divUp(VOLUME_Y, block.y));

  int cols = nmap.cols();
  int rows = nmap.rows() / 3;

  if (noise_components == 0) {
    std::cout << "no noise" << std::endl;
    tsdf23<<<grid, block>>>(
        depthScaled, volume, tranc_dist, Rcurr_inv, tcurr, intr, cell_size);
  }
  if (apply_zivid_noise == 1) {
    if (noise_components == 1) {
      std::cout << "zivid axial" << std::endl;
      tsdf23_axial_zivid<<<grid, block>>>(depthScaled,
                                          nmap,
                                          rows,
                                          cols,
                                          volume,
                                          Rcurr_inv,
                                          tcurr,
                                          intr,
                                          cell_size,
                                          d_data_ptr,
                                          d_dis_ptr,
                                          d_ang_ptr,
                                          data_size,
                                          dis_size,
                                          ang_size);
    }

    if (noise_components == 2) {
      std::cout << "zivid both" << std::endl;
      tsdf23_both_noise_zivid<<<grid, block>>>(depthScaled,
                                               nmap,
                                               rows,
                                               cols,
                                               volume,
                                               Rcurr_inv,
                                               tcurr,
                                               intr,
                                               cell_size,
                                               d_data_ptr,
                                               d_dis_ptr,
                                               d_ang_ptr,
                                               data_size,
                                               dis_size,
                                               ang_size);
    }
  }
  if (apply_zivid_noise == 0) {
    if (noise_components == 1) {
      std::cout << "kinect axial" << std::endl;
      tsdf23_axial_kinect<<<grid, block>>>(
          depthScaled, nmap, rows, cols, volume, Rcurr_inv, tcurr, intr, cell_size);
    }
    if (noise_components == 2) {
      std::cout << "kinect both" << std::endl;
      tsdf23_both_noise_kinect<<<grid, block>>>(
          depthScaled, nmap, rows, cols, volume, Rcurr_inv, tcurr, intr, cell_size);
    }
  }

  // tsdf23normal_hack<<<grid, block>>>(depthScaled, volume, tranc_dist, Rcurr_inv,
  // tcurr, intr, cell_size);

  cudaSafeCall(cudaGetLastError());
  cudaSafeCall(cudaDeviceSynchronize());
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::device::integrateTsdfVolume(const PtrStepSz<ushort>& depth,
                                 const Intr& intr,
                                 const float3& volume_size,
                                 const Mat33& Rcurr_inv,
                                 const float3& tcurr,
                                 float tranc_dist,
                                 PtrStep<short2> volume,
                                 DeviceArray2D<float>& depthScaled)
{
  depthScaled.create(depth.rows, depth.cols);

  dim3 block_scale(32, 8);
  dim3 grid_scale(divUp(depth.cols, block_scale.x), divUp(depth.rows, block_scale.y));

  // scales depth along ray and converts mm -> meters.
  scaleDepth<<<grid_scale, block_scale>>>(depth, depthScaled, intr);
  cudaSafeCall(cudaGetLastError());

  float3 cell_size;
  cell_size.x = volume_size.x / VOLUME_X;
  cell_size.y = volume_size.y / VOLUME_Y;
  cell_size.z = volume_size.z / VOLUME_Z;

  // dim3 block(Tsdf::CTA_SIZE_X, Tsdf::CTA_SIZE_Y);
  dim3 block(16, 16);
  dim3 grid(divUp(VOLUME_X, block.x), divUp(VOLUME_Y, block.y));

  tsdf23<<<grid, block>>>(
      depthScaled, volume, tranc_dist, Rcurr_inv, tcurr, intr, cell_size);
  // tsdf23normal_hack<<<grid, block>>>(depthScaled, volume, tranc_dist, Rcurr_inv,
  // tcurr, intr, cell_size);

  cudaSafeCall(cudaGetLastError());
  cudaSafeCall(cudaDeviceSynchronize());
}

#endif /* TSDF_VOLUME_CU */
