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
      pack_tsdf(0.f, 0.f, *pos);
  }
}
} // namespace device
} // namespace pcl

void
pcl::device::initVolume(PtrStep<short2> volume)
{
  dim3 block(32, 16);
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
            float weight_new = min(weight_prev + Wrk, float(Tsdf::MAX_WEIGHT));

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

__global__ void
tsdf23_axial_noise(const PtrStepSz<float> depthScaled,
                   const PtrStep<float> nmap,
                   int rows,
                   int cols,
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
        coo.y < depthScaled.rows) {
      float Dp_scaled = depthScaled.ptr(coo.y)[coo.x]; // meters

      // calculate sigma_l and sigma_z
      float3 nvect;
      nvect.x = nmap.ptr(coo.y)[coo.x];
      nvect.y = nmap.ptr(coo.y + rows)[coo.x];
      nvect.z = nmap.ptr(coo.y + 2 * rows)[coo.x];
      const float theta = acosf(fabsf(nvect.z));
      if (theta > PI / 2.5f) // ignore normal angle larger than 72 degrees
        continue;
      // const float sigmaL = 0.8f + 0.035f * theta / (PI / 2.0f - theta);
      const float sigmaz = 0.0012f + 0.0019f * (Dp_scaled - 0.4f) * (Dp_scaled - 0.4f) +
                           0.0001f / sqrtf(Dp_scaled) * theta * theta /
                               (PI / 2.f - theta) / (PI / 2.f - theta);

      float sdf = Dp_scaled - sqrtf(v_g_z * v_g_z + v_g_part_norm);

      if (Dp_scaled != 0 && sdf > -6.f * sigmaz) // meters
      {
        float tsdf = sqrtf(1.0f - __expf(-2.f / PI * sdf * sdf / (sigmaz * sigmaz)));
        if (sdf < 0)
          tsdf = -tsdf;

        // read and uppack
        float tsdf_prev, weight_prev;
        unpack_tsdf(*pos, tsdf_prev, weight_prev);

        // assume z_min = 0.4
        const float Wrk = 0.0012f / sigmaz * (0.4f / Dp_scaled) *
                          (0.4f / Dp_scaled); // range from ~0 to 1.0

        float tsdf_new = (tsdf_prev * weight_prev + Wrk * tsdf) / (weight_prev + Wrk);
        float weight_new = min(weight_prev + Wrk, float(2 * Tsdf::MAX_WEIGHT));

        pack_tsdf(tsdf_new, weight_new, *pos);
      }
    }
  }
}

__global__ void
tsdf23_all_noise(const PtrStepSz<float> depthScaled,
                 const PtrStep<float> nmap,
                 int rows,
                 int cols,
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

          // signed distance function
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
}

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
          float Wn = 0.f, Wp = 0.f;
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::device::integrateWeightedTsdfVolume(const PtrStepSz<ushort>& depth,
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
  if (noise_components == 0)
    tsdf23<<<grid, block>>>(
        depthScaled, volume, tranc_dist, Rcurr_inv, tcurr, intr, cell_size);
  else if (noise_components == 1)
    tsdf23_axial_noise<<<grid, block>>>(
        depthScaled, nmap, rows, cols, volume, Rcurr_inv, tcurr, intr, cell_size);
  else // (noise_components = 2)
    tsdf23_all_noise<<<grid, block>>>(
        depthScaled, nmap, rows, cols, volume, Rcurr_inv, tcurr, intr, cell_size);

  cudaSafeCall(cudaGetLastError());
  cudaSafeCall(cudaDeviceSynchronize());
}
