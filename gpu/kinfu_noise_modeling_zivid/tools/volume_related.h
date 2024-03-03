#ifndef VOLUME_RELATED_H
#define VOLUME_RELATED_H

#include <pcl/gpu/kinfu/camera_type.h>

#include <iostream>

// which camera settings
// #define current_camera_type KINECT_CAMERA

#define voxel_tmp 1024 // Original 512

#if current_camera_type == ZIVID_CAMERA
#define volume_size_tmp 500.f              // 950.f. Used till February 14, 2024
#define depth_volume_size_tmp 1000.f / 1000 // 100.f // 150.f
// it is strange that default volume size works even if we put 0.f
#define DEFAULT_VOLUME_SIZE_X 900.f / 1000 // volume_size_tmp - 200 //750.f
#define DEFAULT_VOLUME_SIZE_Y 700.f / 1000 // volume_size_tmp - 50  //900.f
#define DEFAULT_VOLUME_SIZE_Z depth_volume_size_tmp
#define trunc_dist_m volume_size_tmp / (1000.f * 10000.f)
// 0.009f gives fatter face for noise 0 even//volume_size_tmp / (1000.f * 600.f)
#define min_trunc_dist_ 0.2f
#define apply_zivid_noise_ 1

#else // KINECT_CAMERA
#define volume_size_tmp 3000.f
#define depth_volume_size_tmp volume_size_tmp
#define DEFAULT_VOLUME_SIZE_X volume_size_tmp
#define DEFAULT_VOLUME_SIZE_Y volume_size_tmp
#define DEFAULT_VOLUME_SIZE_Z volume_size_tmp
#define trunc_dist_m 0.03f
#define min_trunc_dist_ 0.2f
// 30.f  used in tools / record_tsdfvolume_zivid.cpp. currently not related
#define apply_zivid_noise_ 0
#endif

// Define settings for Kinect camera
#if current_camera_type == KINECT_CAMERA
#define volume_size_tmp 3000.f
#define depth_volume_size_tmp volume_size_tmp
#define DEFAULT_VOLUME_SIZE_X volume_size_tmp
#define DEFAULT_VOLUME_SIZE_Y volume_size_tmp
#define DEFAULT_VOLUME_SIZE_Z volume_size_tmp
#define trunc_dist_m 0.03f
#define min_trunc_dist_ 0.2f
// 30.f  used in tools / record_tsdfvolume_zivid.cpp. currently not related
#define apply_zivid_noise_ 0
#endif

// Define common settings
#define volume_size_tmp_meters volume_size_tmp / 1000.f
#define neg_volume_size_tmp_meters -volume_size_tmp / 1000.f

#define DEFAULT_GRID_RES_X voxel_tmp
#define DEFAULT_GRID_RES_Y voxel_tmp
#define DEFAULT_GRID_RES_Z voxel_tmp

constexpr int volume_xyz = voxel_tmp;
#define default_threshold_dist_ 0.1f

const float Pi = 3.14159265f;
const float PI = 3.14159265f;
#define sigma_color_ 30
#define sigma_space_ 4.5f
#define theta_mean_ Pi / 6.0f
//====================================
#define sig_z_default_value_ 0.1f
#define noise_scaling_factor 1.0
#define reserve_cloud_points 1000000000 // used in tsdf_volume.cpp

#endif /* VOLUME_RELATED_H */
