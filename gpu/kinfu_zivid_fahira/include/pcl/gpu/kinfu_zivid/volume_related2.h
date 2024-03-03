#ifndef VOLUME_RELATED_H
#define VOLUME_RELATED_H

#include <iostream>

//=========================zivid settings=====================================
//-------------------------volume related-------------------used in
// tools/tsdf_volume.h, src/estimate_transform.cu, src/tsdf_volume.cpp
#define voxel_tmp 1024 // 0riginal 512
#define volume_size_tmp 900.f
// 950.f used till feb14, 2024 //850. f does not work and truncates the back
// side of shiva
// original 3000, used in tools/record_tsdfvolume_zivid.cpp
#define volume_size_tmp_meters volume_size_tmp / 1000.f
//#define neg_volume_size_tmp_meters -volume_size_tmp / 1000.f  if
// volume_size_tmp is 900.f
#define neg_volume_size_tmp_meters -volume_size_tmp / 1000.f
//-3.f//-volume_size_tmp/1000 //used in

// src/raycaster.cpp in generateSceneView. not sure
// if valid as original value is -3.f
//------------------------------------------------------
// following default_ defines are used in tools/tsdf_volume.h
// used till feb14, 2024
/*
#define DEFAULT_GRID_RES_X 2 * voxel_tmp
#define DEFAULT_GRID_RES_Y 2 * voxel_tmp
#define DEFAULT_GRID_RES_Z 2 * voxel_tmp
*/
//-----------------------------------
#define DEFAULT_GRID_RES_X 4 * voxel_tmp // pcl::device::VOLUME_X ( and _Y, _Z)
#define DEFAULT_GRID_RES_Y 4 * voxel_tmp
#define DEFAULT_GRID_RES_Z 4 * voxel_tmp
#define depth_volume_size_tmp 150.f // 100.f // 150.f
// it is strange that default volume size works even if we put 0.f
#define DEFAULT_VOLUME_SIZE_X 750.f // volume_size_tmp - 200 //750.f
#define DEFAULT_VOLUME_SIZE_Y 900.f // volume_size_tmp - 50  //900.f
#define DEFAULT_VOLUME_SIZE_Z depth_volume_size_tmp

//---------------------------used in src/internal.h
constexpr int volume_xyz = voxel_tmp;
//-----------------
// control the lower limit of depth we want to see, larger the value, the lower
// depths will not be rendered
// volume_size_tmp_meters/10 works better compared to
// //volume_size_tmp_meters/100
//=======================
// have used this till feb8,2024
//#define trunc_dist_m volume_size_tmp_meters / 600.f
//============================
#define trunc_dist_m volume_size_tmp_meters / 600.f
// 0.95 / 600 // volume_size_tmp_meters / 600.f
// 950.f/(1000.f*300.f)//volume_size_tmp_meters/300 // 0.02 creates pb
// and truncates certain level of depth // volume_size_tmp_meters/100
// //0.1f //volume_size_tmp_meters/10 //0.00001//0.01
// //volume_size_tmp_meters/100 //src/tsdf_volume.cpp as
// setTsdfTruncDist
/*
//======================kinect settings=============================
#define voxel_tmp 1024 // 0riginal 512
#define volume_size_tmp 3000.f
#define volume_size_tmp_meters volume_size_tmp / 1000.f
#define neg_volume_size_tmp_meters -volume_size_tmp / 1000.f

//-----------------------------------
#define DEFAULT_GRID_RES_X 4 * voxel_tmp // pcl::device::VOLUME_X ( and _Y, _Z)
#define DEFAULT_GRID_RES_Y 4 * voxel_tmp
#define DEFAULT_GRID_RES_Z 4 * voxel_tmp
#define depth_volume_size_tmp volume_size_tmp
#define DEFAULT_VOLUME_SIZE_X volume_size_tmp
#define DEFAULT_VOLUME_SIZE_Y volume_size_tmp
#define DEFAULT_VOLUME_SIZE_Z volume_size_tmp

//---------------------------used in src/internal.h
constexpr int volume_xyz = voxel_tmp;
#define trunc_dist_m 0.03f
*/
//=========================common settings=====================================
#define default_threshold_dist_ 0.1f
// 0.1f default used in kinfu_z.cpp for icp filtering correspondences

#define min_trunc_dist_                                                        \
  0.1f; // volume_size_tmp/100.f    //30.0f; used in
        // tools/record_tsdfvolume_zivid.cpp

// 0.2f used till feb 12,2024                      //
//
const float Pi = 3.14159265f;
const float PI = 3.14159265f;
//-------------------------------bilateral_filter
#define sigma_color_ 30   // 0.00003 // 30 * Pi / 180.f
#define sigma_space_ 4.5f // 10.5 // 0.008f // 4 .5 in pixels in bilateral pyro
#define theta_mean_ Pi / 6.0f
//====================================
#define noise_scaling_factor 1.f
#define apply_zivid_noise_ 1
//-------------------------------
#define reserve_cloud_points 1000000000 // used in tsdf_volume.cpp
#endif                                  /* VOLUME_RELATED_H */
