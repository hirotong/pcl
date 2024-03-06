#ifndef VOLUME_RELATED_H
#define VOLUME_RELATED_H

#include "camera_type.h"

#include <iostream>

#define voxel_tmp 1024 // Original 512

#if current_camera_type == ZIVID_CAMERA
#define volume_size_tmp 300 // 450.f till feb26 // 950.f. Used till February 14,
                            // 2024// then 900
#if (getting_correct_pose == 1)
#define volume_size_tmp 900
#endif
// it is strange that default volume size works even if we put 0.f
#define DEFAULT_VOLUME_SIZE_X volume_size_tmp
#define DEFAULT_VOLUME_SIZE_Y volume_size_tmp
#define DEFAULT_VOLUME_SIZE_Z volume_size_tmp
/*// used for upto 450 volume size
#define depth_volume_size_tmp 150.f // 100.f // 150.f
#define DEFAULT_VOLUME_SIZE_X 750.f // volume_size_tmp - 200 //750.f
#define DEFAULT_VOLUME_SIZE_Y 900.f // volume_size_tmp - 50  //900.f
#define DEFAULT_VOLUME_SIZE_Z depth_volume_size_tmp
*/
#define trunc_dist_m volume_size_tmp / (1000.f * 10000.f) // used till feb for
// above used till 26 feb upto 450 volume size. this is default_tranc_dist in
// tsdf inverse of it is used to multiply with sdf to get tsdf in an expression
//#define trunc_dist_m volume_size_tmp / (1000.f * 1000.f) // used for volume
// 300 with shiva
// 0.009f gives fatter face for noise 0 even//volume_size_tmp / (1000.f * 600.f)
#define min_trunc_dist_ 0.2f
#define apply_zivid_noise_ 1
#endif

// Define settings for Kinect camera =========================
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

// Define common settings ==========================
#define volume_size_tmp_meters (volume_size_tmp / 1000.f)
#define neg_volume_size_tmp_meters (-volume_size_tmp / 1000.f)

#define DEFAULT_GRID_RES_X voxel_tmp
#define DEFAULT_GRID_RES_Y voxel_tmp
#define DEFAULT_GRID_RES_Z voxel_tmp

constexpr int volume_xyz = voxel_tmp;
#define default_threshold_dist_ 0.1f

//==========================================covers cases for
// volume_size_tmp_meters<900 which are not 300,450 or 500
#define POSE_FILE_NAME ""
#define shrunk_volume_x (volume_size_tmp_meters / 2)
#define shrunk_volume_y (volume_size_tmp_meters / 2)
#define shrunk_volume_z (volume_size_tmp_meters / 2)

#if shrink_volume_ == 1
#if getting_correct_pose == 0
#if dataset == ganesh
#if volume_size_tmp == 300
#define shrunk_volume_x 0.3 // 0.12
#define shrunk_volume_y 0.4  // 0.12
#define shrunk_volume_z -0.35 // -0.47 // min is 0.38 below this we do not see full frame
//-0.39 works well but drift from back
//-0.45 covers and tracks better compared to -0.39
#define POSE_FILE_NAME "ganesh_vol300.txt"
#endif
#endif
//===========================
#if dataset == shiva
#if volume_size_tmp == 300
#define shrunk_volume_x 0.3
#define shrunk_volume_y 0.3
#define shrunk_volume_z -0.4
//#define POSE_FILE_NAME "noise0_300vol.txt"//obtained with inital volume of 450
//#define POSE_FILE_NAME "noise0_300_gt.txt" // obtained with initial volume of
// 900 in noise 0 setting
//#define POSE_FILE_NAME "ch_generated_poses.txt" // generated by chuong

#define POSE_FILE_NAME "noise1_300_gt.txt" // obtained with initla volume of 600
// fixed pose obtained with noise1 setting
//#endif

#elif volume_size_tmp == 450
#define shrunk_volume_x 0.29
#define shrunk_volume_y 0.24
#define shrunk_volume_z -0.6
#define POSE_FILE_NAME "noise0_450vol.txt"
//#endif

#elif volume_size_tmp == 500
#define shrunk_volume_x 0.3
#define shrunk_volume_y 0.3
#define shrunk_volume_z -0.6
#define POSE_FILE_NAME "noise0_500vol.txt"
//#endif
#endif
#endif
#if dataset == dino
#if volume_size_tmp == 300
#define shrunk_volume_x 0.15 // 0.12
#define shrunk_volume_y 0.22  // 0.12
#define shrunk_volume_z -0.4 // -0.47 // min is 0.38 below this we do not see full frame
//-0.39 works well but drift from back
//-0.45 covers and tracks better compared to -0.39
#define POSE_FILE_NAME "ganesh_vol300.txt"
#endif
#endif
//#endif
//================overwrite defines if getting correct pose
#else
#if volume_size_tmp >= 300
#if dataset == ganesh
//#define shrunk_volume_x 0.12
//#define shrunk_volume_y 0.12
//#define shrunk_volume_z -0.46 // min is 0.38 below this we do not see full
// frame

#define shrunk_volume_x 0.12
#define shrunk_volume_y 0.12
#define shrunk_volume_z -0.42
//-0.39 works well but drift from back
//-0.42 covers and tracks better compared to -0.39

//#endif
#elif dataset == shiva
#define shrunk_volume_x 0.25
#define shrunk_volume_y 0.23
#define shrunk_volume_z -0.7
#endif
#endif
#endif
#endif
//#========================
#if getting_correct_pose == 0
#if volume_size_tmp >= 900
#define shrink_volume_ 0
#define POSE_FILE_NAME "pose_icp_noise0changed.txt"
#endif
#endif

//=======================================================================

const float Pi = 3.14159265f;
const float PI = 3.14159265f;
#define sigma_color_ 30
#define sigma_space_ 4.5f
#define theta_mean_ Pi / 6.0f
//====================================
#define sig_z_default_value_ 0.1f
#define noise_scaling_factor 1.0
#define reserve_cloud_points 1000000000 // used in tsdf_volume.cpp //3000000 //

#endif /* VOLUME_RELATED_H */