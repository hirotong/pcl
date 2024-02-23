#ifndef ZIVID_RELATED_H
#define ZIVID_RELATED_H

#include "camera_type.h"
#include <iostream>

#if current_camera_type == ZIVID_CAMERA
#define hh_1 1200 // height
#define ww_1 1944 // width

#define depth_focal_x 1779.87915039062f
#define depth_focal_y 1780.30529785156f

const float f_x1 = 1779.87915039062;
const float f_y1 = 1780.30529785156;
const float c_x1 = 959.119079589844;
const float c_y1 = 573.3330078125;

const float Height1 = hh_1; // 150 downsampled
const float Width1 = ww_1;  // 243
#endif
//==========================to run kinect
#if current_camera_type == KINECT_CAMERA
#define hh_1 480 // height
#define ww_1 640 // width

#define depth_focal_x 585.f
#define depth_focal_y 585.f

const float f_x1 = 525.0f;
const float f_y1 = 525.0f;
const float c_x1 = 319.5f;
const float c_y1 = 239.5f;

const float Height1 = hh_1;
const float Width1 = ww_1;
#endif
//===========================
/*
const float Height1 = 480;//600;//300;//150;//1200;//150 downsampled
const float Width1 = 640;//972;//486;//243;//1944;//243

constexpr float Height2 = 480;//600;//;300;//150;//1200;
constexpr float Width2 = 640;//972;//486;//243;//1944;
*/

//#define ff_xx 1779.87915039062

#endif /* ZIVID_RELATED_H */
