#ifndef CAMERA_TYPE_H
#define CAMERA_TYPE_H

#define ZIVID_CAMERA 1
#define KINECT_CAMERA 2

#define shiva 3
#define ganesh 4
#define dino 5

#define current_camera_type ZIVID_CAMERA
#define dataset dino // dino // shiva // ganesh

// set shrink volume to 1 if using volume less than 900mm for zivid
// set getting_correct_pose to 1 if want to get better pose for smaller volume
#define shrink_volume_ 1
#define getting_correct_pose 0

#define read_saved_pose_ 1
#define write_pose_to_file_ 1
#define write_pose_file_name "noise2_tracked_noise1_rr.txt"

#define reconst_every_n_frame 0
#define n_th_frame 8

#define block_scale_x 32
#define block_scale_y 8

#define set_icp 1 // sets icp to zero or 1.

#endif /* CAMERA_TYPE_H */
