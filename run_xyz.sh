DATAROOT="/home/ton068/datasets/Noise_model/kinect_data/rgbd_dataset_freiburg2_xyz"

./build/bin/pcl_kinfu_app_noise_modeling_debug -ic -r -pcd $DATAROOT/pcd -z0 0.5 -nc 2 -ts 100
./build/bin/pcl_kinfu_app_noise_modeling_debug -ic -r -pcd $DATAROOT/pcd -z0 0.5 -nc 0 -ts 100
./build/bin/pcl_kinfu_app_noise_modeling_debug -ic -r -pcd $DATAROOT/pcd -z0 0.5 -nc 2 -ts 1000
./build/bin/pcl_kinfu_app_noise_modeling_debug -ic -r -pcd $DATAROOT/pcd -z0 0.5 -nc 0 -ts 1000
./build/bin/pcl_kinfu_app_noise_modeling_debug -ic -r -pcd $DATAROOT/pcd -z0 0.5 -nc 2 -ts 10000
./build/bin/pcl_kinfu_app_noise_modeling_debug -ic -r -pcd $DATAROOT/pcd -z0 0.5 -nc 0 -ts 10000

