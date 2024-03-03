DATAROOT="/home/ton068/datasets/Noise_model/kinect_data"
for i in $(ls $DATAROOT); do
  if [[ -d $DATAROOT/$i ]]; then
    echo "./build/bin/pcl_kinfu_app_noise_modeling -r -pcd $DATAROOT/$i/pcd"
    mv $DATAROOT/$i/kinfu $DATAROOT/$i/kinfu_ori
    ./build/bin/pcl_kinfu_app_noise_modeling_debug -r -pcd $DATAROOT/$i/pcd -z0 0.5 -nc 0
    ./build/bin/pcl_kinfu_app_noise_modeling_debug -r -pcd $DATAROOT/$i/pcd -z0 0.5 -nc 1
    ./build/bin/pcl_kinfu_app_noise_modeling_debug -r -pcd $DATAROOT/$i/pcd -z0 0.5 -nc 2
    ./build/bin/pcl_kinfu_app_noise_modeling_debug -r -pcd $DATAROOT/$i/pcd -z0 0.5 -nc 0 -volume_size 1
    ./build/bin/pcl_kinfu_app_noise_modeling_debug -r -pcd $DATAROOT/$i/pcd -z0 0.5 -nc 1 -volume_size 1
    ./build/bin/pcl_kinfu_app_noise_modeling_debug -r -pcd $DATAROOT/$i/pcd -z0 0.5 -nc 2 -volume_size 1
  fi
done
