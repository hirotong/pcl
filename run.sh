DATAROOT="/home/ton068/datasets/kinect_data"
for i in `ls $DATAROOT`; do
if [[ -d $DATAROOT/$i ]]; then
  echo "./build/bin/pcl_kinfu_app_noise_modeling -r -pcd $DATAROOT/$i/pcd";
./build/bin/pcl_kinfu_app_noise_modeling_debug -r -pcd $DATAROOT/$i/pcd -z0 0.5 -nc 0;
./build/bin/pcl_kinfu_app_noise_modeling_debug -r -pcd $DATAROOT/$i/pcd -z0 0.5 -nc 1;
./build/bin/pcl_kinfu_app_noise_modeling_debug -r -pcd $DATAROOT/$i/pcd -z0 0.5 -nc 2
fi
done

