DATAROOT="/home/ton068/datasets/Noise_model/kinect_data"
for i in `ls $DATAROOT`; do
if [[ -d $DATAROOT/$i ]]; then
  echo "./build/bin/pcl_kinfu_app_noise_modeling_debug -r -pcd $DATAROOT/$i/pcd";
  mv $DATAROOT/$i/kinfu_gtpose $DATAROOT/$i/kinfu_gtpose_ori;
./build/bin/pcl_kinfu_app_noise_modeling_debug -gt_pose $DATAROOT/$i/associate_gt.txt -r -pcd $DATAROOT/$i/pcd -z0 0.5 -volume_size 3.0 -nc 2 -ts 1000;
./build/bin/pcl_kinfu_app_noise_modeling_debug -gt_pose $DATAROOT/$i/associate_gt.txt -r -pcd $DATAROOT/$i/pcd -z0 0.5 -volume_size 3.0 -nc 0 -ts 1000;
fi
done

