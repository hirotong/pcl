#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>

#include <boost/program_options.hpp>

#include <pcl/common/time.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <pcl/features/normal_3d.h>
#include <pcl/surface/poisson.h>

using namespace pcl;
using namespace Eigen;

namespace pc = pcl::console;

using Point_t=pcl::PointXYZ;
//using Point_t=pcl::PointNormal;

using Cloud_t=pcl::PointCloud<Point_t>;

void pcl_reconstruct( Cloud_t::Ptr  cloud, std::string name)
{
  
std::cout<<"computing normals"<<std::endl;
  search::KdTree<Point_t>::Ptr tree(new search::KdTree<Point_t>);
     PointCloud<Normal>::Ptr cloud_normals(new PointCloud<Normal>);
    NormalEstimation<Point_t, Normal> ne;
    ne.setInputCloud(cloud);
    ne.setSearchMethod(tree);
    //ne.setRadiusSearch(0.5);
    ne.setRadiusSearch (0.03);
    ne.compute(*cloud_normals);
  
  PointCloud<PointNormal>::Ptr cloud_with_normals(new PointCloud<PointNormal>);
    concatenateFields(*cloud, *cloud_normals, *cloud_with_normals);
    
    //Poisson<Point_t> recons;
    Poisson<PointNormal> recons;
    std::cout<<"inputting cloud to reconst"<<std::endl;
    recons.setDepth(10);
    recons.setSamplesPerNode(4);
    recons.setPointWeight(2);
    recons.setIsoDivide(10);
    recons.setInputCloud(cloud_with_normals);
   // recons.setInputCloud(cloud);
 std::cout<<"meshing"<<std::endl;
    pcl::PolygonMesh mesh;
    recons.reconstruct(mesh);

   io::savePLYFile(name, mesh);
   std::cout<<"saved mesh"<<name<<std::endl;
}

int main (int argc, char * argv[])
{

     
namespace po = boost::program_options;

     po::options_description desc("Options");
    desc.add_options()
        ("help", "Show the help message")
        ("pcd", po::value<std::string>(), " point cloud")
        ("ply_saved_name", po::value<std::string>(), "the name with full path hwere you want to save the output mesh")
        
    ;

    // Configure positional arguments
    po::positional_options_description pos_args;
    pos_args.add("pcd", 1);
    pos_args.add("ply_saved_name", 1);
    po::variables_map vm;
    po::store(
            po::command_line_parser(argc, argv)
                .options(desc)
                .positional(pos_args)
                .run(),
            vm
    );
    po::notify(vm);

    if(vm.count("help") || !vm.count("pcd") || !vm.count("ply_saved_name")
               )
    {
        std::cout << desc << std::endl;
        return 1;
    }
 Cloud_t::Ptr cloud(new Cloud_t);
if (pcl::io::loadPCDFile<Point_t>(vm["pcd"].as<std::string>(), *cloud)==-1)
    {
        std::cout << "Could not read source file" << std::endl;
        return -1;
    }

pcl_reconstruct( cloud,vm["ply_saved_name"].as<std::string>());
}