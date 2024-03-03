set(SUBSUBSYS_NAME tools)
set(SUBSUBSYS_DESC "Kinfu tools")
set(SUBSUBSYS_DEPS gpu_kinfu_zivid_fahira visualization)
set(SUBSUBSYS_OPT_DEPS opencv)
set(EXT_DEPS glew openni)
set(DEFAULT TRUE)
set(REASON "")

PCL_SUBSUBSYS_OPTION(build ${SUBSYS_NAME} ${SUBSUBSYS_NAME} ${SUBSUBSYS_DESC} ${DEFAULT} ${REASON})
PCL_SUBSYS_DEPEND(build NAME ${SUBSUBSYS_NAME} PARENT_NAME ${SUBSYS_NAME} DEPS ${SUBSUBSYS_DEPS} OPT_DEPS ${SUBSUBSYS_OPT_DEPS} EXT_DEPS ${EXT_DEPS})

if(NOT build)
  return()
endif()

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include" "/usr/include/x86_64-linux-gnu/c++/11")

file(GLOB hdrs "*.h*")
include_directories(SYSTEM ${OPENNI_INCLUDE_DIRS})

## KINECT FUSION
set(the_target pcl_kinfu_app_zivid_fahira)
set(srcs kinfu_app_zivid.cpp capture_zivid.cpp evaluation_zivid.cpp)

source_group("Source Files" FILES ${srcs})

PCL_ADD_EXECUTABLE(${the_target} COMPONENT ${SUBSUBSYS_NAME} SOURCES ${srcs} ${hdrs} BUNDLE)
target_link_libraries(${the_target} pcl_common pcl_io pcl_surface  ${OPENNI_LIBRARIES} ${PCL_LIBRARIES} pcl_visualization pcl_gpu_kinfu_zivid_fahira)

if(OpenCV_FOUND)
  include_directories(${OpenCV_INCLUDE_DIRS})
  target_link_libraries("${the_target}" ${OpenCV_LIBS})
endif()

## mesh the cloud
# set(the_target pcl_kinfu_app_zivid_meshit)
# set(srcs mesh_it.cpp )

# source_group("Source Files" FILES ${srcs})

# PCL_ADD_EXECUTABLE(${the_target} COMPONENT ${SUBSUBSYS_NAME} SOURCES ${srcs} ${hdrs} BUNDLE)
# target_link_libraries(${the_target} pcl_common pcl_io pcl_surface  Boost::program_options pcl_features ${PCL_LIBRARIES})


##bilinear interpolate_check

# set(the_target pcl_kinfu_app_zivid_bilinear_interpolate)
# set(srcs bilinear_interpolate_main.cpp )

# source_group("Source Files" FILES ${srcs})

# PCL_ADD_EXECUTABLE(${the_target} COMPONENT ${SUBSUBSYS_NAME} SOURCES ${srcs} ${hdrs} BUNDLE)
# target_link_libraries(${the_target} Boost::program_options )



# ## RECORD TSDFVOLUME
# set(the_target pcl_record_tsdfvolume_zivid)
# set(srcs record_tsdfvolume_zivid.cpp capture_zivid.cpp)

# PCL_ADD_EXECUTABLE(${the_target} COMPONENT ${SUBSUBSYS_NAME} SOURCES ${srcs} ${hdrs})
# target_link_libraries(${the_target} pcl_common pcl_io pcl_surface ${PCL_LIBRARIES} ${OPENNI_LIBRARIES} pcl_visualization pcl_gpu_kinfu_zivid_fahira)
