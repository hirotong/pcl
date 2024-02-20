set(SUBSUBSYS_NAME tools)
set(SUBSUBSYS_DESC "Kinfu tools noise modeling")
set(SUBSUBSYS_DEPS gpu_kinfu_noise_modeling visualization)
set(SUBSUBSYS_OPT_DEPS opencv)
set(EXT_DEPS glew openni)
set(DEFAULT TRUE)
set(REASON "")

PCL_SUBSUBSYS_OPTION(build ${SUBSYS_NAME} ${SUBSUBSYS_NAME} ${SUBSUBSYS_DESC} ${DEFAULT} ${REASON})
PCL_SUBSYS_DEPEND(build NAME ${SUBSUBSYS_NAME} PARENT_NAME ${SUBSYS_NAME} DEPS ${SUBSUBSYS_DEPS} OPT_DEPS ${SUBSUBSYS_OPT_DEPS} EXT_DEPS ${EXT_DEPS})

if(NOT build)
  return()
endif()

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

file(GLOB hdrs "*.h*")
include_directories(SYSTEM ${OPENNI_INCLUDE_DIRS} /usr/include/x86_64-linux-gnu/c++/9)

## KINECT FUSION
set(the_target pcl_kinfu_app_noise_modeling)
set(srcs kinfu_app.cpp capture.cpp evaluation.cpp)

source_group("Source Files" FILES ${srcs})

PCL_ADD_EXECUTABLE(${the_target} COMPONENT ${SUBSUBSYS_NAME} SOURCES ${srcs} ${hdrs} BUNDLE)
target_link_libraries(${the_target} pcl_common pcl_io ${OPENNI_LIBRARIES} pcl_visualization pcl_gpu_kinfu_noise_modeling)

find_package(OpenCV REQUIRED)
if(OpenCV_FOUND)
  add_compile_definitions(HAVE_OPENCV)
  include_directories(${OpenCV_INCLUDE_DIRS})
  target_link_libraries("${the_target}" ${OpenCV_LIBS})
endif()

## RECORD TSDFVOLUME
# set(the_target pcl_record_tsdfvolume)
# set(srcs record_tsdfvolume.cpp capture.cpp)

# PCL_ADD_EXECUTABLE(${the_target} COMPONENT ${SUBSUBSYS_NAME} SOURCES ${srcs} ${hdrs})
# target_link_libraries(${the_target} pcl_common pcl_io ${OPENNI_LIBRARIES} pcl_visualization pcl_gpu_kinfu)
