cmake_minimum_required (VERSION 2.6)

project (p3pf)

# modify this according to your PC
set( EIGEN_INCLUDE_DIR "D:/eigen" )
find_package (glog 0.3.4 REQUIRED)

include_directories (
${CMAKE_SOURCE_DIR} 
${CMAKE_SOURCE_DIR}/p3pf
${CMAKE_SOURCE_DIR}/p6p
${EIGEN_INCLUDE_DIR}
)

# add our macro directory to cmake search path
if (EXISTS "${CMAKE_SOURCE_DIR}/cmake")
  set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
endif ()

include("${CMAKE_MODULE_PATH}/OptimizeCompilerFlags.cmake")
OptimizeCompilerFlags()

set(OpenCV_DIR "D:/Program Files/OpenCV/build")
find_package(OpenCV REQUIRED)
include_directories(
${OpenCV_CONFIG_PATH}
${OpenCV_CONFIG_PATH}/include
${OpenCV_CONFIG_PATH}/include/opencv)

if(OpenCV_FOUND)
message(STATUS "find OpenCV")
message(STATUS ${OpenCV_CONFIG_PATH})
else()
message(STATUS "not find OpenCV")
endif()


#add_subdirectory (p6p)
add_subdirectory (p3pf)
add_subdirectory (epnp)

# add pnp library
add_library(pnp pnp.cc camera_pose.cc p6p/p6p.cpp)
target_link_libraries(pnp epnp glog::glog p3pf external_pose_algs)

# build example for epnp
add_executable(epnp_example epnp_example.cc)
target_link_libraries(epnp_example pnp)

# build example for p3pf
add_executable(p3pf_example p3pf_example.cc)
target_link_libraries(p3pf_example pnp)

# build example for p6p
add_executable(p6p_example p6p_example.cc p6p/p6p.cpp)
target_link_libraries(p6p_example pnp)