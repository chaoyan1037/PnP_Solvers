cmake_minimum_required (VERSION 2.6)

if (EXISTS "${CMAKE_SOURCE_DIR}/cmake")
  set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
endif ()

if(EIGEN_INCLUDE_DIR)
  message(STATUS "in ./ext/CMakeLists.txt eigen_dir: ${EIGEN_INCLUDE_DIR}")
endif()

find_package (glog 0.3.4 REQUIRED)
if(glog_FOUND)
message (STATUS "in ./ext/CMakeLists.txt glog found")
endif()

include_directories (
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/ext
  ${EIGEN_INCLUDE_DIR}
)


# This should be sufficient to compile the external algorithms. The header files should be found automatically.
set(EXTERNAL_SRC ${CMAKE_SOURCE_DIR}/p3pf/ext/p3p_code_kneip/P3p.cpp)
add_library(external_pose_algs ${EXTERNAL_SRC})
target_link_libraries(external_pose_algs glog::glog)
