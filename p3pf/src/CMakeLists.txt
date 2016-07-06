cmake_minimum_required (VERSION 2.6)

if (EXISTS "${CMAKE_SOURCE_DIR}/cmake")
  set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
endif ()

find_package (glog 0.3.4 REQUIRED)

include_directories (
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/p3pf
  ${EIGEN_INCLUDE_DIR}
)

if(EIGEN_INCLUDE_DIR)
  message("in ./src/CMakeLists.txt ${EIGEN_INCLUDE_DIR}")
endif()

if(GLOG_INCLUDE_DIR)
  message("in ./src/CMakeLists.txt ${GLOG_INCLUDE_DIR}")
endif()


set(P3PF_SRC ${CMAKE_SOURCE_DIR}/camera_pose.cc p3p_solver.cc p3pf.cc)

#FILE(GLOB P3PF_HDRS *.h)
add_library(p3pf ${P3PF_SRC})
target_link_libraries(p3pf glog::glog external_pose_algs)

#add_executable(p3pf_example p3pf_example.cc)
#target_link_libraries(p3pf_example p3pf glog::glog external_pose_algs)

#install(TARGETS p3pf LIBRARY RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
#install(FILES ${P3PF_HDRS} RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
#install(PROGRAMS ${CMAKE_BINARY_DIR}/src/p3pf_example   RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
