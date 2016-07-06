cmake_minimum_required (VERSION 3.0)


include_directories (
${CMAKE_SOURCE_DIR} 
${CMAKE_SOURCE_DIR}/epnp
${OpenCV_CONFIG_PATH}
${OpenCV_CONFIG_PATH}/include
${OpenCV_CONFIG_PATH}/include/opencv)


# add pnp library
add_library(epnp epnp.cpp epnp.h)
target_link_libraries(epnp ${OpenCV_LIBS})

