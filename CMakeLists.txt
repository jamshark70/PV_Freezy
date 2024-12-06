cmake_minimum_required (VERSION 2.8)
project (PV_Freezish)

include_directories(${SC_PATH}/include/plugin_interface)
include_directories(${SC_PATH}/include/common)
include_directories(${SC_PATH}/external_libraries/libsndfile/)

set(CMAKE_CXX_STANDARD 17)


set(CMAKE_SHARED_MODULE_PREFIX "")
if(APPLE OR WIN32)
set(CMAKE_SHARED_MODULE_SUFFIX ".scx")
endif()

add_library(PV_Freezish MODULE source/PVFreezy.cpp)

if (APPLE)
	target_link_libraries(PV_Freezish "-framework Accelerate")
endif()
