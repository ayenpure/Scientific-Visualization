cmake_minimum_required(VERSION 3.10)

PROJECT(proj3)
SET(VTK_DIR /Users/abhishek/repositories/VTK-8.1.0_install/lib/cmake/vtk-8.1)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj3 MACOSX_BUNDLE proj3)

SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(proj3 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj3 vtkHybrid)
endif()


