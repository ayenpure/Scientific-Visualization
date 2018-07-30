cmake_minimum_required(VERSION 3.10)

PROJECT(proj8)
SET(VTK_DIR /Users/abhishek/repositories/VTK-8.1.0_install/lib/cmake/vtk-8.1)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj8 proj8.cxx)

SET(CMAKE_CXX_LINK_FLAGS "-g -framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(proj8 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj8 vtkHybrid)
endif()


