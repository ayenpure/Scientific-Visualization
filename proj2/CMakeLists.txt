cmake_minimum_required(VERSION 3.10)

PROJECT(proj2)
SET(VTK_DIR /Users/abhishek/repositories/VTK-8.1.0_install/lib/cmake/vtk-8.1)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj2 MACOSX_BUNDLE proj2)

SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(proj2 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj2 vtkHybrid)
endif()


