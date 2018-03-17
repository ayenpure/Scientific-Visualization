cmake_minimum_required(VERSION 3.10)

PROJECT(raycaster)
SET(VTK_DIR /Users/abhishek/repositories/VTK-8.1.0_install/lib/cmake/vtk-8.1)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(raycaster raycaster.cxx auxiliary.hxx)

if(VTK_LIBRARIES)
  target_link_libraries(raycaster ${VTK_LIBRARIES})
else()
  target_link_libraries(raycaster vtkHybrid)
endif()


