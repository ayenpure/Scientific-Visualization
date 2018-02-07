cmake_minimum_required(VERSION 3.9)
set (CMAKE_CXX_STANDARD 11)

project(analytical)

find_package(VTKm REQUIRED QUIET
             OPTIONAL_COMPONENTS Serial CUDA
            )

add_executable(analytical analyticaldataset.cxx)
target_include_directories(analytical PRIVATE ${VTKm_INCLUDE_DIRS})
target_link_libraries(analytical ${VTKm_LIBRARIES})
target_compile_options(analytical PRIVATE ${VTKm_COMPILE_OPTIONS})

add_executable(points points.cxx)
target_compile_options(points PRIVATE ${VTKm_COMPILE_OPTIONS})

