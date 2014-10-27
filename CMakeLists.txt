cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

PROJECT(ig2vtk)

find_package(VTK REQUIRED)
find_package(NetCDF REQUIRED)
find_package(Boost COMPONENTS program_options filesystem system)
include(${VTK_USE_FILE})

add_executable(ig2vtk ig2vtk.cpp)

set(IG2VTK_LIBRARIES netcdf ${Boost_LIBRARIES})

if(VTK_LIBRARIES)
  target_link_libraries(ig2vtk ${VTK_LIBRARIES} ${IG2VTK_LIBRARIES})
else()
  target_link_libraries(ig2vtk vtkHybrid vtkWidgets ${IG2VTK_LIBRARIES})
endif()

