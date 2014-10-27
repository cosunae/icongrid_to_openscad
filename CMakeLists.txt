cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

PROJECT(ig2vtk)

find_package(VTK REQUIRED)
find_package(NetCDF REQUIRED)
include(${VTK_USE_FILE})

add_executable(ig2vtk ig2vtk.cpp)

if(VTK_LIBRARIES)
  target_link_libraries(ig2vtk ${VTK_LIBRARIES} netcdf)
else()
  target_link_libraries(ig2vtk vtkHybrid vtkWidgets netcdf)
endif()

