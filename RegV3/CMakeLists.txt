cmake_minimum_required(VERSION 2.8)

file(GLOB_RECURSE srcs ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB_RECURSE incs ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hxx ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(ITK_REG_PROJ_NAME RegistrationV3)

set(VIS_MODE "WO_VISUALIZATION" CACHE STRING "Help message for this variable")
set_property(
  CACHE VIS_MODE
  PROPERTY STRINGS
  "WO_VISUALIZATION" "W_VISUALIZATION"
)

project(${ITK_REG_PROJ_NAME})
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
else()
  find_package(ItkVtkGlue REQUIRED)
  include(${ItkVtkGlue_USE_FILE})
  set(Glue ItkVtkGlue)
endif()
 
add_executable(${ITK_REG_PROJ_NAME} ${srcs} ${incs})
target_link_libraries(${ITK_REG_PROJ_NAME} ${Glue} ${VTK_LIBRARIES} ${ITK_LIBRARIES})

 # Properties-&gt;C/C++-&gt;General-&gt;Additional Include Directories
include_directories(${INCLUDE_DIR})
add_definitions(-D${VIS_MODE})