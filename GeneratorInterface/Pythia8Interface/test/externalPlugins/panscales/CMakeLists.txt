cmake_minimum_required(VERSION 3.15)
project(CMSPanScalesPlugin LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
find_package(PanScales REQUIRED CONFIG)
find_path(PYTHIA_INCLUDE Pythia8/Pythia.h HINTS "${PYTHIA_ROOT}/include" REQUIRED)
find_library(PYTHIA_LIBRARY pythia8 HINTS "${PYTHIA_ROOT}/lib" REQUIRED)
add_library(CMSPanScales SHARED PanScalesPlugin.cc
  "${PANSCALES_SOURCE}/pythia-interface/PanScalesPythiaModule.cc"
  "${PANSCALES_SOURCE}/pythia-interface/PanScalesPythiaProcess.cc")
target_include_directories(CMSPanScales PRIVATE "${PANSCALES_SOURCE}/pythia-interface" "${PYTHIA_INCLUDE}")
target_link_libraries(CMSPanScales PRIVATE PanScales::panscales "${PYTHIA_LIBRARY}")
install(TARGETS CMSPanScales LIBRARY DESTINATION lib)
