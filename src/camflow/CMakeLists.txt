# Specification file for the CMake build system
#
# Robert I A Patterson
# 30 May 2011
cmake_minimum_required(VERSION 2.6)

project (CAMFLOW)

add_library(camflow source/cam_read.cpp
                    source/cam_setup.cpp
                    source/cam_geometry.cpp
                    source/cam_control.cpp
                    source/cam_conveter.cpp
                    source/cam_math.cpp
                    source/cam_residual.cpp
                    source/cam_conc.cpp
                    source/cam_admin.cpp
                    source/cam_boundary.cpp
                    source/cam_configuration.cpp
                    source/cam_profile.cpp
                    source/cam_reporter.cpp
                    source/cam_radiation.cpp
                    source/cam_models.cpp
                    source/cam_soot.cpp
                    source/cvode_wrapper.cpp
                    source/radau_wrapper.cpp
                    source/cam_sdr.cpp
                    source/flamelet.cpp
                    source/batch.cpp)

include_directories(include 
                    ../utils/include
                    ../io/chemkinReader/include
                    ../io/comostrings/include
                    ../io/camxml/include
                    ../sprogc/include
                    ../odesolvers/cvode/include
                    ../odesolvers/radau/include)

# To link againast camflow the following libraries will also
# be required to resolve references from within camflow.
target_link_libraries(camflow camxml comostrings sprog cvode radau chemkinReader)
