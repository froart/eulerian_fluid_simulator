cmake_minimum_required(VERSION 3.22.1) 
project(eulerian_fluid_simulator)
add_library(eulerian_fluid_simulator src/opengl_setup.cpp src/fluid_dynamics.cpp)
install(TARGETS eulerian_fluid_simulator DESTINATION lib)
