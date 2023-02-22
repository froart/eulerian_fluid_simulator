cmake_minimum_required(VERSION 3.22.1) 
project(eulerian_fluid_simulator)
add_library(eulerian_fluid_simulator opengl_setup.cpp fluid_dynamics.cpp)
