cmake_minimum_required(VERSION 3.22.1) 
project(fluid_dynamics)
add_library(${PROJECT_NAME} fluid_dynamics.cpp)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER fluid_dynamics.hpp)
install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION lib
				PUBLIC_HEADER DESTINATION include	
				)
