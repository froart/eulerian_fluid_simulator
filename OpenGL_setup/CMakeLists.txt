cmake_minimum_required(VERSION 3.22.1) 
project(opengl_setup)
add_library(${PROJECT_NAME} opengl_setup.cpp)
target_link_libraries(${PROJECT_NAME} GL GLU glut)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER opengl_setup.hpp)
install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION lib
				PUBLIC_HEADER DESTINATION include	
				)
