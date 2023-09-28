cmake_minimum_required( VERSION 3.22.1 )

project( eulerian_fluid_dynamics )

set( CMAKE_CXX_STANDARD 20 )

add_library( eulerian_fluid_dynamics SHARED eulerian_fluid_dynamics.cpp)

target_link_libraries( eulerian_fluid_dynamics PUBLIC screen )

target_compile_options( eulerian_fluid_dynamics PUBLIC -fsanitize=leak -g )

set_target_properties( eulerian_fluid_dynamics PROPERTIES PUBLIC_HEADER eulerian_fluid_dynamics )

install( TARGETS eulerian_fluid_dynamics
         LIBRARY DESTINATION lib 
         PUBLIC_HEADER DESTINATION inlcude )
