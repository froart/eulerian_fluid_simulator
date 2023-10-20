cmake_minimum_required( VERSION 3.22.1 )

project( eulerian_fluid_dynamics )

file( GLOB LIBRARY_SOURCES "./src/eulerian_fluid_dynamics.cpp")
file( GLOB LIBRARY_HEADERS "./include/eulerian_fluid_dynamics.hpp")

add_library( eulerian_fluid_dynamics SHARED ${LIBRARY_SOURCES} )

set( CMAKE_CXX_STANDARD 20 )
target_compile_options( eulerian_fluid_dynamics PUBLIC -fsanitize=leak -g )
find_package( SDL2pp REQUIRED )
target_link_libraries( eulerian_fluid_dynamics PUBLIC VT SDL2pp::SDL2pp )

# set_target_properties( eulerian_fluid_dynamics PROPERTIES PUBLIC_HEADER eulerian_fluid_dynamics )
