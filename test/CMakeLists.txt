cmake_minimum_required( VERSION 3.22.1 )

project( test VERSION 1.0 LANGUAGES CXX )

add_executable( test test.cpp )

set( CMAKE_CXX_STANDARD 20 )
target_compile_options( test PUBLIC -fsanitize=leak -g )
# link_directories( test "../../build/")
find_package( SDL2pp REQUIRED)
find_library( EFD 
              NAMES eulerian_fluid_dynamics
              PATHS ~/codes/eulerian_fluid_simulator/build/ )

target_link_libraries( test VT ${EFD} )
