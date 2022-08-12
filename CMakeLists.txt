cmake_minimum_required(VERSION 3.22)
project(RADMC)

set(CMAKE_CXX_STANDARD 14)

add_executable(RADMC src/radmc.cc include/grid.hh src/grid.cc include/frequencies.hh src/frequencies.cc src/stars.cc include/stars.hh src/star.cc include/star.hh src/dust.cc include/dust.hh src/dust_species.cc include/dust_species.hh)
