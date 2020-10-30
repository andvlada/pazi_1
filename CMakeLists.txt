cmake_minimum_required(VERSION 3.5)

project(pazi_1 LANGUAGES C)

find_package(PkgConfig REQUIRED)
pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)

add_executable(pazi_1 main.c jacobi.h jacobi.c)

target_link_libraries(pazi_1 PRIVATE PkgConfig::gmp)
