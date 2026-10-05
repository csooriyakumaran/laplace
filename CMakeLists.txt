cmake_minimum_required(VERSION 3.25)

project(laplace
  LANGUAGES C CXX Fortran
)

# Single-config generators (Makefiles, Ninja) leave this empty by default,
# which applies none of the per-config flags below (e.g. -O3) -- silently
# building unoptimized. Multi-config generators (Visual Studio, Ninja Multi-
# Config) use CMAKE_CONFIGURATION_TYPES instead, so leave those alone.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- KERNEL LIBRARY (dependency-free: stdlib + ISO_C_BINDING only) ---
add_library(laplacelib STATIC
    src/laplace_types.f90
    src/laplace_core.f90
    src/laplace_f_api.f90
    src/laplace_c_api.f90
    src/laplace_atomic.c
)

set_target_properties(laplacelib PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod
)

target_include_directories(laplacelib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_BINARY_DIR}/mod
)

# --- FORTRAN APP (native reference driver, no deps) ---
add_executable(fortran_app
    application/fortran/main.f90
)

target_link_libraries(fortran_app
    PRIVATE laplacelib
)

# --- CPP APP (demo driver, uses aether for arena-managed buffers) ---
include(FetchContent)

FetchContent_Declare(
    aether
    GIT_REPOSITORY https://github.com/csooriyakumaran/aether.git
    GIT_TAG v0.0.18
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(aether)

add_executable(cpp_app
    application/cpp/main.cpp
)

target_link_libraries(cpp_app
    PRIVATE laplacelib aether::aether
)

add_executable(c_app
    application/c/main.c
)

target_link_libraries(c_app
    PRIVATE laplacelib aether::aether
)

# --- TESTS ---
enable_testing()

add_executable(test_straight_duct
    tests/c/test_straight_duct.c
)

target_link_libraries(test_straight_duct
    PRIVATE laplacelib aether::aether
)

add_test(NAME straight_duct COMMAND test_straight_duct)

add_executable(test_straight_duct_f
    tests/fortran/test_straight_duct.f90
)

target_link_libraries(test_straight_duct_f
    PRIVATE laplacelib
)

add_test(NAME straight_duct_fortran COMMAND test_straight_duct_f)

add_executable(test_thomas_solve
    tests/fortran/test_thomas_solve.f90
)

target_link_libraries(test_thomas_solve
    PRIVATE laplacelib
)

add_test(NAME thomas_solve COMMAND test_thomas_solve)

add_executable(test_sweep
    tests/fortran/test_sweep.f90
)

target_link_libraries(test_sweep
    PRIVATE laplacelib
)

add_test(NAME sweep COMMAND test_sweep)

add_executable(test_solve
    tests/fortran/test_solve.f90
)

target_link_libraries(test_solve
    PRIVATE laplacelib
)

add_test(NAME solve COMMAND test_solve)

add_executable(test_output
    tests/fortran/test_output.f90
)

target_link_libraries(test_output
    PRIVATE laplacelib
)

add_test(NAME output COMMAND test_output)

add_executable(test_atomic
    tests/fortran/test_atomic.f90
    tests/fortran/test_atomic_probe.c
)

target_link_libraries(test_atomic
    PRIVATE laplacelib
)

add_test(NAME atomic COMMAND test_atomic)
