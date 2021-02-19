cmake_minimum_required(VERSION 3.16)

project(GPR5204)

set (CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    set(COMPILER_FLAGS " /arch:AVX2 /Oi /GR- /EHs-c- /fp:fast")
elseif(UNIX)
    set(COMPILER_FLAGS " -fno-rtti -fno-exceptions -ffast-math -save-temps -fno-omit-frame-pointer -masm=intel")
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL x86_64 OR CMAKE_SYSTEM_PROCESSOR STREQUAL amd64)
        set(COMPILER_FLAGS "${COMPILER_FLAGS} -march=haswell")
    endif()
endif()
message("${CMAKE_HOST_SYSTEM_PROCESSOR}")

file(GLOB_RECURSE SRC_FILES src/*.cpp include/*.h)
add_library(Common STATIC ${SRC_FILES})
target_include_directories(Common PUBLIC "include/")
set_property(TARGET Common APPEND_STRING PROPERTY COMPILE_FLAGS " ${COMPILER_FLAGS} ")

find_package(benchmark CONFIG REQUIRED)

file(GLOB_RECURSE BENCH_FILES bench/*.cpp)
foreach(bench_file ${BENCH_FILES})
    get_filename_component(bench_name ${bench_file} NAME_WE)
    add_executable(${bench_name} ${bench_file})
    target_link_libraries(${bench_name} PRIVATE Common)
    target_link_libraries(${bench_name} PRIVATE benchmark::benchmark benchmark::benchmark_main)
    set_property(TARGET ${bench_name} APPEND_STRING PROPERTY COMPILE_FLAGS " ${COMPILER_FLAGS} ")
endforeach()
