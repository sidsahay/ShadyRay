cmake_minimum_required(VERSION 3.5)
project(ShadyRay)

set(CMAKE_CXX_STANDARD 14)

#       Raytracer Executable
add_executable(ShadyRay
        ShadyRay/src/defs.hpp
        ShadyRay/src/shadyray.cpp
        ShadyRay/src/camera.cpp
        ShadyRay/src/camera.hpp
        ShadyRay/src/film.cpp
        ShadyRay/src/film.hpp
        ShadyRay/src/filter.cpp
        ShadyRay/src/filter.hpp
        ShadyRay/src/integrator.hpp
        ShadyRay/src/integratorth.hpp
        ShadyRay/src/interaction.hpp
        ShadyRay/src/light.cpp
        ShadyRay/src/light.hpp
        ShadyRay/src/linalg.cpp
        ShadyRay/src/linalg.hpp
        ShadyRay/src/material.cpp
        ShadyRay/src/material.hpp
        ShadyRay/src/pathintegrator.hpp
        ShadyRay/src/pathintegratorth.hpp
        ShadyRay/src/sampler.cpp
        ShadyRay/src/sampler.hpp
        ShadyRay/src/scene.cpp
        ShadyRay/src/scene.hpp
        ShadyRay/src/sdlfilm.cpp
        ShadyRay/src/sdlfilm.hpp
        ShadyRay/src/spectrum.cpp
        ShadyRay/src/spectrum.hpp
        ShadyRay/src/geometry.cpp
        ShadyRay/src/geometry.hpp
        ShadyRay/src/memorypool.cpp
        ShadyRay/src/memorypool.hpp
        ShadyRay/src/texture.hpp
        ShadyRay/src/image.cpp
        ShadyRay/src/image.hpp
        ShadyRay/src/ppm.cpp
        ShadyRay/src/ppm.hpp
        ShadyRay/src/texturemap.cpp
        ShadyRay/src/texturemap.hpp
        ShadyRay/src/texture.cpp
        ShadyRay/src/ray.cpp
        ShadyRay/src/ray.hpp
        ShadyRay/src/luaintegration.hpp
        ShadyRay/src/particletracer.cpp
        ShadyRay/src/particletracer.hpp
        ShadyRay/dkm/dkm.hpp
        ShadyRay/dkm/dkm_utils.hpp)

find_library(EMBREE_LIB embree)
if(NOT EMBREE_LIB)
    message(WARNING "Intel Embree 2.x library not found explicitly. Build may still work if the linker can find it.")
    target_link_libraries(ShadyRay embree)
else()
    target_link_libraries(ShadyRay ${EMBREE_LIB})
endif()

find_library(ASSIMP_LIB assimp)
if(NOT ASSIMP_LIB)
    message(WARNING "Assimp library not found explicitly. Build may still work if the linker can find it.")
    target_link_libraries(ShadyRay assimp)
else()
    target_link_libraries(ShadyRay ${ASSIMP_LIB})
endif()

find_library(SDL2_LIB SDL2)
if(NOT SDL2_LIB)
    message(WARNING "SDL2 library not found explicitly. Build may still work if the linker can find it.")
    target_link_libraries(ShadyRay SDL2)
else()
    target_link_libraries(ShadyRay ${SDL2_LIB})
endif()

find_library(TBB_LIB tbb)
if(NOT TBB_LIB)
    message(WARNING "Intel Threaded Building Blocks library not found explicitly. Build may still work if the linker can find it.")
    target_link_libraries(ShadyRay tbb)
else()
    target_link_libraries(ShadyRay ${TBB_LIB})
endif()

find_library(LUA_LIB lua)
if(NOT LUA_LIB)
    message(WARNING "Lua 5.3.x library not found explicitly. Build may still work if the linker can find it.")
    target_link_libraries(ShadyRay lua)
else()
    target_link_libraries(ShadyRay ${LUA_LIB})
endif()

target_link_libraries(ShadyRay dl)
