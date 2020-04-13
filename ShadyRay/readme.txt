ShadyRay Raytracer
------------------

Build instructions
------------------
ShadyRay uses CMake as the build system. There is a CMakeLists.txt file in the main directory. Create a build
directory and run CMake just as you would for any other program.

Currently, it supports building only on Linux systems. Support for Windows is complicated by some of the libraries used.
All library licenses are as specified on their GitHub repositories. This project uses the dkm library from GitHub user
genbattle.

The dependencies required to build ShadyRay are:
    1. Intel Embree 2.7.x : get this from Embree's GitHub repo. IMPORTANT: Build Embree from source with support for
    ray masking! ShadyRay will not work without this.

    2. Intel Threaded Building Blocks : get this from Intel Parallel Studio or from your distro's package manager.

    3. Assimp 3.3 or above : Grab the latest version of the code from GitHub or SourceForge. IMPORTANT: Assimp 4
    will not work!

    4. Lua 5.3.x : get it from the Lua repositories.

    5. dkm : get it from genbattle's GitHub repository.

    6. SDL2 : get the latest development build files from the SDL2 website.

How to use it
-------------
ShadyRay uses a Lua script to set up parameters for a scene. Consult the sample scene.lua file on how to configure a
scene. Once you have a Lua file, just run the executable with the path to the Lua script as the ONLY command line arg.