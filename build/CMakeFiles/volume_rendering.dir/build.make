# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = E:/CMake/bin/cmake.exe

# The command to remove a file.
RM = E:/CMake/bin/cmake.exe -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:/MinGW/CG/Tetrahedron-Volume-Rendering

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:/MinGW/CG/Tetrahedron-Volume-Rendering/build

# Include any dependencies generated for this target.
include CMakeFiles/volume_rendering.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/volume_rendering.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/volume_rendering.dir/flags.make

CMakeFiles/volume_rendering.dir/src/compositor.cpp.obj: CMakeFiles/volume_rendering.dir/flags.make
CMakeFiles/volume_rendering.dir/src/compositor.cpp.obj: CMakeFiles/volume_rendering.dir/includes_CXX.rsp
CMakeFiles/volume_rendering.dir/src/compositor.cpp.obj: ../src/compositor.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:/MinGW/CG/Tetrahedron-Volume-Rendering/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/volume_rendering.dir/src/compositor.cpp.obj"
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/volume_rendering.dir/src/compositor.cpp.obj -c D:/MinGW/CG/Tetrahedron-Volume-Rendering/src/compositor.cpp

CMakeFiles/volume_rendering.dir/src/compositor.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/volume_rendering.dir/src/compositor.cpp.i"
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:/MinGW/CG/Tetrahedron-Volume-Rendering/src/compositor.cpp > CMakeFiles/volume_rendering.dir/src/compositor.cpp.i

CMakeFiles/volume_rendering.dir/src/compositor.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/volume_rendering.dir/src/compositor.cpp.s"
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:/MinGW/CG/Tetrahedron-Volume-Rendering/src/compositor.cpp -o CMakeFiles/volume_rendering.dir/src/compositor.cpp.s

CMakeFiles/volume_rendering.dir/src/main.cpp.obj: CMakeFiles/volume_rendering.dir/flags.make
CMakeFiles/volume_rendering.dir/src/main.cpp.obj: CMakeFiles/volume_rendering.dir/includes_CXX.rsp
CMakeFiles/volume_rendering.dir/src/main.cpp.obj: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:/MinGW/CG/Tetrahedron-Volume-Rendering/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/volume_rendering.dir/src/main.cpp.obj"
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/volume_rendering.dir/src/main.cpp.obj -c D:/MinGW/CG/Tetrahedron-Volume-Rendering/src/main.cpp

CMakeFiles/volume_rendering.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/volume_rendering.dir/src/main.cpp.i"
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:/MinGW/CG/Tetrahedron-Volume-Rendering/src/main.cpp > CMakeFiles/volume_rendering.dir/src/main.cpp.i

CMakeFiles/volume_rendering.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/volume_rendering.dir/src/main.cpp.s"
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:/MinGW/CG/Tetrahedron-Volume-Rendering/src/main.cpp -o CMakeFiles/volume_rendering.dir/src/main.cpp.s

CMakeFiles/volume_rendering.dir/src/opticsData.cpp.obj: CMakeFiles/volume_rendering.dir/flags.make
CMakeFiles/volume_rendering.dir/src/opticsData.cpp.obj: CMakeFiles/volume_rendering.dir/includes_CXX.rsp
CMakeFiles/volume_rendering.dir/src/opticsData.cpp.obj: ../src/opticsData.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:/MinGW/CG/Tetrahedron-Volume-Rendering/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/volume_rendering.dir/src/opticsData.cpp.obj"
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/volume_rendering.dir/src/opticsData.cpp.obj -c D:/MinGW/CG/Tetrahedron-Volume-Rendering/src/opticsData.cpp

CMakeFiles/volume_rendering.dir/src/opticsData.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/volume_rendering.dir/src/opticsData.cpp.i"
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:/MinGW/CG/Tetrahedron-Volume-Rendering/src/opticsData.cpp > CMakeFiles/volume_rendering.dir/src/opticsData.cpp.i

CMakeFiles/volume_rendering.dir/src/opticsData.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/volume_rendering.dir/src/opticsData.cpp.s"
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:/MinGW/CG/Tetrahedron-Volume-Rendering/src/opticsData.cpp -o CMakeFiles/volume_rendering.dir/src/opticsData.cpp.s

CMakeFiles/volume_rendering.dir/src/renderer.cpp.obj: CMakeFiles/volume_rendering.dir/flags.make
CMakeFiles/volume_rendering.dir/src/renderer.cpp.obj: CMakeFiles/volume_rendering.dir/includes_CXX.rsp
CMakeFiles/volume_rendering.dir/src/renderer.cpp.obj: ../src/renderer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:/MinGW/CG/Tetrahedron-Volume-Rendering/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/volume_rendering.dir/src/renderer.cpp.obj"
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/volume_rendering.dir/src/renderer.cpp.obj -c D:/MinGW/CG/Tetrahedron-Volume-Rendering/src/renderer.cpp

CMakeFiles/volume_rendering.dir/src/renderer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/volume_rendering.dir/src/renderer.cpp.i"
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:/MinGW/CG/Tetrahedron-Volume-Rendering/src/renderer.cpp > CMakeFiles/volume_rendering.dir/src/renderer.cpp.i

CMakeFiles/volume_rendering.dir/src/renderer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/volume_rendering.dir/src/renderer.cpp.s"
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:/MinGW/CG/Tetrahedron-Volume-Rendering/src/renderer.cpp -o CMakeFiles/volume_rendering.dir/src/renderer.cpp.s

CMakeFiles/volume_rendering.dir/src/volume.cpp.obj: CMakeFiles/volume_rendering.dir/flags.make
CMakeFiles/volume_rendering.dir/src/volume.cpp.obj: CMakeFiles/volume_rendering.dir/includes_CXX.rsp
CMakeFiles/volume_rendering.dir/src/volume.cpp.obj: ../src/volume.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:/MinGW/CG/Tetrahedron-Volume-Rendering/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/volume_rendering.dir/src/volume.cpp.obj"
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/volume_rendering.dir/src/volume.cpp.obj -c D:/MinGW/CG/Tetrahedron-Volume-Rendering/src/volume.cpp

CMakeFiles/volume_rendering.dir/src/volume.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/volume_rendering.dir/src/volume.cpp.i"
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:/MinGW/CG/Tetrahedron-Volume-Rendering/src/volume.cpp > CMakeFiles/volume_rendering.dir/src/volume.cpp.i

CMakeFiles/volume_rendering.dir/src/volume.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/volume_rendering.dir/src/volume.cpp.s"
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:/MinGW/CG/Tetrahedron-Volume-Rendering/src/volume.cpp -o CMakeFiles/volume_rendering.dir/src/volume.cpp.s

# Object files for target volume_rendering
volume_rendering_OBJECTS = \
"CMakeFiles/volume_rendering.dir/src/compositor.cpp.obj" \
"CMakeFiles/volume_rendering.dir/src/main.cpp.obj" \
"CMakeFiles/volume_rendering.dir/src/opticsData.cpp.obj" \
"CMakeFiles/volume_rendering.dir/src/renderer.cpp.obj" \
"CMakeFiles/volume_rendering.dir/src/volume.cpp.obj"

# External object files for target volume_rendering
volume_rendering_EXTERNAL_OBJECTS =

volume_rendering.exe: CMakeFiles/volume_rendering.dir/src/compositor.cpp.obj
volume_rendering.exe: CMakeFiles/volume_rendering.dir/src/main.cpp.obj
volume_rendering.exe: CMakeFiles/volume_rendering.dir/src/opticsData.cpp.obj
volume_rendering.exe: CMakeFiles/volume_rendering.dir/src/renderer.cpp.obj
volume_rendering.exe: CMakeFiles/volume_rendering.dir/src/volume.cpp.obj
volume_rendering.exe: CMakeFiles/volume_rendering.dir/build.make
volume_rendering.exe: CMakeFiles/volume_rendering.dir/linklibs.rsp
volume_rendering.exe: CMakeFiles/volume_rendering.dir/objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:/MinGW/CG/Tetrahedron-Volume-Rendering/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable volume_rendering.exe"
	E:/CMake/bin/cmake.exe -E rm -f CMakeFiles/volume_rendering.dir/objects.a
	D:/MinGW/bin/ar.exe cr CMakeFiles/volume_rendering.dir/objects.a @CMakeFiles/volume_rendering.dir/objects1.rsp
	D:/MinGW/bin/x86_64-w64-mingw32-g++.exe  -fopenmp -g    -Wl,--whole-archive CMakeFiles/volume_rendering.dir/objects.a -Wl,--no-whole-archive  -o volume_rendering.exe -Wl,--out-implib,libvolume_rendering.dll.a -Wl,--major-image-version,0,--minor-image-version,0 @CMakeFiles/volume_rendering.dir/linklibs.rsp

# Rule to build all files generated by this target.
CMakeFiles/volume_rendering.dir/build: volume_rendering.exe

.PHONY : CMakeFiles/volume_rendering.dir/build

CMakeFiles/volume_rendering.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/volume_rendering.dir/cmake_clean.cmake
.PHONY : CMakeFiles/volume_rendering.dir/clean

CMakeFiles/volume_rendering.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" D:/MinGW/CG/Tetrahedron-Volume-Rendering D:/MinGW/CG/Tetrahedron-Volume-Rendering D:/MinGW/CG/Tetrahedron-Volume-Rendering/build D:/MinGW/CG/Tetrahedron-Volume-Rendering/build D:/MinGW/CG/Tetrahedron-Volume-Rendering/build/CMakeFiles/volume_rendering.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/volume_rendering.dir/depend

