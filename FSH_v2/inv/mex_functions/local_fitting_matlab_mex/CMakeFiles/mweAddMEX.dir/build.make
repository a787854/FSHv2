# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/yanglei/Downloads/mweCmakeMexCppCuda-master

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/yanglei/Downloads/mweCmakeMexCppCuda-master

# Include any dependencies generated for this target.
include CMakeFiles/mweAddMEX.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mweAddMEX.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mweAddMEX.dir/flags.make

CMakeFiles/cuda_compile_1.dir/src/lib/cuda/cuda_compile_1_generated_gpuadd.cu.o: CMakeFiles/cuda_compile_1.dir/src/lib/cuda/cuda_compile_1_generated_gpuadd.cu.o.depend
CMakeFiles/cuda_compile_1.dir/src/lib/cuda/cuda_compile_1_generated_gpuadd.cu.o: CMakeFiles/cuda_compile_1.dir/src/lib/cuda/cuda_compile_1_generated_gpuadd.cu.o.cmake
CMakeFiles/cuda_compile_1.dir/src/lib/cuda/cuda_compile_1_generated_gpuadd.cu.o: src/lib/cuda/gpuadd.cu
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/yanglei/Downloads/mweCmakeMexCppCuda-master/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building NVCC (Device) object CMakeFiles/cuda_compile_1.dir/src/lib/cuda/cuda_compile_1_generated_gpuadd.cu.o"
	cd /home/yanglei/Downloads/mweCmakeMexCppCuda-master/CMakeFiles/cuda_compile_1.dir/src/lib/cuda && /usr/bin/cmake -E make_directory /home/yanglei/Downloads/mweCmakeMexCppCuda-master/CMakeFiles/cuda_compile_1.dir/src/lib/cuda/.
	cd /home/yanglei/Downloads/mweCmakeMexCppCuda-master/CMakeFiles/cuda_compile_1.dir/src/lib/cuda && /usr/bin/cmake -D verbose:BOOL=$(VERBOSE) -D build_configuration:STRING= -D generated_file:STRING=/home/yanglei/Downloads/mweCmakeMexCppCuda-master/CMakeFiles/cuda_compile_1.dir/src/lib/cuda/./cuda_compile_1_generated_gpuadd.cu.o -D generated_cubin_file:STRING=/home/yanglei/Downloads/mweCmakeMexCppCuda-master/CMakeFiles/cuda_compile_1.dir/src/lib/cuda/./cuda_compile_1_generated_gpuadd.cu.o.cubin.txt -P /home/yanglei/Downloads/mweCmakeMexCppCuda-master/CMakeFiles/cuda_compile_1.dir/src/lib/cuda/cuda_compile_1_generated_gpuadd.cu.o.cmake

CMakeFiles/mweAddMEX.dir/src/mex/mainMex.cpp.o: CMakeFiles/mweAddMEX.dir/flags.make
CMakeFiles/mweAddMEX.dir/src/mex/mainMex.cpp.o: src/mex/mainMex.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanglei/Downloads/mweCmakeMexCppCuda-master/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/mweAddMEX.dir/src/mex/mainMex.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mweAddMEX.dir/src/mex/mainMex.cpp.o -c /home/yanglei/Downloads/mweCmakeMexCppCuda-master/src/mex/mainMex.cpp

CMakeFiles/mweAddMEX.dir/src/mex/mainMex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mweAddMEX.dir/src/mex/mainMex.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanglei/Downloads/mweCmakeMexCppCuda-master/src/mex/mainMex.cpp > CMakeFiles/mweAddMEX.dir/src/mex/mainMex.cpp.i

CMakeFiles/mweAddMEX.dir/src/mex/mainMex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mweAddMEX.dir/src/mex/mainMex.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanglei/Downloads/mweCmakeMexCppCuda-master/src/mex/mainMex.cpp -o CMakeFiles/mweAddMEX.dir/src/mex/mainMex.cpp.s

CMakeFiles/mweAddMEX.dir/src/lib/add.cpp.o: CMakeFiles/mweAddMEX.dir/flags.make
CMakeFiles/mweAddMEX.dir/src/lib/add.cpp.o: src/lib/add.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanglei/Downloads/mweCmakeMexCppCuda-master/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/mweAddMEX.dir/src/lib/add.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mweAddMEX.dir/src/lib/add.cpp.o -c /home/yanglei/Downloads/mweCmakeMexCppCuda-master/src/lib/add.cpp

CMakeFiles/mweAddMEX.dir/src/lib/add.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mweAddMEX.dir/src/lib/add.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanglei/Downloads/mweCmakeMexCppCuda-master/src/lib/add.cpp > CMakeFiles/mweAddMEX.dir/src/lib/add.cpp.i

CMakeFiles/mweAddMEX.dir/src/lib/add.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mweAddMEX.dir/src/lib/add.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanglei/Downloads/mweCmakeMexCppCuda-master/src/lib/add.cpp -o CMakeFiles/mweAddMEX.dir/src/lib/add.cpp.s

# Object files for target mweAddMEX
mweAddMEX_OBJECTS = \
"CMakeFiles/mweAddMEX.dir/src/mex/mainMex.cpp.o" \
"CMakeFiles/mweAddMEX.dir/src/lib/add.cpp.o"

# External object files for target mweAddMEX
mweAddMEX_EXTERNAL_OBJECTS = \
"/home/yanglei/Downloads/mweCmakeMexCppCuda-master/CMakeFiles/cuda_compile_1.dir/src/lib/cuda/cuda_compile_1_generated_gpuadd.cu.o"

lib/mweAddMEX.mexa64: CMakeFiles/mweAddMEX.dir/src/mex/mainMex.cpp.o
lib/mweAddMEX.mexa64: CMakeFiles/mweAddMEX.dir/src/lib/add.cpp.o
lib/mweAddMEX.mexa64: CMakeFiles/cuda_compile_1.dir/src/lib/cuda/cuda_compile_1_generated_gpuadd.cu.o
lib/mweAddMEX.mexa64: CMakeFiles/mweAddMEX.dir/build.make
lib/mweAddMEX.mexa64: /home/yanglei/Polyspace/R2020b/bin/glnxa64/libmex.so
lib/mweAddMEX.mexa64: /home/yanglei/Polyspace/R2020b/bin/glnxa64/libmx.so
lib/mweAddMEX.mexa64: /usr/local/cuda/lib64/libcudart_static.a
lib/mweAddMEX.mexa64: /usr/lib/x86_64-linux-gnu/librt.so
lib/mweAddMEX.mexa64: /usr/local/cuda/lib64/libcublas.so
lib/mweAddMEX.mexa64: /usr/local/cuda/lib64/libcufft.so
lib/mweAddMEX.mexa64: Matlabdef.def
lib/mweAddMEX.mexa64: CMakeFiles/mweAddMEX.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/yanglei/Downloads/mweCmakeMexCppCuda-master/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX shared library lib/mweAddMEX.mexa64"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mweAddMEX.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mweAddMEX.dir/build: lib/mweAddMEX.mexa64

.PHONY : CMakeFiles/mweAddMEX.dir/build

CMakeFiles/mweAddMEX.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mweAddMEX.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mweAddMEX.dir/clean

CMakeFiles/mweAddMEX.dir/depend: CMakeFiles/cuda_compile_1.dir/src/lib/cuda/cuda_compile_1_generated_gpuadd.cu.o
	cd /home/yanglei/Downloads/mweCmakeMexCppCuda-master && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yanglei/Downloads/mweCmakeMexCppCuda-master /home/yanglei/Downloads/mweCmakeMexCppCuda-master /home/yanglei/Downloads/mweCmakeMexCppCuda-master /home/yanglei/Downloads/mweCmakeMexCppCuda-master /home/yanglei/Downloads/mweCmakeMexCppCuda-master/CMakeFiles/mweAddMEX.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mweAddMEX.dir/depend

