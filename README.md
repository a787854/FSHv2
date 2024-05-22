# Forest Stand Height (FSH) v2  
This software is developed to generate the forest height mosaic by using ALOS-1 or ALOS-2 repeat-pass SAR Interferometry and GEDI spaceborne LiDAR data after InSAR preprocessing by JPL's ISCE software. This project is the extended version of original [FSH](https://github.com/leiyangleon/FSH) scripts.


## Requirements
* Python (>=3.6) Enviroment (mandatory)
* [MATLAB](https://de.mathworks.com/) (mandatory)
* [CMake](https://cmake.org/) (mandatory)
* [CUDA](https://developer.nvidia.com/cuda-downloads) (mandatory)


## Tested environment:
Ubuntu 20.04 and 22.04; CMake 3.16.3; CUDA 12.1; MATLAB R2020b

## Getting started
You have to compile the mex related files to obtain mex functions based on [CMake template] (https://github.com/BjoernHaefner/mweCmakeMexCppCuda):
1) go to the dectories FSH_v2/inv/mex_functions/global_fitting_matlab_mex and FSH_v2/inv/mex_functions/local_fitting_matlab_mex, and do
   `cmake ..`

   `make`

   `the generated .mexa64 file is located in the .../build/libs, and put them into the place where MATLAB can invoke`





 




