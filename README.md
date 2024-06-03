[forest height mosaic for northeastern China]<img src="https://github.com/a787854/FSHv2/blob/dev/preview_china_FSH_mosaic.jpeg" width="100%">
![Forest Height Mosaic for northeastern China]([image_url](https://github.com/a787854/FSHv2/blob/dev/preview_china_FSH_mosaic.jpeg))
# Forest Stand Height (FSH) v2  
This software is developed to generate the forest height mosaic by using ALOS-1 or ALOS-2 repeat-pass SAR Interferometry and GEDI spaceborne LiDAR data after InSAR preprocessing by JPL's ISCE software. This project is the extended version of original [FSH](https://github.com/leiyangleon/FSH) scripts.


## Requirements
* Python (>=3.6) Enviroment (mandatory)
* [MATLAB](https://de.mathworks.com/) (mandatory)
* [CMake](https://cmake.org/) (mandatory)
* [CUDA](https://developer.nvidia.com/cuda-downloads) (mandatory)


## Tested environment:
Ubuntu 20.04 and 22.04; CMake 3.16.3; CUDA 12.1; MATLAB R2020b

## Installation
You have to compile the mex related files to obtain mex functions based on [CMake template] (https://github.com/BjoernHaefner/mweCmakeMexCppCuda):
1) go to the dectories FSH_v2/inv/mex_functions/global_fitting_matlab_mex and FSH_v2/inv/mex_functions/local_fitting_matlab_mex, and do   
   `cmake ..`  
   `make`    
   `the generated .mexa64 file is located in the .../build/libs, and put them into the place where MATLAB terminal can invoke`

## How to use:
An example dataset and a quick start guide is provided [here](https://github.com/a787854/FSHv2/edit/dev/quick_start_example.md).

## Relevant products:
Two forest height mosaics were generated in the northeast of U.S. and China using these scripts and can be downloaded via [this link](https://doi.org/10.5281/zenodo.11142224)

   



## Contributors:
Yanghai Yu (Email: yuyanghai@nssc.ac.cn, or filippoyu0717@gmail.com)        
Yang Lei (Email: leiyang@nssc.ac.cn, leiyangfrancis@gmail.com)         
Paul Siqueira 







 




