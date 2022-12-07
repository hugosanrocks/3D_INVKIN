# 3D_INKIN

This repository contains the main code to perform a kinematic inversion model following Sánchez-Reyes et al. (2018). Some auxiliary plotting tools and one example (Source Inversion Validation case 1) are supplied.

The following list presents the subdirectories of this repository

* [code](https://github.com/hugosanrocks/3D_INVKIN/tree/main/code):
  + adjoint
  + bin
  + covariance/
  + delaunay/
  + filtro/
  + forward/
  + general/
  + hessian/
  + include/
  + install/
  + inversion/
  + obj/
  + other_tools
  + plotting_tools/
  + preprocess/
  + prior/
  + sensitivity/
  + uncertain/
  + windows/
* [run](https://github.com/hugosanrocks/3D_INVKIN/tree/main/run):
  + [siv1](https://github.com/hugosanrocks/3D_INVKIN/tree/main/run/siv1)

## Installation

### Requirements:

This version of 3D_INVKIN requires to be compiled with ifort and the main code has to be linked to the intel MKL library. It is suggested to download and install the latest ifort version using the FREE OneAPI interl distribution of ifort and MKL. You can use the following links to install this tools once you have your free intel account.

* [ifort](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran) 
* [MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html):

Once you have installed the ifort compiler and the MKL library, you will have to compile the OPTIMIZATION libray from [SEISCOPE](https://seiscope2.osug.fr/SEISCOPE-OPTIMIZATION-TOOLBOX) and the FFT library [FFTW3](https://www.fftw.org/download.html). A version of these two libraries is already provided inside the main repository 3D_INVKIN under the subdirectory [code>other_tools](https://github.com/hugosanrocks/3D_INVKIN/tree/main/code/other_tools).

Please follow the next steps in order to correctly compile the OPTIMIZATION and FFTW3 libraries.

## OPTIMIZATION TOOLS_BOX:

" $cd code/other_tools/TOOLS_BOX/trunk/OPTIMIZATION/


