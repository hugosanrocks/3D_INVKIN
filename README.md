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

Do not forget to correctly set the environment variables after the installation. For instance:

`$source /home/user/intel/oneapi/setvars.sh`

or export the paths and variables inside your local .bashrc

Please follow the next steps in order to correctly compile the OPTIMIZATION and FFTW3 libraries.

#### OPTIMIZATION TOOLS_BOX:

First, modify the file `Makefile.inc` in such a way that ifort and icc are the default compiler options. This file is located at `code/other_tools/TOOLS_BOX/trunk/Makefile.inc`

Then, go to the correct path and compile everything and build the libraries:

`$cd code/other_tools/TOOLS_BOX/trunk/OPTIMIZATION/`

`$make`

`$make lib`

#### FFTW3

Go to the correct path and configure all the files before the compilation. 

`$cd code/other_tools/fftw-3.3.10/`

`$chmod +x configure`

`$./configure --enable-single --enable-shared`

`$sudo make`

`$sudo make install`


More information about the installation can be found at the FFTW3 official website. It is needed to copile fftw3 with the option `--enable-single` to be compatigle with 3D_INVKIN.

### Install 3D_INVKIN



