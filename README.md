# 3D_INKIN


![Figure text](https://github.com/hugosanrocks/hugosanrocks.github.com/blob/master/assets/img/model_time_space_windows.jpg?raw=true)

This repository contains the main code to perform a kinematic inversion model following Sánchez-Reyes et al. (2018). Some auxiliary plotting tools and one example (Source Inversion Validation case 1) are supplied.

The following list presents the subdirectories of this repository

* [GF_stress](https://github.com/hugosanrocks/3D_INVKIN/tree/main/GF_stress): Matlab code used to compute (previous to any forward or inverse modeling) the pseudo Green's functions. This code computes the six independent components of the stress tensor at every receiver location given three different uni-axial forces (fx, fy, fz). The input files must be inside `GF_stress/ins` folder, while the output files are written inside `GF_stress/out` folder. The main code is inside `GF_stress/multiibemdwn3dgeneral`. More information about the use of this code is given inside the main code folder.
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

#### Plotting tools

The post processing steps after performing a direct forward modeling or a kinematic source inversion are performed using [GMT](https://www.generic-mapping-tools.org/) (Generic Mapping Tools) and [CWP](https://wiki.seismic-unix.org/start) (Seismic Unix SU). Please be sure to have install GMT and CWP in order to run the plotting tools inside the `code/plotting_tools`. Here below you will find some steps that you can follow to install GMT 6.4 and CWP

##### GMT 6.4 Installation

The following are steps to install GMT 6.4 under a conda (Anaconda 3) environment.

First, define the channels where to look for available software:

`conda config --prepend channels conda-forge`

Now, create a conda environment with all the required software:

`conda create --name gmtenv python=3.9 numpy pandas xarray netcdf4 packaging gmt`

That's all, now you have an environment with GMT 6 installed. To use it, activate the environment:

`conda activate gmtenv`

##### CWP Installation

First, download the CWP version that you want:

[CWP Download site](https://nextcloud.seismic-unix.org/index.php/s/LZpzc8jMzbWG9BZ)

Once downloaded, you can follow the instructions given [here](https://wiki.seismic-unix.org/sudoc:su_installation) to correctly install CWP on your precise system. The next are generic steps that can help you to rapidly install CWP if you are on Ubuntu 20.04.

Create the directory for installation and give user permissions:

`$sudo mkdir /usr/local/cwp`

`$sudo chown USERNAME /usr/local/cwp`

Include the next lines inside your `.bashrc` file:

`# CWP`

`export CWPROOT=/usr/local/cwp`

`# set the PATH variable`

`export PATH=$PATH:$CWPROOT/bin`

Source the new `.bashrc` file

`source ~/.bashrc`

Go to the installation folder and untar the downloaded files

`$cd /usr/local/cwp`

`$tar -xvf Downloads/cwp_su_all_44R26.tgz`

Configure the installation files and compile everything

`$cp configs/Makefile.config_Linux_Ubuntu_20.04 /usr/local/cwp/src/Makefile.config`

`$make install`



### Install 3D_INVKIN

The most important thing to do in order to successfully compile 3D_INVKIN is to correctly set the paths for all the libraries. To do that, please modify or create your own `make.user.inc` file. Some examples are available at the location `code/install/`. Remember that the compiler options have to be set to "ifort".

Once you have correctly set the paths, 3D_INVKIN has to be compiled:

`$cd code/obj/`

`$make all`

Other partial options of compilation are available:

`$make PREPROCESS`: builds bin/PREPROCESS code that is used to prepared the Green's function bank.

`$make FORWARD`: builds bin/FORWARD code that only estimates synthetic seismograms given a source model and the Green's function bank.

`$make INV3DKIN`: builds bin/INV3DKIN code that performs the kinematic inversion. Forward and inverse problem are performed given the observations, inversion options and Green function bank.

`$make HESSIAN`: builds bin/HESSIAN code that can be used to build an approximation of the Hessian (not fully tested).

`$make SENSITIVITY`: builds bin/SENSITIVITY which can be used for sensitivity analysis (under construction still)

`$make objs`: builds objectcs and modules

`$make clean`: cleans all the binary files


## Example (SIV1)


Due to the large size of the files corresponding to the Green functions, it is necessary to downlothese files from another site besides GitHub. The preprocessed Green function bank can be downloaded from the following link:

SIV1 Preprocessed Green's function bank: [Download](https://www.dropbox.com/scl/fo/6kno3ar6ukxu1y6oc2b7i/h?dl=0&rlkey=cbdlfv5sbv0gxu2a6s9wkum9z)

It is necessary to download the file `TRACT_time.bin` and place it inside the directory `run/siv1/dat`.

Once the Green's function bank is placed in the correct directory, the scripts forward.sh and inversion.sh can be used to perform a forward modeling or a kinematic inversion.

NEXT STEPS TO RUN PROGRESSIVE INVERSIONS WILL BE SOON UPLOADED.



