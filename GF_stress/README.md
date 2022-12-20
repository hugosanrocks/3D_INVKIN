# multiibemdwn3dgeneral

![Stress tensor illustration](https://github.com/hugosanrocks/hugosanrocks.github.com/blob/master/assets/img/stress_tensor.jpg?raw=true)

This sub-repository contains the main code to perform wave propagation simulations inside layered media using the IBEM-DWN method following [Perton et al. (2016)](https://academic.oup.com/gji/article/205/3/1832/657753?login=false). Some auxiliary post-processing and I/O tools are provided in order to write the output files (stress tensor) in a format which is compatible with the forward and inverse modeling code 3D_INVKIN.

The following list presents the subdirectories of this repository

* [multiibem3dgeneral](https://github.com/hugosanrocks/3D_INVKIN/tree/main/GF_stress/multiibemdwn3dgeneral): All Matlab codes used to perform wave propagation simulations and stress tensor estimations.
* [ins](https://github.com/hugosanrocks/3D_INVKIN/tree/main/GF_stress/ins): Folder with input files: `sta.dat`, `fault.dat`, `vmodelV2.dat`. The `sta.dat` contains station locations where the stress tensor has to be computed. `fault.dat` has the subfault locations where the uni-lateral forces have to be applied. `vmodelV2.dat` sets the layered velocity model.
* [out](https://github.com/hugosanrocks/3D_INVKIN/tree/main/GF_stress/out): Folder where the output files will be saved. The important output files are named after the six indempendent stress tensor components and the direction on which the uni-lateral force has been applied. For the six compoents of the stress tensor given a force applied following the "+X" direction we have: `SIGMA_XX_C1`, `SIGMA_YY_C1`, `SIGMA_ZZ_C1`, `SIGMA_XY_C1`,`SIGMA_XZ_C1`, `SIGMA_YZ_C1`.
