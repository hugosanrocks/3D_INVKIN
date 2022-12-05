\rm ../pre_pmcl3d_kine

ifort pre_pmcl3d_kine.f90 geometry.f90 subcell.f90 xapiir.f -o ../pre_pmcl3d_kine_3h

#gfortran pre_pmcl3d_kine.f90 geometry.f90 subcell.f90 xapiir.f -o ../../bin/pre_pmcl3d_kine
