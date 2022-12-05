ifort matrix_svd.f90 -L/opt/intel/compilers_and_libraries_2017.1.132/linux/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -qopenmp -lmkl_lapack95_lp64 -o matrix
#rm matrix ~/MEGA/INV3DKIN/run/toysmall8x4/matrix
cp matrix ~/MEGA/INV3DKIN/run/toysmall8x4/
