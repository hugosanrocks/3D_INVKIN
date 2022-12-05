rm *_C1 *_C2 *_C3

gfortran filemanage.f90 -o file1
gfortran filemanage2.f90 -o file2
gfortran filemanage3.f90 -o file3
./file1
./file2
./file3
