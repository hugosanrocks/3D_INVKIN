#!/bin/bash
#
gfortran -c -fno-underscoring fftest.f90
if [ $? -ne 0 ]; then
  echo "Errors compiling fftest.f90."
  exit
fi
#
gfortran fftest.o -lfftw3f -lm -lc -o fftest
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fftest.o."
  exit
fi
#
#

./fftest > fft_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fftest."
  exit
fi
#
echo "Program output written to fft_output.txt"
