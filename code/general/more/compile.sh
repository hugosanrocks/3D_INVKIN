#!/bin/bash
#
gfortran -c -fno-underscoring fftw3_prb.f -o test.o
if [ $? -ne 0 ]; then
  echo "Errors compiling test.f."
  exit
fi
#
gfortran test.o -lfftw3 -lm -lc -o test
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test.o."
  exit
fi
#
#
./test > test_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test."
  exit
fi
#
echo "Program output written to test_output.txt"
