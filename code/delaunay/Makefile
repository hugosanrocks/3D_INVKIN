#
#---------------------------------------------------------------
#
#	Makefile for simple example driver program 
#
#  Makefile for 2-D NN-routines for calculating the delaunay triangulation
#  and locating teh Voronoi cell containing a point.
#
F77 = mpif90 #ifort 
CC = gcc -w 
#
#---------------------------------------------------------------

all  : nn qhull driver 


driver  : test.f90 
	$(F77) test.f90 -L./ -lnn2d -lqh -o test


nn	:: 
	$(F77) -c nn/del_sub.f
	$(F77) -c nn/nn.f
	$(F77) -c nn/delaun.f
	$(F77) -c nn/nnplot.f
	$(CC) -c nn/stack.c
	$(CC) -c nn/stackpair.c
	$(CC) -c nn/volume.c
	$(CC) -c nn/utils.c
	ar -r libnn2d.a  del_sub.o nn.o delaun.o nnplot.o stack.o stackpair.o volume.o utils.o
	\rm del_sub.o nn.o delaun.o nnplot.o stack.o stackpair.o volume.o utils.o

qhull	:: 
	$(CC) -c qhull/geom.c
	$(CC) -c qhull/globals.c
	$(CC) -c qhull/io.c
	$(CC) -c qhull/poly.c
	$(CC) -c qhull/qhull.c
	$(CC) -c qhull/set.c
	$(CC) -c qhull/mem.c
	ar -r libqh.a  geom.o globals.o io.o poly.o qhull.o set.o mem.o
	\rm geom.o globals.o io.o poly.o qhull.o set.o mem.o

clean   ::
	\rm -f *.o *.a nn/*.o qhull/*.o 


