#complile
ifort test_delaunay_mod.f90 delaunay_mod.f90 -o test

#run
./test

#plot
octave plot_mesh_delaunay.m

#position of nodes inside peruchile.xyz

#final vertex of triangles in triangle_vert.out
