========================================
  3D_INVKIN SIV1 Green's function bank
========================================

This folder contains the following files

TRACT_time.bin   242 MB

This file contains the traction time 
history for 648 subfaults (36 strike X
18 dip) that where estimated using the
40 station locations from SIV1 excercise
as uniaxial point sources (along X, Y and
Z). 

The file contains 
648 x 271 x 40 x 3 x 3 x 4 = 252875520 = 242 MB 
 |     |    |    |   |   |
 |     |    |    |   | bytes
 |     |    |    | components
 |     |    |  forces
 |     |   stations
 |   time samples
Subfaults

In order to be able to use this file 
for computing synthetic seismograms,
or running the kinematic inversion,
this file has to be inside the "dat"
directory.



