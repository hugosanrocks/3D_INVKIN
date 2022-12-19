% D. A. Dunavant, 
% High degree efficient symmetrical Gaussian quadrature rules for the triangle, 
% Int. J. Num. Meth. Engng, 21(1985):1129-1148.

 % n = 3
para.cubature ...
   = [0.33333333333333 0.33333333333333 -0.56250000000000 ;
      0.20000000000000 0.20000000000000 0.52083333333333 ; 
      0.20000000000000 0.60000000000000 0.52083333333333 ;
      0.60000000000000 0.20000000000000 0.52083333333333];

para.gaussian.ngau  = size(para.cubature,1);