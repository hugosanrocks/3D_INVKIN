clear all
close all
clc

      A(1:6,1) = [ 7.52,-0.76, 5.13,-4.75, 1.33,-2.40];
      A(1:6,2) = [-1.10, 0.62, 6.62, 8.52, 4.91,-6.77];
      A(1:6,3) = [-7.95, 9.34,-5.66, 5.75,-5.49, 2.34];
      A(1:6,4) = [ 1.08,-7.10, 0.87, 5.30,-3.52, 3.95];

[U,S,V] = svd(A)