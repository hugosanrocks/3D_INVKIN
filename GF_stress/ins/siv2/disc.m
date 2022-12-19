clear all
close all
clc

f=load('inv2_.dat');
nx=81;
ny=41;


k=1;
for i=2:2:ny
 for j=2:2:nx
  fault(k,:) = f((i-1)*nx+j,:);
  k=k+1;
 end
end

save('-ascii','fault_2.dat','fault');
