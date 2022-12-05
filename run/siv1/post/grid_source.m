clear all
close all
clc

%load data
source = load('');
%time samples
nt=39;

%number of nodes
[siz,a] = size(source);
nnodes = siz/nt;

%separate slip-rate functions per node
source_r = reshape(source,nt,nnodes);

%read grid points
grid = load('dat/grid.in');
x = grid(:,1);
y = grid(:,2);
v = source_r(7,:)';

%points where to interpolate
x_want = 500:100:34250;
y_want = 750:100:16250;
[xq,yq]= meshgrid(x_want,y_want);
v_want = griddata(x,y,v,xq,yq);


