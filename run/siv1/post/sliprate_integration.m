clear all
close all
clc

%clean output files
system('rm file_*.bin');

opt = load('opt_grid.dat');

%option = 1 (regular grid), 2 (iregular grid)
%option = input('grid option (1 regular) (2 iregular):')
option = opt(1)

%input parameters
%tmax = input('Maximum time for integration (sec): ');
tmax = opt(2)
%dt = input('dt of files: ');
dt = opt(3)

%nodes along dip and strike
if ( option == 1 )
 ly = 18;
 lx = 36;
 %load files for target and inversion results
 %fileinv=input('inversion input file: ');
 %filetar=input('target input file: ');
 fileinv='source.out'
 filetar='../dat/model_target.dat'
else
 mesh_lines
 ltot = 595;
 %load files for target and inversion results
 %fileinv=input('inversion input file: ');
 %filetar=input('target input file: ');
 fileinv='source.out'
 filetar='../dat/model_target_triang.dat'
end


%total time samples
ntot = floor(tmax/dt)+1
%time axis
t = 0:dt:dt*(ntot-1);

%output files
%fileoutinv=input('output: ');
%fileouttar=input('output: ');
fileouttar='file_1.bin'
fileoutinv='file_2.bin'

%bin=input('bin(1) not bin(2): ');
bin=2;

%load files
if (bin==1)
 fid=fopen(file,'r');
 source_in=fread(fid,'single');
 source=reshape(source_in,39,648);
elseif (bin==2)
 source = load(fileinv);
 target = load(filetar);
else
 tsam=39;
 lx=36; ly=18; 
 source=reshape(weig,lx*ly,tsam)';
end

%regular grid
if ( option == 1 ) 
m = 1;
for i=1:ly
  for j=1:lx
   slipinv = source(:,m);
   sliptar = target(:,m);
     dispoldinv = 0;
     dispoldtar = 0;
     for it=1:ntot
        % Compute displacements (integrate velocities)
        if it > 1
            dispinv(it) = slipinv(it) * dt + dispoldinv;
            disptar(it) = sliptar(it) * dt + dispoldtar;
        else
            dispinv(it) = 0.0;
            disptar(it) = 0.0;
        end
        dispoldinv = dispinv(it);
        dispoldtar = disptar(it);
     end
   disi(:,m) = dispinv(:);
   dist(:,m) = disptar(:);
   maxdispinv(m) = max(disi(:,m));
   maxdisptar(m) = max(dist(:,m));
   m = m+1;
  end
end
%iregular grid
else
m = 1;
for i=1:ltot
   slipinv = source(:,m);
   sliptar = target(:,m);
     dispoldinv = 0;
     dispoldtar = 0;
     for it=1:ntot
        % Compute displacements (integrate velocities)
        if it > 1
            dispinv(it) = slipinv(it) * dt + dispoldinv;
            disptar(it) = sliptar(it) * dt + dispoldtar;
        else
            dispinv(it) = 0.0;
            disptar(it) = 0.0;
        end
        dispoldinv = dispinv(it);
        dispoldtar = disptar(it);
     end
   disi(:,m) = dispinv(:);
   dist(:,m) = disptar(:);
   maxdispinv(m) = max(disi(:,m));
   maxdisptar(m) = max(dist(:,m));
   m = m+1;
end
end


%slip solution and target
if ( option == 1 )
 Bsol = reshape(maxdispinv,lx,ly)';
 Btar = reshape(maxdisptar,lx,ly)';


 %coordinates of slip array
 maxslip=max(max(Bsol))
 %along strike and dip increment
 dx=1;
 dy=1;
 dxint=0.1;
 dyint=0.1;
 x=0:dx:dx*(lx-1);
 y=0:dy:dx*(ly-1);
 %coordinates to interpolate
 xint=x(1):dxint:x(end);
 yint=y(1):dxint:y(end);
 %tables to interpolate
 [X,Y] = meshgrid(x,y);
 [XI,YI] = meshgrid(xint,yint);

 %interpolate to correctly plot
 Bsolint = interp2(X,Y,Bsol,XI,YI);
 Btarint = interp2(X,Y,Btar,XI,YI);
else
 %read grid points
% grid = load('../dat/nodes.in');
% dx = 0.1; dy=0.1;
% dx1= 0.3; dy1=0.3;
% x = grid(:,1);
% y = grid(:,2);
 vinv = maxdispinv';
 vtar = maxdisptar';
 a = load('slipc.out');
 b = load('slipc.out');
 ar= reshape(a(:,1),340,160)';
 br= reshape(a(:,1),340,160)';
 Bsolint = ar;
 Btarint = br;
 %points where to interpolate
% x_want = min(x)+dx1:dx:max(x)-dx1;
% y_want = min(y)+dy1:dy:max(y)-dy1;
% [xq,yq]= meshgrid(x_want,y_want);
% Bsolint = griddata(x,y,vinv,xq,yq);
% Btarint = griddata(x,y,vtar,xq,yq);
end

message1 = fileoutinv;
  fid=fopen(message1,'w');
  fwrite(fid,Bsolint,'single');
  fclose(fid)
message2 = fileouttar;
  fid=fopen(message2,'w');
  fwrite(fid,Btarint,'single');
  fclose(fid)


