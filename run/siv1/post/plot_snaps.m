clear all
close all
clc

%option = input('grid option (1=regular) (2=iregular):')
opt = load('opt_grid.dat');
option = opt(1)

%time samples per snapshot
samp = [1 6 11 16 21 26 31 36];
%time samples for whole history
tsam=41; %for weights
nt=39;   %for source

if ( option == 1 )
%length along sides
 lx = 36;
 ly = 18;
else
%number of nodes
 ltot = 595;
end

%load files
 %rupture times used for contour lines
 trup=load('../dat/rupture_times.dat');
 %weighting matrix
 weig=load('../dat/weight_prior.ascii');
 %inverted sliprate
 if (option ==1 )
  tarsrc=load('../dat/model_target.dat');
 else
  tarsrc=load('../dat/model_target_triang.dat');
 end
 %target sliprate model
 invsrc=load('source.out');
 %prior model --- comment target and uncomment prior
 %tarsrc=load('../dat/prior_model.dat');

if ( option == 1 )
 %coordinates of nodes along strike and dip
 x=0.5:lx;
 y=0.5:ly;

 %interpolated coordinates
 xint=x(1):0.1:x(end);
 yint=y(1):0.1:y(end);

 %tables used to interpolate
 [X,Y]=meshgrid(x,y);
 [XI,YI]=meshgrid(xint,yint);
 %length of tables
 lyi=length(yint);
 lxi=length(xint);

 %rupture times used for contour lines
 trupr=reshape(trup,lx,ly)';
 %weighting matrix
 weisrc=reshape(weig,lx*ly,tsam)';

 %prepare snapshots and files
 for i=1:length(samp)
   %arrange arrays
   tar = reshape(tarsrc(samp(i),:),lx,ly)';
   inv = reshape(invsrc(samp(i),:),lx,ly)';
   w =   reshape(weisrc(samp(i),:),lx,ly)';
   %interpolate arrays
   tari = interp2(X,Y,tar,XI,YI);
   invi = interp2(X,Y,inv,XI,YI);
   wi = interp2(X,Y,w,XI,YI);
   %print maximum values of snapshots
   text1=sprintf('snapshot %01d', i);
   text2=sprintf('maxslipT maxslipI maxweight');
   display(text1)
   display(text2)
   maxvals =[max(max(tar)) max(max(inv)) max(max(w))]
   k = 1;
%  for j=1:lxi
%    p = lyi;
%    for jj=1:lyi
%     sr1(k,1:3) = [XI(1,j) YI(jj,1) invi(jj,j)];
%     sr2(k,1:3) = [XI(1,j) YI(jj,1) prii(jj,j)];
%     sr3(k,1:3) = [XI(1,j) YI(jj,1) wi(jj,j)];
%     k = k + 1;
%     p = p - 1;
%    end
%  end
  %save arrays in output files
   message1 = sprintf('file_inv_%d.bin',i);
   message2 = sprintf('file_pri_%d.bin',i);
   message3 = sprintf('file_wei_%d.bin',i);
   fid=fopen(message1,'w');
   fwrite(fid,tari,'single');
   fclose(fid);
   fid=fopen(message2,'w');
   fwrite(fid,invi,'single');
   fclose(fid);
   fid=fopen(message3,'w');
   fwrite(fid,wi,'single');
   fclose(fid);
   k = k + 1;
 end
 %save contour lines of estimated rupture times
 fid=fopen('contour.bin','w');
 fwrite(fid,trupr,'single');
 fclose(fid);
else
 %read grid points
 grid = load('../dat/grid.in');
 dx = 100; dy=100;
 dx1= 300; dy1=300;
 x = grid(:,1);
 y = grid(:,2);
 %prepare snapshots and files
 for i=1:length(samp)
  vinv = reshape(invsrc(samp(i),:),ltot,1);
  vtar = reshape(tarsrc(samp(i),:),ltot,1);
  %points where to interpolate
  x_want = min(x)+dx1:dx:max(x)-dx1;
  y_want = min(y)+dy1:dy:max(y)-dy1;
  [xq,yq]= meshgrid(x_want,y_want);
  invi = griddata(x,y,vinv,xq,yq);
  tari = griddata(x,y,vtar,xq,yq);
   %save arrays in output files
    message1 = sprintf('file_inv_%d.bin',i);
    message2 = sprintf('file_pri_%d.bin',i);
    message3 = sprintf('file_wei_%d.bin',i);
    fid=fopen(message1,'w');
    fwrite(fid,tari,'single');
    fclose(fid);
    fid=fopen(message2,'w');
    fwrite(fid,invi,'single');
    fclose(fid);
    fid=fopen(message3,'w');
    fwrite(fid,invi,'single');
    fclose(fid);
 end
%save contour lines of estimated rupture times
%fid=fopen('contour.bin','w');
%fwrite(fid,trupr,'single');
%fclose(fid);
end


%hypocenter location along strike and dip (km)
hyp = [27; 12];

%save contour lines of estimated rupture times
fid=fopen('epi.bin','w');
fwrite(fid,hyp,'single');
fclose(fid);

