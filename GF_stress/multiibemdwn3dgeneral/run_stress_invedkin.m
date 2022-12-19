

%run the code
[res,para]=batch_runFromFiles



%stress tensor removing samples at the end
result.stc(1:5000,:,:,:)=res.stc(1:5000,:,:,:);
stc=result.stc;

%displacement field removing samples at the end
utc=res.utc(1:5000,:,:,:);
save('../out/stress.mat','stc','-v7.3');
%save('-ascii','../out/stress.ascii','sol')
save('../out/displacement.mat','utc','-v7.3');
%escalamiento de los resultados en tiempo
%s=load('../out/init_Gij3Dlayers.mat');

%clear sxx_x sxy_x sxz_x syy_x szz_x syz_x
%clear sxx_y sxy_y sxz_y syy_y szz_y syz_y
%clear sxx_z sxy_z sxz_z syy_z szz_z syz_z
%result=r;

%TO WRITE FILES AS NEEDED BY PREPROCESS
%CHANGE THE NEXT LINES ACCORDING TO YOUR
%GEOMETRY
sampt = size(result.stc);

%stations are the sources
%faults are the receivers
nsta=15;
nsub=400;

for j=1:nsub
%For S-x
for i=1:nsta;
   cont = i+(j-1)*nsta;
    sxx_x(:,cont) = result.stc(:,j,i,1);
    syy_x(:,cont) = result.stc(:,j,i,2);
    szz_x(:,cont) = result.stc(:,j,i,3);
    sxy_x(:,cont) = result.stc(:,j,i,4);
    sxz_x(:,cont) = result.stc(:,j,i,5);
    syz_x(:,cont) = result.stc(:,j,i,6);
end    
end

tx1 = ['../out/SIGMA_XX_C1']; % nombr del archivo
fileID = fopen(tx1,'w');
fwrite(fileID,sxx_x,'single');
fclose(fileID);
tx2 = ['../out/SIGMA_YY_C1']; % nombr del archivo
fileID = fopen(tx2,'w');
fwrite(fileID,syy_x,'single');
fclose(fileID);
tx3 = ['../out/SIGMA_ZZ_C1']; % nombr del archivo
fileID = fopen(tx3,'w');
fwrite(fileID,szz_x,'single');
fclose(fileID);
tx4 = ['../out/SIGMA_XY_C1']; % nombr del archivo
fileID = fopen(tx4,'w');
fwrite(fileID,sxy_x,'single');
fclose(fileID);
tx5 = ['../out/SIGMA_XZ_C1']; % nombr del archivo
fileID = fopen(tx5,'w');
fwrite(fileID,sxz_x,'single');
fclose(fileID);
tx6 = ['../out/SIGMA_YZ_C1']; % nombr del archivo
fileID = fopen(tx6,'w');
fwrite(fileID,syz_x,'single');
fclose(fileID);

%Jump to next stress components (force on y direction)
p=nsta;

for j=1:nsub
%For S-y
for i=1:nsta;
   cont = i+(j-1)*nsta;
    sxx_y(:,cont) = result.stc(:,j,i+p,1);
    syy_y(:,cont) = result.stc(:,j,i+p,2);
    szz_y(:,cont) = result.stc(:,j,i+p,3);
    sxy_y(:,cont) = result.stc(:,j,i+p,4);
    sxz_y(:,cont) = result.stc(:,j,i+p,5);
    syz_y(:,cont) = result.stc(:,j,i+p,6);
end
end

tx1 = ['../out/SIGMA_XX_C2']; % nombr del archivo
fileID = fopen(tx1,'w');
fwrite(fileID,sxx_y,'single');
fclose(fileID);
tx2 = ['../out/SIGMA_YY_C2']; % nombr del archivo
fileID = fopen(tx2,'w');
fwrite(fileID,syy_y,'single');
fclose(fileID);
tx3 = ['../out/SIGMA_ZZ_C2']; % nombr del archivo
fileID = fopen(tx3,'w');
fwrite(fileID,szz_y,'single');
fclose(fileID);
tx4 = ['../out/SIGMA_XY_C2']; % nombr del archivo
fileID = fopen(tx4,'w');
fwrite(fileID,sxy_y,'single');
fclose(fileID);
tx5 = ['../out/SIGMA_XZ_C2']; % nombr del archivo
fileID = fopen(tx5,'w');
fwrite(fileID,sxz_y,'single');
fclose(fileID);
tx6 = ['../out/SIGMA_YZ_C2']; % nombr del archivo
fileID = fopen(tx6,'w');
fwrite(fileID,syz_y,'single');
fclose(fileID);

%Jump to next stress components (force on z direction)
p=p+nsta;

for j=1:nsub
%For S-x
for i=1:nsta;
   cont = i+(j-1)*nsta;
    sxx_z(:,cont) = result.stc(:,j,i+p,1);
    syy_z(:,cont) = result.stc(:,j,i+p,2);
    szz_z(:,cont) = result.stc(:,j,i+p,3);
    sxy_z(:,cont) = result.stc(:,j,i+p,4);
    sxz_z(:,cont) = result.stc(:,j,i+p,5);
    syz_z(:,cont) = result.stc(:,j,i+p,6);
end
end

tx1 = ['../out/SIGMA_XX_C3']; % nombr del archivo
fileID = fopen(tx1,'w');
fwrite(fileID,sxx_z,'single');
fclose(fileID);
tx2 = ['../out/SIGMA_YY_C3']; % nombr del archivo
fileID = fopen(tx2,'w');
fwrite(fileID,syy_z,'single');
fclose(fileID);
tx3 = ['../out/SIGMA_ZZ_C3']; % nombr del archivo
fileID = fopen(tx3,'w');
fwrite(fileID,szz_z,'single');
fclose(fileID);
tx4 = ['../out/SIGMA_XY_C3']; % nombr del archivo
fileID = fopen(tx4,'w');
fwrite(fileID,sxy_z,'single');
fclose(fileID);
tx5 = ['../out/SIGMA_XZ_C3']; % nombr del archivo
fileID = fopen(tx5,'w');
fwrite(fileID,sxz_z,'single');
fclose(fileID);
tx6 = ['../out/SIGMA_YZ_C3']; % nombr del archivo
fileID = fopen(tx6,'w');
fwrite(fileID,syz_z,'single');
fclose(fileID);



