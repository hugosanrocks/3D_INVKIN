function [RESULT,para] = batch_runFromFiles
%batch_runFromFiles Resuelve el modelo definido en los archivos de datos.
%   Resuelve la función de Green de esfuerzos y desplazamientos en un medio
%   estratificado dadas fuentes y receptores definidos en los archivos:
%
% 	../out/init_Gij3Dlayers.mat   inicializar variables
% 	../ins/vmodelV2.dat           estructura de velocidades
% 	../ins/sta.dat                donde se requiere aplicar las fuerzas
% 	../ins/fault.dat            	donde se requiere grabar los tensores
%
%   Los resutlados de esta función son la variable de parámetros 'para' y
%   la variable de resultados 'RESULT'
%
%% Inicializar variables
clear para
thisDir = pwd; disp(['at: ' thisDir])
load('../out/init.mat')
para.zeropad = 2^13; %hugo pone esto
%para.pulso.tipo = 5;
%para.pulso.a = 0.25;
%paral.pulso.b = 1;
%para.nf = 2048*2;
%para.fmax = 14;
%para.Nbptkx = 21000;

%rt =0.5
para.nf  =2048;
para.fmax = 8;

para.DWNxl = 1000; % espacio entre fuentes virtuales ( L )

% ajustar variables dependientes del sistema
      [pathstr1,pathstr2,~] = fileparts(pwd);
      [~,name,ext] = fileparts(para.name);
      para.name = [pathstr1,pathstr1(1),pathstr2,pathstr1(1),name,ext];
      clear ext
      para.nametmp = para.name;
      para.nomcarpeta = pwd;
      para.nomrep = [pathstr1,pathstr1(1),'out'];

disp('1) Inicializar variables ... done')
% Usar modelo de velocidades del archivo
indat = importdata('../ins/vmodelV2.dat');
para.nsubmed = size(indat.data,1); % cantidad de estratos + semiespacio
% cargar valores 
para.reg(1).sub = [];
for i=1:para.nsubmed
  para.reg(1).sub(i).rho      = indat.data(i,4);
  para.reg(1).sub(i).alpha    = indat.data(i,2);
  para.reg(1).sub(i).bet      = indat.data(i,3);
  para.reg(1).sub(i).qd       = 10000;  % factor de calidad :: elástico
  para.reg(1).sub(i).tipoatts = 1; % tipo de modelo de atenuación :: Q
  para.reg(1).sub(i).lambda   = para.reg(1).sub(i).rho*(para.reg(1).sub(i).alpha^2-2*para.reg(1).sub(i).bet^2);
  para.reg(1).sub(i).mu       = para.reg(1).sub(i).rho*para.reg(1).sub(i).bet^2;
  para.reg(1).sub(i).h        = indat.data(i,1);
end
para.reg(1).sub(i).h        = 0; %ultimo estrato=semi espacio

disp('2) Usar modelo de velocidades del archivo ...done')
% Usar receptores indicados en archivo
indat = importdata('../ins/fault.dat');
nFault = size(indat,1); %130
para.recpos=3; %receptores en posicion libre
para.chgrec = 0;
para.rec.nrecx= size(indat,1); % cantidad de receptores
para.rec.xr= indat(:,1);
para.rec.yr= indat(:,2);
para.rec.zr= indat(:,3);
disp('3) Usar receptores indicados en archivo ...done')

% Usar fuentes indicadas en archivo
indat = importdata('../ins/sta.dat'); 
n = size(indat,1); %40
para.ninc = n*3; % cantida de fuentes * 3 direcciones
para.fuente = 2; % son Fuentes Puntuales

% posiciones:
para.xs=[indat(:,1);indat(:,1);indat(:,1)];
para.ys=[indat(:,2);indat(:,2);indat(:,2)];
para.zs=[indat(:,3);indat(:,3);indat(:,3)];

% dirección +x
para.gam(1:n) = 90;
para.phi(1:n) = 0;
% dirección +y
para.gam(n+1:2*n) = 90;
para.phi(n+1:2*n) = 90;
% dirección +z
para.gam(2*n+1:3*n) = 180;
para.phi(2*n+1:3*n) = 0;

clear indat i
disp('4) Usar fuentes indicadas en archivo ...done')

% dibujar configuracion geometrica

bouton = [];
bouton.iinc = 2;
figure(1001)
set(gcf,'Name','Configuracion geométrica','numberTitle','off');
dibujo_conf_geo(para,gca)
%% Ejecución del análisis
para.espyinv=1;
[RESULT,para]=calculo(para);
% los tensores están en
%    RESULT.sw(:,:,:,:)
%              | | | '--- componente: sxx syy szz sxy sxz syz
%              | | '----- fuente: [1:n dirx  1:n diry  1:n dirz] 
%              | '------- receptor: 1:nFault
%              '--------- datos en frecuenica positiva 1:(para.nf/2+1) 

%% Escritura de resultados
disp('5) Escribiendo archivos binarios')
% Archivos independientes para cada componente del tensor de esfuerzos y dirección
txDir = ['x','y','z'];
txComp={'sxx', 'syy', 'szz', 'sxy', 'sxz', 'syz'};
nFrec = size(RESULT.sw,1);

iniSta = 1;
for iDir = 1:3 % cada direccion de la fuerza
finSta = n*iDir;
for iComp = 1:6 % cada componente
%  un renglón por cada receptor*fuente  :: iv
%  una columna por cada frecuencia
%           5200  , 201 frecuencias
v = zeros(nFault*n,nFrec);
iv = 1;

for iFault = 1:nFault %receptores
for iSta = iniSta:finSta %fuentes
%[record1fsta1,rfault1,...,record40fsta40,rfault1,record41fsta1,rfault2,...,record80fsta40,rfault2,...
%...,record5161fsta1,rfault130,...,record5200fsta40,rfault130]
v(iv,:) = RESULT.sw(:,iFault,iSta,iComp);
iv = iv+1;
end %iSta
end %iFault
% guardar archivo en formato binario
tx = ['../out/' txComp{iComp} '-' txDir(iDir) '.bin']; % nombr del archivo
%COMMENTED NOT TO WRITE STRESS ON FREQUENCY DOMAIN
%fileID = fopen(tx,'w');
v_real = real(v);
v_imag = imag(v);
adjacent = [v_real v_imag];
%DO NOT WRITE IN FREQUENCY DOMAIN
%fwrite(fileID,adjacent,'double');
%fclose(fileID);

disp (['File ' tx ' variable size:']);
disp (size(adjacent));

end %iComp
iniSta= (n*iDir)+1;
disp (['done saving tensors for direction ' num2str(iDir)])

end %iDir

% archivo de información
fileID = fopen('../out/info.txt','w');
fprintf(fileID, 'Number of frequencies:\n');
fprintf(fileID, [num2str(nFrec) '\n']);
fprintf(fileID, 'Number of rows:\n');
fprintf(fileID, [num2str(size(adjacent,1)) '\n']);
fprintf(fileID, 'Number of columns:\n');
fprintf(fileID, [num2str(size(adjacent,2)) '\n']);
fprintf(fileID, 'Number of Fault points(receibers):\n');
fprintf(fileID, [num2str(nFault) '\n']);
fprintf(fileID, 'Number of Stations (sources):\n');
fprintf(fileID, [num2str(n) '\n']);
fclose(fileID);
disp('done')
end

