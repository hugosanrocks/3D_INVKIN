%% importar datos
clear
cd ..
cd ins
rawdata = load('MAPA.DAT');
cd ..
cd multi-dwn-ibem.matlab
% quitar los que están fuera de la zona lacustre
%del = [540 547 550 551 552 553 570 593 595];
del = [541 547 555 570 593 595 616 617 618 630 639 640 641 642 645 659 675];
% del = [];
rawdata(del,:) = [];
ndel = length(del);
data = rawdata(1:1010-ndel,:);      % lat lon periodo
borde = rawdata(1011-ndel:2516-ndel,:); % lat lon periodo
allda = rawdata(1:2516-ndel,:);

figure('Color','white'); hold on
axis on; grid on;
plot(borde(:,1),borde(:,2),'r.')
% plot(data(:,1),data(:,2),'b.')
plot3(data(:,1),data(:,2),data(:,3),'k.')
%
f = scatteredInterpolant(allda(:,1:2),allda(:,3));
%
n = 90;
xlin = linspace(min(rawdata(:,1)),max(rawdata(:,1)),n);
ylin = linspace(min(rawdata(:,2)),max(rawdata(:,2)),n);
[X,Y] = meshgrid(xlin,ylin);
Z = f(X,Y);
% figure
% mesh(X,Y,Z) %interpolated
axis tight; hold on
contour(X,Y,Z,50)
% plot3(X,Y,Z,'.','MarkerSize',15) %nonuniform


%%
%contour
%deg2km
%%
clear
close all
cd /Users/marshall/Documents/DOC/coco/multi-dwn-ibem.matlab
cd ..
cd ins
I = imread('profs.png');
I = flipud(I);
rawdata = load('MAPA.DAT');
cd ..
cd multi-dwn-ibem.matlab


% tamaño de imagen en pixeles
HorDimPx = size(I,2);
VerDimPx = size(I,1);
% WorldLimits
pasitoD = 10/60/16;
pasito = deg2km(pasitoD); % distancia en km entre crucecitas
escala = 1;
% en km
% xWorldLimits = [deg2km(-(99+10/60))-6.2*pasito deg2km(-(99))+15.4*pasito];%[0 HorDimPx*escala];%
% yWorldLimits = [deg2km(19+20/60)-13.8*pasito deg2km(19+30/60)+9.2*pasito];%[0 VerDimPx*escala];%
% en grados
xWorldLimits = [(-(99+10/60))-6.2*pasitoD (-(99))+15.4*pasitoD];%[0 HorDimPx*escala];%
yWorldLimits = [(19+20/60)-13.85*pasitoD (19+30/60)+9.2*pasitoD];%[0 VerDimPx*escala];%
R = imref2d(size(I),xWorldLimits,yWorldLimits) ;

% reducir ruido
% h = ones(6,6)/36; % <--- parametro que se ha calibrado
% h = ones(8,8)/64; % <--- parametro que se ha calibrado
% I = imfilter(I,h);

% convertirla en una imagen binaria a partir de un umbral de intensidad.
umbral = 0.4; 
BW = im2bw(I, umbral); % <--- parametro que se ha calibrado

% extraer contornos de la imagen
figure('Name','Datos de curvas de nivel'); hold on
% identificar las fronteras:
[B,L,N] = bwboundaries(BW,4,'noholes'); % en coordenadas de pixeles
imshow(BW,R,'InitialMagnification','fit');
hold on
% del analisis de los contornos elegimos solo as curvas de nivel:
clear ik
i = 1; ik(i) = 3 ; la(i) = 10; ii(i) = 1; ff(i) = 9161;
i = 2; ik(i) = 3 ; la(i) = 20; ii(i) = 9287; ff(i) = 15662;
i = 3; ik(i) = 3 ; la(i) = 10; ii(i) = 15723; ff(i) = length(B{ik(end),1}(:,2));
i = 4; ik(i) = 6 ; la(i) = 30; ii(i) = 1; ff(i) = length(B{ik(end),1}(:,2));
i = 5; ik(i) = 7 ; la(i) = 40; ii(i) = 1; ff(i) = length(B{ik(end),1}(:,2));
i = 6; ik(i) = 8 ; la(i) = 20; ii(i) = 1; ff(i) = 1677;
i = 7; ik(i) = 8 ; la(i) = 30; ii(i) = 1678; ff(i) = 5453;
i = 8; ik(i) = 8 ; la(i) = 20; ii(i) = 5454; ff(i) = length(B{ik(end),1}(:,2));
i = 9; ik(i) = 11 ; la(i) = 40; ii(i) = 1; ff(i) = length(B{ik(end),1}(:,2));
i = 10;ik(i) = 13 ; la(i) = 10; ii(i) = 1; ff(i) = 545;
i = 11;ik(i) = 13 ; la(i) = 10; ii(i) = 1278; ff(i) = length(B{ik(end),1}(:,2));
i = 12;ik(i) = 15 ; la(i) = 50; ii(i) = 1; ff(i) = length(B{ik(end),1}(:,2));
i = 13;ik(i) = 22 ; la(i) = 30; ii(i) = 1; ff(i) = 161;
i = 14;ik(i) = 22 ; la(i) = 20; ii(i) = 162; ff(i) = 568;
i = 15;ik(i) = 22 ; la(i) = 30; ii(i) = 569; ff(i) = length(B{ik(end),1}(:,2));
i = 16;ik(i) = 25 ; la(i) = 10; ii(i) = 1; ff(i) = 100;
i = 17;ik(i) = 25 ; la(i) = 10; ii(i) = 318; ff(i) = length(B{ik(end),1}(:,2));
i = 18;ik(i) = 40 ; la(i) = 20; ii(i) = 1; ff(i) = length(B{ik(end),1}(:,2));
i = 19;ik(i) = 41 ; la(i) = 10; ii(i) = 1; ff(i) = 122;
i = 20;ik(i) = 41 ; la(i) = 10; ii(i) = 482; ff(i) = length(B{ik(end),1}(:,2));
i = 21;ik(i) = 52 ; la(i) = 20; ii(i) = 1; ff(i) = length(B{ik(end),1}(:,2));
i = 22;ik(i) = 56 ; la(i) = 10; ii(i) = 1; ff(i) = 469;
i = 23;ik(i) = 56 ; la(i) = 10; ii(i) = 813; ff(i) = length(B{ik(end),1}(:,2));
i = 24;ik(i) = 57 ; la(i) = 10; ii(i) = 1; ff(i) = length(B{ik(end),1}(:,2));
i = 25;ik(i) = 61 ; la(i) = 20; ii(i) = 1; ff(i) = length(B{ik(end),1}(:,2));
i = 26;ik(i) = 66 ; la(i) = 20; ii(i) = 1; ff(i) = length(B{ik(end),1}(:,2));
i = 27;ik(i) = 67 ; la(i) = 20; ii(i) = 1; ff(i) = length(B{ik(end),1}(:,2));

% RR para referenciar al sistema coordenado del mapa
RR = makerefmat(R.YWorldLimits(1), R.XWorldLimits(1), ...
(R.YWorldLimits(2) - R.YWorldLimits(1))/(R.ImageSize(1)),...
(R.XWorldLimits(2) - R.XWorldLimits(1))/(R.ImageSize(2))...
);

% trazar curvas de nivel
col = linspace(0.5,0.9,50);
latS = []; lonS = []; proF = [];
for k = 1:length(ik)
   [lat,lon] = pix2latlon(RR,B{ik(k),1}(ii(k):ff(k),2).*escala,...
                             B{ik(k),1}(ii(k):ff(k),1).*escala);
   latS = [latS;lat];
   lonS = [lonS;lon];
   proF = [proF;lat.*0 + la(k)];
   plot(lat,lon,'--','Color',col(la(k)).*[0 0 1], 'LineWidth', 0.1,'DisplayName',num2str(la(k)));
end

% trazar el borde de la cuenca de la otra serie de datos
borde = rawdata(1011:2516,1:2); % lat lon borde de la cuenca
plot(borde(:,1),borde(:,2),'r.')
set(gca,'yDir', 'normal')
latS = [latS;borde(:,1)];
lonS = [lonS;borde(:,2)];
proF = [proF;borde(:,1).*0];

% bordes extra para mejorar la interpolación
latS = [latS;linspace(min(borde(:,1)),max(borde(:,1)),10)';...
             linspace(min(borde(:,1)),min(borde(:,1)),10)';...
             linspace(max(borde(:,1)),min(borde(:,1)),10)';...
             linspace(max(borde(:,1)),max(borde(:,1)),10)'];
lonS = [lonS;linspace(min(borde(:,2)),min(borde(:,2)),10)';...
             linspace(min(borde(:,2)),max(borde(:,2)),10)';...
             linspace(max(borde(:,2)),max(borde(:,2)),10)';...
             linspace(max(borde(:,2)),min(borde(:,2)),10)'];
proF = [proF;linspace(min(borde(:,1)),max(borde(:,1)),10)'.*0;...
             linspace(min(borde(:,1)),min(borde(:,1)),10)'.*0;...
             linspace(max(borde(:,1)),min(borde(:,1)),10)'.*0;...
             linspace(max(borde(:,1)),max(borde(:,1)),10)'.*0];

latS = [latS;[-99.03;-99.03;-99.07;-99.05;-99.09;-99.05;-99.032;-99.01;-98.99]];
lonS = [lonS;[19.45 ;19.43 ;19.41 ;19.42 ;19.28 ;19.28 ;19.275 ;19.27; 19.26 ]];
proF = [proF;[60;60;60;60;45;40;40;40;40]];


%% interpolar una malla de datos
n = 40; %30 % El tamaño de la malla

f = scatteredInterpolant(latS,lonS,proF);
xlin = linspace(min(rawdata(:,1)),max(rawdata(:,1)),n);
ylin = linspace(min(rawdata(:,2)),max(rawdata(:,2)),n);
[X,Y] = meshgrid(xlin,ylin);
Z = f(X,Y);
figure('Name','Curvas de nivel interpoladas');
axis tight; hold on
contour(X,Y,Z,7)
plot(borde(:,1),borde(:,2),'r.')

figure('Name','Malla interpolada');
axis tight; hold on
h = mesh(X,Y,Z); %interpolated
plot(borde(:,1),borde(:,2),'r.')
%
bol = h.VertexNormals(:,:,3)<0.99 * h.ZData(:,:) > 0;
XData = h.XData(bol);%XData = [h.XData(bol);borde(:,1)];
YData = h.YData(bol);%YData = [h.YData(bol);borde(:,2)];
ZData = h.ZData(bol);%ZData = [h.ZData(bol);borde(:,1).*0];
figure('Name','Vértices seleccionados');
plot3(XData,YData,ZData,'.')
hold on; plot(borde(:,1),borde(:,2),'r.')

%% Triangular
DT = delaunayTriangulation([XData,YData]); %lista de conectividad
DT3 = delaunayTriangulation([XData,YData,ZData]);
P = DT3.Points;
TR = triangulation(DT.ConnectivityList,P);
CN = incenter(TR);
inside = true(length(CN),1); 
area = zeros(length(CN),1);
for i=1:length(CN)
  t = TR.Points(TR.ConnectivityList(i,:),:);
  % lados del tríangulo al cuadrado
  a = (t(1,1:2) - t(2,1:2)).^2; a2 = a(1)+a(2);% v1 a v2
  b = (t(2,1:2) - t(3,1:2)).^2; b2 = b(1)+b(2); % v2 a v3
  c = (t(3,1:2) - t(1,1:2)).^2; c2 = c(1)+c(2); % v3 a v1 
sum = a2+b2+c2;
area(i) = 0.25 * sqrt(4*a2*b2 - (a2+b2-c2)^2);
dis = 0.2;
if (a2 < dis*sum) || (b2 < dis*sum) || (c2 < dis*sum)
  inside(i) = false;
end
end
aProm = mean(area);
inside(area>2*aProm) = false;

T = DT.ConnectivityList(inside,:);

% Ajustar el borde de la malla a la superficie libre
TR = triangulation(T,P);
E = freeBoundary(TR);
P(E(:,1),3) = 0;

% la malla final
TR = triangulation(T,P);
FN = faceNormal(TR);
% CN = incenter(TR);

figure('Name',['Malla de ' num2str(length(FN)) ' triángulos']);
trimesh(TR);
hold on; plot(borde(:,1),borde(:,2),'r.')
% plot(TR.Points(E(:,1),1),TR.Points(E(:,1),2),'g.') %
view([0,90])


% grados a metros
centro(1) = mean(P(:,1)); centro(2) = mean(P(:,2)); centro(3) = 0;
Pc(:,1) = deg2km(P(:,1) - centro(1))*1000;
Pc(:,2) = deg2km(P(:,2) - centro(2))*1000;
Pc(:,3) = P(:,3); 
% TR = triangulation(T,Pc);
% trimesh(TR);

%% escribir STL de frontera de continuidad
clear t
t = sprintf('%s \n','solid ValleMex');
for i = 1:length(FN)
t = sprintf('%s%s %f %f %f \n',t,'facet normal ',FN(i,1),FN(i,2),FN(i,3));
t = sprintf('%s%s \n',t,'outer loop');
t = sprintf('%s%s %f %f %f \n',t,'vertex ',Pc(T(i,1),1),Pc(T(i,1),2),Pc(T(i,1),3));
t = sprintf('%s%s %f %f %f \n',t,'vertex ',Pc(T(i,2),1),Pc(T(i,2),2),Pc(T(i,2),3));
t = sprintf('%s%s %f %f %f \n',t,'vertex ',Pc(T(i,3),1),Pc(T(i,3),2),Pc(T(i,3),3));
t = sprintf('%s%s \n',t,'endloop');
t = sprintf('%s%s \n',t,'endfacet');
end
t = sprintf('%s%s',t,'endsolid ValleMex');
t = sprintf('%s%s \n',t,'endfacet');
cd ..
cd out
fileID = fopen('VallemexCont.stl','w');
fprintf(fileID,'%s',t);
fclose(fileID);
cd ..
cd multi-dwn-ibem.matlab

%% croissant
% Sánchez-Sesma & Luzon 1995
% a > b  && b = 0.7a && h = 0.4/a && a = 4km
% f(x,y) = h(b^2 - P^2)[1 - 2a(a-x)/P^2]
Bmin = 1; %km/s
fmax = 1; %Hertz
lamdamin = Bmin/fmax/6; %km
numelems = 4/lamdamin/2; %cantidad de triangulos 

n = 37; %(número non) %35
xlin = linspace(-5,5,n);
ylin = linspace(-5,5,n);
[x,y] = meshgrid(xlin,ylin);
z = x.*0;
for i=1:length(xlin)
for j=1:length(ylin)
z(i,j) = Croissant(x(i,j),y(i,j));
end
end
fi = figure('Name','Curvas de nivel Croissant');
axis tight; hold on
h = mesh(x,y,z);
bol = h.VertexNormals(:,:,3)~=1;
XData = h.XData(bol);%XData = [h.XData(bol);borde(:,1)];
YData = h.YData(bol);%YData = [h.YData(bol);borde(:,2)];
ZData = h.ZData(bol);%ZData = [h.ZData(bol);borde(:,1).*0];
DT = delaunayTriangulation([XData,YData]); disp('lista de conectividad')
DT3 = delaunayTriangulation([XData,YData,ZData]);
Pc = DT3.Points;
TR = triangulation(DT.ConnectivityList,Pc);
E = freeBoundary(TR);
Pc(E(:,1),3) = 0;
FN = faceNormal(TR); disp(['Numero de caras del cuernito: ', num2str(length(FN))])
bol = FN(:,3) ~= 1;
T = TR.ConnectivityList(bol,:);
TR = triangulation(T,Pc);
FN = faceNormal(TR);
Cen = zeros(length(FN),1);
for i = 1:length(FN)
    Cen(i) = mean(Pc(T(i,1:3),2));
end
% hacerlo simétrico en y
iflp = Cen(:,1)<0; % espejear
NN = sum(iflp)*2;
PcT = zeros(NN,3,3); % i,vertice,componente
i = 1;
norm = zeros(NN,3);
for ind = 1:length(iflp)
    if iflp(ind) 
        % copiar original
norm(i,:) = FN(ind,:);
PcT(i,1,1)=Pc(T(ind,1),1); PcT(i,1,2)=Pc(T(ind,1),2); PcT(i,1,3)=Pc(T(ind,1),3);
PcT(i,2,1)=Pc(T(ind,2),1); PcT(i,2,2)=Pc(T(ind,2),2); PcT(i,2,3)=Pc(T(ind,2),3);
PcT(i,3,1)=Pc(T(ind,3),1); PcT(i,3,2)=Pc(T(ind,3),2); PcT(i,3,3)=Pc(T(ind,3),3);  
        % espejear en y
norm(i+1,1:2:3) =  norm(i,1:2:3); 
norm(i+1,  2  ) = -norm(i,  2  );
PcT (i+1,3:-1:1,1:2:3) =  PcT(i,1:3,1:2:3); 
PcT (i+1,3:-1:1,  2  ) = -PcT(i,1:3,  2  );
i = i+2;
    end
end
Cen = zeros(NN,3);
for i = 1:NN
    Cen(i,1) = mean(PcT(i,1:3,1));
    Cen(i,2) = mean(PcT(i,1:3,2));
    Cen(i,3) = mean(PcT(i,1:3,3));
end
%
% figure('Name',['Malla de ' num2str(length(FN)) ' triángulos']);
% % trimesh(TR);
% hold on
% for i = 1:NN
%     patch(PcT(i,1:3,1),PcT(i,1:3,2),PcT(i,1:3,3),'blue');
%     quiver3(Cen(i,1), Cen(i,2), Cen(i,3),...
%            norm(i,1),norm(i,2),norm(i,3),0.15,'Color','r');
% end
% axis equal
% %
clear nam
nam = ['Cuerno_',num2str(n),'_A'];

clear t
t = sprintf('%s \n',['solid ',nam]);
for i = 1:NN % length(FN)
t = sprintf('%s%s %f %f %f \n',t,'facet normal ',norm(i,1),norm(i,2),norm(i,3)); % FC
t = sprintf('%s%s \n',t,'outer loop');
% t = sprintf('%s%s %f %f %f \n',t,'vertex ',Pc(T(i,1),1),Pc(T(i,1),2),Pc(T(i,1),3));
% t = sprintf('%s%s %f %f %f \n',t,'vertex ',Pc(T(i,2),1),Pc(T(i,2),2),Pc(T(i,2),3));
% t = sprintf('%s%s %f %f %f \n',t,'vertex ',Pc(T(i,3),1),Pc(T(i,3),2),Pc(T(i,3),3));
t = sprintf('%s%s %f %f %f \n',t,'vertex ',PcT(i,1,1),PcT(i,1,2),PcT(i,1,3));
t = sprintf('%s%s %f %f %f \n',t,'vertex ',PcT(i,2,1),PcT(i,2,2),PcT(i,2,3));
t = sprintf('%s%s %f %f %f \n',t,'vertex ',PcT(i,3,1),PcT(i,3,2),PcT(i,3,3));
t = sprintf('%s%s \n',t,'endloop');
t = sprintf('%s%s \n',t,'endfacet');
end
t = sprintf('%s%s',t,['endsolid ',nam]);
t = sprintf('%s%s \n',t,'endfacet');
cd ..
cd out
fileID = fopen([nam,'.stl'],'w');
fprintf(fileID,'%s',t);
fclose(fileID);
cd ..
cd multi-dwn-ibem.matlab

clear nam
nam = ['Cuerno_',num2str(n),'_B'];
clear t
t = sprintf('%s \n',['solid ',nam]);
for i = 1:NN% length(FN)
t = sprintf('%s%s %f %f %f \n',t,'facet normal ',0,0,-1);
t = sprintf('%s%s \n',t,'outer loop');
% t = sprintf('%s%s %f %f %f \n',t,'vertex ',Pc(T(i,3),1),Pc(T(i,3),2),0);
% t = sprintf('%s%s %f %f %f \n',t,'vertex ',Pc(T(i,2),1),Pc(T(i,2),2),0);
% t = sprintf('%s%s %f %f %f \n',t,'vertex ',Pc(T(i,1),1),Pc(T(i,1),2),0);
t = sprintf('%s%s %f %f %f \n',t,'vertex ',PcT(i,3,1),PcT(i,3,2),0);
t = sprintf('%s%s %f %f %f \n',t,'vertex ',PcT(i,2,1),PcT(i,2,2),0);
t = sprintf('%s%s %f %f %f \n',t,'vertex ',PcT(i,1,1),PcT(i,1,2),0);
t = sprintf('%s%s \n',t,'endloop');
t = sprintf('%s%s \n',t,'endfacet');
end
t = sprintf('%s%s',t,['endsolid ',nam]);
t = sprintf('%s%s \n',t,'endfacet');
cd ..
cd out
fileID = fopen([nam,'.stl'],'w');
fprintf(fileID,'%s',t);
fclose(fileID);
cd ..
cd multi-dwn-ibem.matlab
disp('end')

