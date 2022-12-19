function [phi_fv,coord,DWN,para]=sol3D_dfv(para,fj,DWN,J)
% calcula la densidad de fuerzas virtuales (dfv)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% discretizacion de los contornos %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if para.dim == 3 %  3D axisimétrico
  coord   = malla_geom(para,fj);
  coord   = malla_geom_axi(coord,para);
elseif para.dim == 4% 3D goemetría general
  para = cargardelista(para,J);
  para = malla_geom3Dgen(para,fj);
  coord = para.coord;
end

if para.GraficarCadaDiscretizacion % a cada frecuencia
  figure(100);cla
  dibujo_conf_geo(para,gca)
  hold on
  if para.dim == 3
    for i_elem = 1:coord.nbpt % cada disco
      plotCircle3D([coord.x(i_elem) coord.y(i_elem) coord.z(i_elem)],...
        [coord.vnx(i_elem) coord.vny(i_elem) coord.vnz(i_elem)],...
        sqrt(coord.dA(i_elem)/pi),2);
    end
  elseif para.dim == 4
    for m=2:para.nmed
      for p = 1:para.cont(m,1).NumPieces % para cada pieza del contorno
        kind = para.cont(m,1).piece{p}.kind;
        if kind ~= 3 % si no es una frontera auxiliar
          if size(para.cont(m,1).piece{p}.subdibData,2)>0 % y hay datos
            % pintar todos los triángulos subdivididos
            sdD = para.cont(m,1).piece{p}.subdibData;
            for t = 1:size(sdD.triangles,3)
              plot3(sdD.triangles(1,[1 2 3 1],t),...
                sdD.triangles(2,[1 2 3 1],t),...
                sdD.triangles(3,[1 2 3 1],t),'w--','LineWidth',1.0)
            end
          end
        end
      end
    end
  end
  hold off
  drawnow
  % para guardar la imagen:
  %   [~,~]=mkdir('../out/subdib');
  %   cd ../out/subdib
  %   AxesH = gca;
  %   hf = figure('Visible','on','NumberTitle','off',...
  %       'Name',['Discretization at fj ' num2str(fj)]);
  %   copyobj(AxesH,hf);
  %   savefig(hf,['Disc_' num2str(fj) '.fig'],'version2'); %MATLAB 2014b or higher!
  %   close(hf)
  %   cd ../../multi-dwn-ibem.matlab
end

% % % Para cuando nos queremos saltar el ibem:
% % nbeq    = coord.nbeq;
% % nbeq2   = 3*nbeq;
% % phi_fv =zeros(nbeq2,para.ninc);
% % warning('debug skip');
% % return

%el dr varia en funcion de la velocidad y el dA en funcion de v^2
nbeq    = coord.nbeq;
nbpt    = coord.nbpt;
ju      = coord.ju;
nbeq2   = 3*nbeq;       %hay ecuaciones para phix, phiy y phiz
% disp(['ifreq = ',num2str(para.j),' f=',num2str(fj),', system size: ',num2str(nbeq2),'^2'])
% str = sprintf([ '                               system size: ',num2str(nbeq2),'^2']);
% disp(str);

% phi_fv =zeros(nbeq2,para.ninc);
% warning('l 66')
% return

% lStr = length(str);
% lPrompt = 10;
%% IBEM-DWN
if para.geo(1)==3
  %%%%%%%%%%%%
  % IBEM-DWN %
  %%%%%%%%%%%%
  % en esta parte, se calcula la matriz del DWN, una primera vez con paso
  % constante, otra con paso que permite mejor discretizar los polos.
  % Csc: acceleracion de los calculos del DWN para rellenar la matriz del
  % IBEM.
  % Y tambien se busca los puntos en los cuales hay que calcular G^DWN,
  % tratando conjuntamente todos los receptores que tienen una misma z,
  
  %se propone a la primera iteracion un vector de la componente
  %horizontal de kr con paso constante lo que corresponde a
  %xmax = para.DX*para.DWNnbptkx/2;
  %los puntos en polar estan en el centro del segmento de longitud DK
  nk          = para.DWNnbptkx;
  DK      	= para.DWNkmax/nk;%kmax posicion central del ultimo
  %xmax doit etre > vmax*tmax=alpha/dfe
  kr          = (0.5+(0:nk))*DK;
  DWN.kr      = kr;
  DWN.dkr     = kr*0+DK;
  DWN.k2      = kr;%copia para rebuildk2, sino inutil
  DWN.omegac  = 2*pi*fj;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  % calculo DWN - parte 1 %
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  % calculo de la matriz del DWN
  DWN         = calcul_A_DWN_3D_polar_Ncapas_HS(para,DWN);
  %     DWN         = rebuildk2(DWN);
  %     DWN.kr      = DWN.k2;
  %     DWN.dkr     = DWN.dk2;
  %     DWN.dkr(1)  = DWN.kr(2)-DWN.dkr(2)/2;
  %     DWN         = calcul_A_DWN_3D_polar_Ncapas_HS(para,DWN);
  % inversion de las matrices
  
  notInvA = DWN.A_DWN;
  notInvB = DWN.B_DWN;
  parfor i=1:length(DWN.kr)
    %       disp(['****', num2str(i), '****']);
    invA(:,:,i)=inv(notInvA(:,:,i));
    %       disp('A');
    invB(:,:,i)=inv(notInvB(:,:,i));
    %       disp('B'); disp('');
  end
  DWN.A_DWN = invA;
  DWN.B_DWN = invB;
  clear notInvA notInvB invA invB
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % identificacion de los puntos en los cuales se tiene que calcular DWN -parte 2 %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% para el calculo de la matriz del IBEM:
  %   posicion puntos de colocacion: receptores y fuentes (virtuales)
  
  indpt               = 1:coord.nbpt;         %indice de todos los puntos
  indi                = coord.indm(1,1).ind;  %logical de los puntos perteneciendo a m=1
  indir               = indpt(indi);          %indice de estos puntos
  
  xrv                 = coord.x(indi);
  yrv                 = coord.y(indi);
  zrv                 = coord.z(indi);
  nxrv                = length(xrv);
  vn(1,1:nxrv)        = coord.vnx(indi);
  vn(2,1:nxrv)        = coord.vny(indi);
  vn(3,1:nxrv)        = coord.vnz(indi);
  salu                = ones(nxrv,1);
  sals                = ones(nxrv,1);
  
  %reduccion du numero de profundidades de las funtes virtuales, moviendo
  %los puntos en z de manera a agrupar puntos cercanos (<lambda/10)
  [zrfv,izrfv,zrref]	= pseudo_unique(zrv,para);
  % zrfv=zrv.';zrref=zrv;izrfv=(1:length(zrv)).';
  coord.z(indi)       = zrref; %actualizacion de la posicion modificada de los puntos de colocacion
  
  %%% se incluye de una ves los receptores reales para el calculo de los
  %   campos difractados.
  %   El campo difractado es indepediente de las
  %   amplitudes de las fuerzas en cada punto de colocacion.
  %   Eso se toma en cuenta  hasta tomar en cuenta los phi cf inversion_3D_k
  xrr                 = para.rec.xr(para.rec.m(1).ind).';
  yrr                 = para.rec.yr(para.rec.m(1).ind).';
  zrr                 = para.rec.zr(para.rec.m(1).ind).';
  nrr                 = para.rec.m(1).nr;
  
  [zrfr,izr0,~]       = pseudo_unique(zrr,para);
  %         zrr=zrr.';izr0=(1:length(zrr)).';
  
  % inicialisacion del vector de los diffractados en los receptores
  % reales y salidas
  sal                 = para.sortie;
  ns                  =(sal.Ux + sal.Uy + sal.Uz)*sal.Ut;
  nss                 = sal.sxx + sal.syy + sal.szz + sal.sxy + sal.sxz + sal.syz;
  DWN.Udiff           = zeros(ns ,nbeq2,para.rec.m(1).nr);
  DWN.Sdiff           = zeros(nss,nbeq2,para.rec.m(1).nr);
  salur               = logical(ones(nrr,1)*sal.Ut);
  salsr               = logical(ones(nrr,1)*(nss>0));
  
  %concatenacion de los puntos virtuales y reales y concatenacion de las salidas
  xr                  = [xrv,xrr];
  yr                  = [yrv,yrr];
  zr                  = [zrv,zrr];
  zr0                 = [zrfv;zrfr];
  izr0                = [izrfv;izr0+length(zrfv)];
  salu                = [salu;salur];
  sals                = [sals;salsr];
  
  nzr0                = length(zr0);
  zricrall            = zeros(1,nzr0);
  icrall              = zeros(1,nzr0);
  for ir=1:nzr0
    icrv     = 1; % Estrato del receptor
    zricrv   = zr0(ir); % profundidad relativa a la interface de la capa del receptor
    while zricrv>para.reg(1).sub(icrv).h && icrv<para.nsubmed
      zricrv   = zricrv-para.reg(1).sub(icrv).h;
      icrv     = icrv+1;
    end
    zricrall(ir) = zricrv;
    icrall(ir)   = icrv;
  end
  
  %reduccion du numero de profundidades de los receptores
  DWN.xr              = xr;
  DWN.yr              = yr;
  DWN.zr              = zr;
  DWN.zr0             = zr0;
  DWN.izr0            = izr0;
  DWN.zricrall     	= zricrall;
  DWN.icrall        	= icrall;
  DWN.salu            = salu;
  DWN.sals            = sals;
  DWN.nxrv            = nxrv;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculo de los terminos independientes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% termino fuente conocido %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% esfuerzos nulos en la superficie libre
% continuidad de los esfuerzos y desplazamientos normales
% en los contornos fuera de la superficie libre

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vectores de condiciones  a las fronteras %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[B,DWN]	= vector_fuente3D(para,coord,DWN);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% matriz de coeficientes %
%%%%%%%%%%%%%%%%%%%%%%%%%%


% Elementos de la matriz correspondientes a regiones homogéneas:
A     = zeros(nbeq2,nbeq2);
% figure(192);
for i=1:coord.nbpt
  if i==4
    disp (i)
  end
  [A(i,:),...       %continuidad tx
   A(i+nbeq,:),...  %continuidad ty
   A(i+2*nbeq,:)]...%continuidad tz
   = Aij_Tn_3D(i,coord,para);
  if sum(coord.Xim(i,:)~=0)==2
    [A(ju(i)+nbpt,:),...       %continuidad ux
     A(ju(i)+nbpt+nbeq,:),...  %continuidad uy
     A(ju(i)+nbpt+2*nbeq,:)]...%continuidad uz
     = Aij_Gn_3D(i,coord,para);
  end
%   figure(192);spy(A);
end
%   % si Xi pertenece a un solo medio, entonces solo se considera la
%   % ecuacion de esfuerzos en superficie libre, si no, se considera
%   % la continuidad de los esfuerzos normales y desplazamientos.
% Adirect = A;
% A     = zeros(nbeq2,nbeq2);
%   % La matriz se llena para cada receptor (renglón i) a la vez.
% for i=1:coord.nbpt
%   A = Aij_Tn_3D_recip(A,nbeq,i,coord,para);
% %   figure(172);spy(A)
%   
% %   if sum(coord.Xim(i,:)~=0)==2
% %     [A(ju(i)+nbpt,:),...       %continuidad ux
% %      A(ju(i)+nbpt+nbeq,:),...  %continuidad uy
% %      A(ju(i)+nbpt+2*nbeq,:)]...%continuidad uz
% %      = Aij_Gn_3D(i,coord,para);
% %    figure(172);spy(A)
% %   end
% end
%   figure(172);spy(A)
% Adif = A-Adirect;
% a = 1:12; b = 11; [a.' A(a,b) Adirect(a,b) Adif(a,b)]
% a = 11:22; b=4; [a ; A(b,a) ; Adirect(b,a) ; Adif(b,a)]

% A     = zeros(nbeq2,nbeq2);
% AijxT = zeros(nbeq, nbeq2);
% AijyT = zeros(nbeq, nbeq2);
% AijzT = zeros(nbeq, nbeq2);
% AijxG = zeros(nbeq, nbeq2);
% AijyG = zeros(nbeq, nbeq2);
% AijzG = zeros(nbeq, nbeq2);
% 
% for i=1:coord.nbpt
%   [AijxT(i,:),AijyT(i,:),AijzT(i,:)]= Aij_Tn_3D(i,coord,para);
%   if sum(coord.Xim(i,:)~=0)==2
%     [AijxG(i,:),AijyG(i,:),AijzG(i,:)]= Aij_Gn_3D(i,coord,para);
%   end
% end
% 
% for i=1:coord.nbpt
%   %si Xi pertenece a un solo medio, entonces solo se considera la
%   %ecuacion de esfuerzos en superficie libre
%   %sino se considera la continuidad de los esfuerzos normales y de los
%   %desplazamientos
%   %     [Aijx,Aijy,Aijz]=Aij_Tn_3D(i,coord,para);
%   %continuidad tx
%   A(i,:)        = AijxT(i,:);% Aijx
%   %continuidad ty
%   A(i+nbeq,:)   = AijyT(i,:); %Aijy
%   %continuidad tz
%   A(i+2*nbeq,:) = AijzT(i,:); %Aijz
%   if sum(coord.Xim(i,:)~=0)==2
%     %         [Aijx,Aijy,Aijz]        = Aij_Gn_3D(i,coord,para);
%     %continuidad ux
%     A(ju(i)+nbpt,:)       = AijxG(i,:);% Aijx;
%     %continuidad uy
%     A(ju(i)+nbpt+nbeq,:)	= AijyG(i,:);% Aijy;
%     %continuidad uz
%     A(ju(i)+nbpt+2*nbeq,:)= AijzG(i,:);% Aijz;
%   end
% end
% clear AijxT AijyT AijzT AijxG AijyG AijzG

if para.geo(1)==3
  %en el caso del IBEM-DWN, es mas ventajoso, rellenar la matrice
  %considerando los phi, es a decir, columna por columna
  %tambien por cada punto de colocacion perteneciendo al multi-estratos
  %hay que considerar una continuidad de traciones y desplazamientos
  
  % cada fuente esta considerada por separada, la integracion se hace
  % bien, se tendria mejorar por fuentes de misma orientation, mismo z,
  % y misma longitud
  
  %hay que recordar que el medio 1 siempre es el medio de mas bajo
  %indice asi que el signo de la matriz ((m==im0) - (m~=im0))*t22(jj)
  %es obvio
  for j1=1:length(zrfv)
    j=indir(izrfv==j1); % indice de los puntos de colocación a la j1 profundidad
    nj=length(j);
    %phi = [1,0,0]; [0,1,0]; [0,0,1]
    rec.xr          = DWN.xr;
    rec.yr          = DWN.yr;
    rec.zr          = DWN.zr;
    rec.zr0         = DWN.zr0;
    rec.izr0        = DWN.izr0;
    rec.zricrall    = DWN.zricrall;
    rec.icrall      = DWN.icrall;
    nrec            = length(DWN.xr);
    
    ind_superpos    = false(1,nrec);
    ind_superpos(izrfv==j1)= true;
    signo1          = (coord.Xim(indir(ind_superpos),1)==1) - (coord.Xim(indir(ind_superpos),1)==2);
    
    %     str = toc; disp(['Elapsed time is ', num2str(str), ' seconds. : Pre_for  calcul_']) ;tic
    %     disp(' ');disp(['j1=',num2str(j1),'/',num2str(length(zrfv)),' |j|=',num2str(nj)])
    
    % repartir j entre varios procesos
    %     disp('j original:');disp(j)
    
    l = length(j); % cantidad de fuentes a esta profundidad
    
    % Variables locales (completas):
    U_f1      = zeros(3,  nrec,l);
    U_f2      = zeros(3,  nrec,l);
    U_f3      = zeros(3,  nrec,l);
    S_f1      = zeros(3,3,nrec,l);
    S_f2      = zeros(3,3,nrec,l);
    S_f3      = zeros(3,3,nrec,l);
    clear coordf
%     myCluster = parcluster('local');
    if l==1 % 
%     disp('single j');disp(j)
      coordf.xs       = coord.x(j);
      coordf.ys       = coord.y(j);
      coordf.zs       = coord.z(j);
      coordf.vnx      = coord.vnx(j);
      coordf.vny      = coord.vny(j);
      coordf.vnz      = coord.vnz(j);
      coordf.dr       = coord.drxz(j);
      if para.dim ~= 4
        coordf.x0       = coord.x0(j);
        coordf.th       = coord.th(j);
        coordf.drxz     = coord.drxz(j);
        coordf.nbptc    = coord.nbptc(j);
      else % para.dim == 4
        coordf.CubPtsex = coord.CubPtsex(:,:,j);
        coordf.CubPts   = coord.CubPts  (:,:,j);
      end
      [U_f1,S_f1,U_f2,S_f2,U_f3,S_f3] = calcul_US_DWN_3D_polar_Ncapas_HS2(...
                                         para,rec,salu,sals,...
                                         coordf,...
                                         ones(nrec,3),...
                                         DWN);
      
    else % el cálculo para varias fuentes se distribuye
      jDistD = distributed(j);
      %     disp(['numworks=' num2str(myCluster.NumWorkers)])
%       disp('j distribuida:');disp(jDistD)
      spmd %entonces coordf queda repartido y U_fic se apilaría por la dimension 3
        if labindex <= l
          jDist = getLocalPart(jDistD); %disp(['j at worker ' num2str(labindex) ':']); disp(jDist)
          coordf.xs       = coord.x(jDist);
          coordf.ys       = coord.y(jDist);
          coordf.zs       = coord.z(jDist);
          coordf.vnx      = coord.vnx(jDist);
          coordf.vny      = coord.vny(jDist);
          coordf.vnz      = coord.vnz(jDist);
          coordf.dr       = coord.drxz(jDist);
          if para.dim ~= 4
            coordf.x0       = coord.x0(jDist);
            coordf.th       = coord.th(jDist);
            coordf.drxz     = coord.drxz(jDist);
            coordf.nbptc    = coord.nbptc(jDist);
          else % para.dim == 4
            coordf.CubPtsex = coord.CubPtsex(:,:,jDist);
            coordf.CubPts   = coord.CubPts  (:,:,jDist);
          end
          %       disp('coordf:');disp(coordf);disp(['nsub:',num2str(para.nsubmed)]);disp('zs');disp(coordf.zs);disp('MAT');disp(para.reg(1).sub)
          [U_f1c,S_f1c,U_f2c,S_f2c,U_f3c,S_f3c] = calcul_US_DWN_3D_polar_Ncapas_HS2(para,rec,salu,sals,coordf,ones(nrec,3),DWN);
        end
      end %spmd
      % gather all composites (se transmite desde los procesadores)
      l = size(U_f1c,2); %cantidad de piezas
      l = min(l,length(j)); %con datos
      di3I = 1;
      for ic = 1:l
        di3F = di3I-1 + size(U_f1c{ic},3);
        U_f1(:,:  ,di3I:di3F) = U_f1c{ic};
        U_f2(:,:  ,di3I:di3F) = U_f2c{ic};
        U_f3(:,:  ,di3I:di3F) = U_f3c{ic};
        S_f1(:,:,:,di3I:di3F) = S_f1c{ic};
        S_f2(:,:,:,di3I:di3F) = S_f2c{ic};
        S_f3(:,:,:,di3I:di3F) = S_f3c{ic};
        di3I = di3F+1;
      end
    end
    
    U_f1 = permute(U_f1,[1 3 2]);
    S_f1 = permute(S_f1,[1 2 4 3]);
    U_f2 = permute(U_f2,[1 3 2]);
    S_f2 = permute(S_f2,[1 2 4 3]);
    U_f3 = permute(U_f3,[1 3 2]);
    S_f3 = permute(S_f3,[1 2 4 3]);
    %     str = toc; disp(['Elapsed time is ', num2str(str), ' seconds. : Calculate calcul_']) ;tic
    DWN.Udiff(:,coord.phi(j,1),:)         = U_f1(:,:,nxrv+1:nxrv+nrr);%(3  ,nbeq,para.rec.m(1).nr)
    if length(sals)>nxrv && sals(nxrv+1)
      S0  = squeeze(reshape(S_f1(:,:,:,nxrv+1:nrec),1,9,nj,nrec-nxrv));
      S0  = S0(DWN.inds,:,:);
      DWN.Sdiff(:,coord.phi(j,1),:)= S0;%(nss,nbeq,para.rec.m(1).nr)
    end
    U_f1    = squeeze(U_f1(:,:,1:nxrv));
    S_f1    = squeeze(S_f1(:,:,:,1:nxrv));
    
    DWN.Udiff(:,coord.phi(j,1)+nbeq,:)    = U_f2(:,:,nxrv+1:nxrv+nrr);
    if length(sals)>nxrv && sals(nxrv+1)
      S0	= squeeze(reshape(S_f2(:,:,:,nxrv+1:nrec),1,9,nj,nrec-nxrv));
      S0 	= S0(DWN.inds,:,:);
      DWN.Sdiff(:,coord.phi(j,1)+nbeq,:)= S0;%(nss,nbeq,para.rec.m(1).nr)
    end
    U_f2    = squeeze(U_f2(:,:,1:nxrv));
    S_f2	= squeeze(S_f2(:,:,:,1:nxrv));
    
    DWN.Udiff(:,coord.phi(j,1)+2*nbeq,:)  = U_f3(:,:,nxrv+1:nxrv+nrr);
    if length(sals)>nxrv && sals(nxrv+1)
      S0	= squeeze(reshape(S_f3(:,:,:,nxrv+1:nrec),1,9,nj,nrec-nxrv));
      S0 	= S0(DWN.inds,:,:);
      DWN.Sdiff(:,coord.phi(j,1)+2*nbeq,:)= S0;%(nss,nbeq,para.rec.m(1).nr)
    end
    U_f3    = squeeze(U_f3(:,:,1:nxrv));
    S_f3    = squeeze(S_f3(:,:,:,1:nxrv));
    
    vn1=repmat(vn(1,:).',1,nj);
    vn2=repmat(vn(2,:).',1,nj);
    vn3=repmat(vn(3,:).',1,nj);
    dAj=repmat(coord.dA(j),nxrv,1);
    
    if length(j)==1
      U_f1 = permute(U_f1,[1 3 2]);
      S_f1 = permute(S_f1,[1 2 4 3]);
      U_f2 = permute(U_f2,[1 3 2]);
      S_f2 = permute(S_f2,[1 2 4 3]);
      U_f3 = permute(U_f3,[1 3 2]);
      S_f3 = permute(S_f3,[1 2 4 3]);
      vn1=vn1.';
      vn2=vn2.';
      vn3=vn3.';
      dAj=dAj.';
    end
    
    diagindsup      =zeros(nj,2);
    ind0            =1:nxrv;
    diagindsup(:,1) =ind0(ind_superpos(ind0));
    diagindsup(:,2) =1:nj;
    diagindsup=sub2ind([nxrv,nj], diagindsup(:,1), diagindsup(:,2));
    
    %phix
    tn1         = squeeze(S_f1(1,1,:,:)).'.*vn1+squeeze(S_f1(2,1,:,:)).'.*vn2+squeeze(S_f1(3,1,:,:)).'.*vn3;
    tn2         = squeeze(S_f1(1,2,:,:)).'.*vn1+squeeze(S_f1(2,2,:,:)).'.*vn2+squeeze(S_f1(3,2,:,:)).'.*vn3;
    tn3         = squeeze(S_f1(1,3,:,:)).'.*vn1+squeeze(S_f1(2,3,:,:)).'.*vn2+squeeze(S_f1(3,3,:,:)).'.*vn3;
    %     tn1         = squeeze(S_f1(1,1,:,:)).*vn1+squeeze(S_f1(2,1,:,:)).*vn2+squeeze(S_f1(3,1,:,:)).*vn3;
    %     tn2         = squeeze(S_f1(1,2,:,:)).*vn1+squeeze(S_f1(2,2,:,:)).*vn2+squeeze(S_f1(3,2,:,:)).*vn3;
    %     tn3         = squeeze(S_f1(1,3,:,:)).*vn1+squeeze(S_f1(2,3,:,:)).*vn2+squeeze(S_f1(3,3,:,:)).*vn3;
    tn1(diagindsup)=signo1*0.5./coord.dA(j).';
    tn2(diagindsup)=0;
    tn3(diagindsup)=0;
    %continuidad tx
    A(indir                 ,coord.phi(j,1))= tn1.*dAj;
    %continuidad ty
    A(indir+nbeq            ,coord.phi(j,1))= tn2.*dAj;
    %continuidad ty
    A(indir+nbeq*2          ,coord.phi(j,1))= tn3.*dAj;
    if min(ju(indir))>0
      %continuidad ux
      A(ju(indir)+nbpt        ,coord.phi(j,1))= squeeze(U_f1(1,:,:)).'.*dAj;
      %continuidad uy
      A(ju(indir)+nbpt+nbeq   ,coord.phi(j,1))= squeeze(U_f1(2,:,:)).'.*dAj;
      %continuidad uy
      A(ju(indir)+nbpt+nbeq*2 ,coord.phi(j,1))= squeeze(U_f1(3,:,:)).'.*dAj;
    end
    
    %phiy
    tn1         = squeeze(S_f2(1,1,:,:)).'.*vn1+squeeze(S_f2(2,1,:,:)).'.*vn2+squeeze(S_f2(3,1,:,:)).'.*vn3;
    tn2         = squeeze(S_f2(1,2,:,:)).'.*vn1+squeeze(S_f2(2,2,:,:)).'.*vn2+squeeze(S_f2(3,2,:,:)).'.*vn3;
    tn3         = squeeze(S_f2(1,3,:,:)).'.*vn1+squeeze(S_f2(2,3,:,:)).'.*vn2+squeeze(S_f2(3,3,:,:)).'.*vn3;
    tn1(diagindsup)=0;
    tn2(diagindsup)=signo1*0.5./coord.dA(j).';
    tn3(diagindsup)=0;
    %continuidad tx
    A(indir                 ,coord.phi(j,1)+nbeq)= tn1.*dAj;
    %continuidad ty
    A(indir+nbeq            ,coord.phi(j,1)+nbeq)= tn2.*dAj;
    %continuidad ty
    A(indir+nbeq*2          ,coord.phi(j,1)+nbeq)= tn3.*dAj;
    if min(ju(indir))>0
      %continuidad ux
      A(ju(indir)+nbpt        ,coord.phi(j,1)+nbeq)= squeeze(U_f2(1,:,:)).'.*dAj;
      %continuidad uy
      A(ju(indir)+nbpt+nbeq   ,coord.phi(j,1)+nbeq)= squeeze(U_f2(2,:,:)).'.*dAj;
      %continuidad uy
      A(ju(indir)+nbpt+nbeq*2 ,coord.phi(j,1)+nbeq)= squeeze(U_f2(3,:,:)).'.*dAj;
    end
    
    %phiz
    tn1         = squeeze(S_f3(1,1,:,:)).'.*vn1+squeeze(S_f3(2,1,:,:)).'.*vn2+squeeze(S_f3(3,1,:,:)).'.*vn3;
    tn2         = squeeze(S_f3(1,2,:,:)).'.*vn1+squeeze(S_f3(2,2,:,:)).'.*vn2+squeeze(S_f3(3,2,:,:)).'.*vn3;
    tn3         = squeeze(S_f3(1,3,:,:)).'.*vn1+squeeze(S_f3(2,3,:,:)).'.*vn2+squeeze(S_f3(3,3,:,:)).'.*vn3;
    tn1(diagindsup)=0;
    tn2(diagindsup)=0;
    tn3(diagindsup)=signo1*0.5./coord.dA(j).';
    %continuidad tx
    A(indir                 ,coord.phi(j,1)+2*nbeq)= tn1.*dAj;
    %continuidad ty
    A(indir+nbeq            ,coord.phi(j,1)+2*nbeq)= tn2.*dAj;
    %continuidad ty
    A(indir+nbeq*2          ,coord.phi(j,1)+2*nbeq)= tn3.*dAj;
    if min(ju(indir))>0
      %continuidad ux
      A(ju(indir)+nbpt        ,coord.phi(j,1)+2*nbeq)= squeeze(U_f3(1,:,:)).'.*dAj;
      %continuidad uy
      A(ju(indir)+nbpt+nbeq   ,coord.phi(j,1)+2*nbeq)= squeeze(U_f3(2,:,:)).'.*dAj;
      %continuidad uy
      A(ju(indir)+nbpt+nbeq*2 ,coord.phi(j,1)+2*nbeq)= squeeze(U_f3(3,:,:)).'.*dAj;
    end
    %     str = toc; disp(['Elapsed time is ', num2str(str), ' seconds. : Distribute_calc'])
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% resolucion del sistema %
%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_fv =zeros(nbeq2,para.ninc);

if para.ninc>1
  A1=inv(A);
  for iinc=1:para.ninc
    phi_fv(:,iinc)=A1*B(:,iinc); %#ok<*MINV>
    %regularizacion de phi_fv
  end
else
  
  phi_fv=A\B;
  %         regularizacion_phi;
  
  %     if issparse(A), R = qr(A);
  %     else R = triu(qr(A)); end
  %     for iinc=1:para.ninc
  %         phi_fv(:,iinc) = R\(R'\(A'*B(:,iinc)));
  %         r = B(:,iinc) - A*phi_fv(:,iinc);
  %         e = R\(R'\(A'*r));
  %         phi_fv(:,iinc) = phi_fv(:,iinc) + e;
  %         regularizacion_phi;
  %     end
  %   error('DebugBrake:  Fin de sol3D_dfv')
  % [char(8)*ones(1,lStr+lPrompt),str]  % borrar el texto sobre el tamaño del sistema
end
end