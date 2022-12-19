function [A]=Aij_Tn_3D_recip(A,nbeq,i,coord,para)

% coeficientes de la matriz A correspondiente a la continuidad de
% los esfuerzos normales
% i : indice de l element (de surface) sur lequel se calcul la contribution de toutes
% les autres sources virtuelles j (ponctuelle)


%%%%%%%%%%%%%
% unpacking %
%%%%%%%%%%%%%

x 	=coord.x;
y 	=coord.y;
z 	=coord.z;
vnx	=coord.vnx;
vny	=coord.vny;
vnz	=coord.vnz;
dr  =coord.drxz;
dA	=coord.dA;
phi =coord.phi;


j       = 1:coord.nbpt;
jphi    = 1:coord.nbeq;
tij     = zeros(3,3,coord.nbeq);
tijR    = zeros(3,3,coord.nbeq);

gaussian   = para.gaussian;

%indice de los medios en contactos
im  = find(coord.Xim(i,:)~=0);

% longitud de onda minima


%min de los indice de los medios en contacto que se lleva el signo +
im0 = min(im);

for m=im;
  if m==1 && para.geo(1)==3 % si es un fondo estratificado dwn
    continue
  end
  if isfield(para,'subm1')
    mmat= para.subm1(m);
  else
    mmat = m;
  end
  Ci  = para.reg(mmat).Ci;
  ksi = para.reg(mmat).ksi;
  kpi = para.reg(mmat).kpi;
  
  %indice (logic y natural) de los puntos de colocacion perteneciendo a m
  jjx     = coord.indm(m).ind;
  jr      = j(jjx);
  
  %indice (logic y natural) de los phi que hay que contemplar (columnas
  %de la matriz)
  jjphi   = false(coord.nbeq,1);
  jjphi(phi(j(jjx),m)) = true(1);
  jj      = jphi(jjphi);
  njj     = length(jj);
  
  drj     = dr(jjx);
  dAj     = dA(jjx);
  xij     = x(i)-x(jjx);
  yij     = y(i)-y(jjx);
  zij     = z(i)-z(jjx);
  rij     = sqrt(xij.^2+yij.^2+zij.^2);
  g       = zeros(3,njj);
  g(1,:)  = xij./rij;
  g(2,:)  = yij./rij;
  g(3,:)  = zij./rij;
  %                                                             índice
  % vector normal del punto de colocación (x)                  ::  i
  % vectores normales de los puntos receptores (x') recíprocos :: jjxR
  vn(1,:) = vnx;
  vn(2,:) = vny;
  vn(3,:) = vnz;
  
  % Los que requerirán integración gaussiana:
  if para.dim == 4 % 3Dgeneral
    jjtmp2  = find((rij<=1.5*drj(i)).*(rij~=0)); % a menos de 1.5 radios
  else %3D axisimétrico
    jjtmp2  = find((rij<para.npplo/2*drj).*(rij~=0));
  end
  if i>1
    jjtmp2(jjtmp2<=(i-1)) = []; % tampoco incluir las fuentes que ya se calcularon
  end
  %  ver:
  %  a =i; hold on; h_dot = plot3(x(a),y(a),z(a),'r*');
  %  a =jjtmp2; hold on; h_dot = plot3(x(a),y(a),z(a),'w*');
  %  delete(h_dot)
  
  jjxFar = jjx;
  jjxFar(rij==0) = false; % no incluir el que sabemos que es +-1/2
  if i>1
    jjxFar(1:i-1) = false; % tampoco incluir las fuentes que ya se calcularon
  end
  % no hacer los que requerirán integración gaussiana:
  if ~isempty(jjtmp2); jjxFar(jjtmp2) = false; end
  jjphi(phi(j(~jjxFar),m)) = false(1);
  
  % Uso el teorema de reciprocidad de las funciones de Green
  %    Aprovechamos los resultados para el renglón i, columnas jjphi
  %    para los elementos de los renglones jjphiR, columna i:
  
  jjxR = jjxFar; jjxR(1:i) = false; % no calcular recíproco de lo que ya está calculado
  %     jrR      = j(jjxR);
  jjphiR   = false(coord.nbeq,1);
  jjphiR(phi(j(jjxR),m)) = true(1); % (renglones de la matriz)
  jjphiRA = jjphiR;
  jjphiRA(jj(jjtmp2)) = true;
  % [(1:22).' jjphi jjphiR jjphiRA]
  [tij(:,:,jjphi),tijR(:,:,jjphiR)]=Tij_3D(ksi,kpi,rij(jjxFar),g(:,jjxFar),Ci,vn,i,jjxR);
  for i0=1:3 % multiplicar por el área de las fuentes
      for j0=1:3
        tij (i0,j0,jjphi)  = squeeze(tij (i0,j0,jjphi)).*dAj(jjxFar).';
        tijR(i0,j0,jjphiR) = squeeze(tijR(i0,j0,jjphi))*dAj(i);
      end
  end
    
  %integracion gaussiana donde se requiere        ,--- receptor
  if para.dim == 4 % 3Dgeneral                    |  ,--- fuentes
    tij(:,:,jj(jjtmp2)) = Tij_3D_r_smallGEN(coord,i,jr(jjtmp2),ksi,kpi,para,Ci);
    for i0=1:3
      for j0=1:3
        tij(i0,j0,jj(jjtmp2)) = squeeze(tij(i0,j0,jj(jjtmp2))).*dAj(jjtmp2).';
      end
    end%                                         ,--- receptor
    for h = jjtmp2 % los recíprocos              |    ,--- fuente
      if i == 2
        disp('here')
      end
      tijR(:,:,jj(h)) = Tij_3D_r_smallGEN(coord,jr(h),i,ksi,kpi,para,Ci);
      for i0=1:3
        for j0=1:3
          tijR(i0,j0,jj(h)) = squeeze(tijR(i0,j0,jj(h)))*dAj(i);
        end
      end
    end
  else %3D axisimétrico
    tij(:,:,jj(jjtmp2)) = Tij_3D_r_small(coord,i,jr(jjtmp2),ksi,kpi,gaussian,Ci);
    for i0=1:3
      for j0=1:3
        tij(i0,j0,jj) = squeeze(tij(i0,j0,jj)).*dAj.';
      end
    end
  end
  
  
  % cuando rij=0
  % en este punto, hay discontinuidad de las tractiones ti-ti'=phi
  % el signo1 es igual a +1 (resp. -1) si la normal apunta hacia el exterior
  % (resp. a interior) del medio a que pertenece el phi del punto de
  % colocacion en la posicion singular rij=0
  % como la convencion es que las normales esten siempre dirigidas hacia
  % abajo, hay que conocer si el phi pertenece a un contorno arriba o
  % abajo. Esta informacion esta en Xim
  signo1= (coord.Xim(i,m)==1) - (coord.Xim(i,m)==2);
  tij(:,:,jj(rij==0))  = signo1*0.5*eye(3);
  
  % en un punto de colocacion, hay siempre 2 medios en contacto,
  % aunque a veces este medio es el vacio
  % las ecuaciones de continuidad de traccion o de desplazamiento
  % involucran siempre 2 medios  : sigma(m)=sigma(m1) o u(m)=u(m1)
  % y como cada campo corresponde a (por ejemplo) u=u_diff+u_inc
  % se reorganiza entonces como :
  % sigma_diff(m)-sigma_diff(m1)=u_inc(m1)-u_inc(m)
  % los signos se fijan deacuerdo con el indice del medio,
  % + si el indice es el mas pequeno de los medios en contacto
  % - en lo contrario
  tij(:,:,jj) = ((m==im0) - (m~=im0))*tij(:,:,jj);
  tijR(:,:,jj) = ((m==im0) - (m~=im0))*tijR(:,:,jj);
end
% Con las funciones de Green calculadas directamente,
% rellenar renglón de la matriz:
colx = phi(j(i),m):phi(j(i),m)+sum(jjxFar)+length(jjtmp2);
coly = phi(j(i),m)+nbeq:phi(j(i),m)+sum(jjxFar)+length(jjtmp2)+nbeq;
colz = phi(j(i),m)+2*nbeq:phi(j(i),m)+sum(jjxFar)+length(jjtmp2)+2*nbeq;

A(i,colx)=squeeze(tij(1,1,phi(j(i),m):end));%t_x=t_x1.Phi_j
A(i,coly)=squeeze(tij(1,2,phi(j(i),m):end));%t_x=t_x2.Phi_j
A(i,colz)=squeeze(tij(1,3,phi(j(i),m):end));%t_x=t_x3.Phi_j

A(i+nbeq,colx)=squeeze(tij(2,1,phi(j(i),m):end));%t_y=t_y1.Phi_j
A(i+nbeq,coly)=squeeze(tij(2,2,phi(j(i),m):end));%t_y=t_y2.Phi_j
A(i+nbeq,colz)=squeeze(tij(2,3,phi(j(i),m):end));%t_y=t_y3.Phi_j

A(i+2*nbeq,colx)=squeeze(tij(3,1,phi(j(i),m):end));%t_z=t_z1.Phi_j
A(i+2*nbeq,coly)=squeeze(tij(3,2,phi(j(i),m):end));%t_z=t_z2.Phi_j
A(i+2*nbeq,colz)=squeeze(tij(3,3,phi(j(i),m):end));%t_z=t_z3.Phi_j

% Con las funciones de Green recíprocas,
% rellenar columnas de la matriz:
rens = i+1:i+sum(jjxR)+length(jjtmp2);
col = phi(j(i),m);
A(rens,col)        = squeeze(tijR(1,1,jjphiRA));
A(rens,col+nbeq)   = squeeze(tijR(1,2,jjphiRA));
A(rens,col+2*nbeq) = squeeze(tijR(1,3,jjphiRA));

rens = i+nbeq+1:i+nbeq+sum(jjxR)+length(jjtmp2);
A(rens,col)        = squeeze(tijR(2,1,jjphiRA));
A(rens,col+nbeq)   = squeeze(tijR(2,2,jjphiRA));
A(rens,col+2*nbeq) = squeeze(tijR(2,3,jjphiRA));

rens = i+2*nbeq+1:i+2*nbeq+sum(jjxR)+length(jjtmp2);
A(rens,col)        = squeeze(tijR(3,1,jjphiRA));
A(rens,col+nbeq)   = squeeze(tijR(3,2,jjphiRA));
A(rens,col+2*nbeq) = squeeze(tijR(3,3,jjphiRA));
%    figure(172);spy(A)
end