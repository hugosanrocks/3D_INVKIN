function [Aijx,Aijy,Aijz]=Aij_Gn_3D(i,coord,para)

% coeficientes de la matriz A correspondiente a la continuidad de los desplazamientos
% i : indice de l element (de surface) sur lequel se calcul la contribution de toutes
% les autres sources virtuelles j (ponctuelle)

%%%%%%%%%%%%%
% unpacking %
%%%%%%%%%%%%%

x 	= coord.x;
y 	= coord.y;
z 	= coord.z;

xr  = coord.x(i);
yr  = coord.y(i);
zr  = coord.z(i);

dr	= coord.drxz;
dA	= coord.dA;
phi = coord.phi;

j       = 1:coord.nbpt;
jphi    = 1:coord.nbeq;
Gij     = zeros(3,3,coord.nbeq);

gaussian   = para.gaussian;
if para.dim == 3
gaussex = para.gaussex;
end

%indice de los medios en contactos
im  = find(coord.Xim(i,:)~=0);

%min de los indice de los medios en contacto que se lleva el signo +
im0 = min(im);

for m=im
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
    xij     = xr-x(jjx);
    yij     = yr-y(jjx);
    zij     = zr-z(jjx);
    rij     = sqrt(xij.^2+yij.^2+zij.^2);
    g       = zeros(3,njj);
    g(1,:)  = xij./rij;
    g(2,:)  = yij./rij;
    g(3,:)  = zij./rij;
    
    Gij(:,:,jj)         = Gij_3D(ksi,kpi,rij,g,Ci,njj);
    
    %integracion gaussiana donde se requiere
    if para.dim == 4 % 3Dgeneral
%     jjtmp2              = find((rij<=para.npplo/2*drj).*(rij~=0));
    jjtmp2  = find((rij<=1.5*drj).*(rij~=0)); % a menos de 1.5 radios
    ex = false;
    Gij(:,:,jj(jjtmp2)) = Gij_3D_r_smallGEN(coord,[xr,yr,zr],jr(jjtmp2),ksi,kpi,para,Ci,ex);
    
    jjtmp2              = find(rij==0);
    ex = true;
    Gij(:,:,jj(jjtmp2)) = Gij_3D_r_smallGEN(coord,[xr,yr,zr],jr(jjtmp2),ksi,kpi,para,Ci,ex);
    else %3D axisimétrico
    jjtmp2              = find((rij<=para.npplo/2*drj).*(rij~=0));
    Gij(:,:,jj(jjtmp2)) = Gij_3D_r_small(coord,xr,yr,zr,jr(jjtmp2),ksi,kpi,gaussian,Ci);
    
    jjtmp2              = find(rij==0);
    Gij(:,:,jj(jjtmp2)) = Gij_3D_r_small(coord,xr,yr,zr,jr(jjtmp2),ksi,kpi,gaussex,Ci);
    end 
    %segun la normal
    for i0=1:3
        for j0=1:3
            Gij(i0,j0,jj) = ((m==im0) - (m~=im0))*squeeze(Gij(i0,j0,jj)).*dAj.';
        end
    end
    Aijx=[squeeze(Gij(1,1,:));squeeze(Gij(1,2,:));squeeze(Gij(1,3,:))];%u_i=G_ij.Phi_j
    Aijy=[squeeze(Gij(2,1,:));squeeze(Gij(2,2,:));squeeze(Gij(2,3,:))];%u_i=G_ij.Phi_j
    Aijz=[squeeze(Gij(3,1,:));squeeze(Gij(3,2,:));squeeze(Gij(3,3,:))];%u_i=G_ij.Phi_j
end