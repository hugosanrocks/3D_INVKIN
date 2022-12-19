function [Aijx,Aijz]=Aij_Gn_F(i,coord,para)

% coeficientes de la matriz A correspondiente a la continuidad de los desplazamientos
% i : indice de l element (de surface) sur lequel se calcul la contribution de toutes
% les autres sources virtuelles j (ponctuelle)

%%%%%%%%%%%%%
% unpacking %
%%%%%%%%%%%%%

x 	= coord.x;
z 	= coord.z;
vnx	= coord.vnx;
vnz	= coord.vnz;
dr	= coord.dr;
phi = coord.phi;

j       = 1:coord.nbpt;
jphi    = 1:coord.nbeq;
Gij     = zeros(2,2,coord.nbeq);

gaussian   = para.gaussian;

%indice de los medios en contactos
im  = find(coord.Xim(i,:)~=0);

%min de los indice de los medios en contacto que se lleva el signo +
im0 = min(im);

for m=im
    if m==1 && para.geo(1)==3
        continue
    end
    mmat= para.subm1(m);
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
    xij     = x(i)-x(jjx);
    zij     = z(i)-z(jjx);
    rij     = sqrt(xij.^2+zij.^2);
    g       = zeros(2,njj);
    g(1,:)  = xij./rij;
    g(2,:)  = zij./rij;
    
    if Ci(6,6)==0
        Gij(1,1,jj)  = G22_SH(kpi,rij,1);
        Gij(2,2,jj)  = Gij(1,1,jj);
    else
        Gij(:,:,jj)  = Gij_PSV(ksi,kpi,rij,g,Ci,njj);
    end
    
    %integracion gaussiana donde se requiere
    jjtmp2=find((rij<=para.npplo*drj).*(rij>0));
    if Ci(6,6)==0
        Gij(1,1,jj(jjtmp2))  = G22_SH_r_small(coord,x(i),z(i),jr(jjtmp2),kpi,gaussian,1);
        Gij(2,2,jj(jjtmp2))  = Gij(1,1,jj);
    else
        Gij(:,:,jj(jjtmp2)) = Gij_PSV_r_small(coord,x(i),z(i),jr(jjtmp2),ksi,kpi,gaussian,Ci);
    end
    %cuando rij=0
    jjtmp2=find(rij==0);
    %     vn(1)=vnx(jr(jjtmp2));
    %     vn(2)=vnz(jr(jjtmp2));
    %     Gij(:,:,jj(jjtmp2)) = Gij_PSV_r_eq_0(ksi,kpi,0,drj(jjtmp2),Ci,vn);
    
    gn(1)=-vnz(jr(jjtmp2));
    gn(2)= vnx(jr(jjtmp2));
    if Ci(6,6)==0
        Gij(1,1,jj(jjtmp2)) = G22_SH_r_eq_0(kpi,drj(jjtmp2),1);
        Gij(2,2,jj)  = Gij(1,1,jj);
    else
        Gij(:,:,jj(jjtmp2)) = Greenex_PSV(ksi,kpi,gn,Ci,drj(jjtmp2));
    end
    
    %segun la normal
    for i0=1:2
        for j0=1:2
            Gij(i0,j0,jj) = ((m==im0) - (m~=im0))*squeeze(Gij(i0,j0,jj)).*drj.';
        end
    end
    Aijx=[squeeze(Gij(1,1,:));squeeze(Gij(1,2,:))];%u_i=G_ij.Phi_j
    Aijz=[squeeze(Gij(2,1,:));squeeze(Gij(2,2,:))];%u_i=G_ij.Phi_j
end