function [Udiff,Sdiff]=int_fv_Gij_3D(phi_fv,coordf,para,m,xr,yr,zr,gaussian,kpi,ksi,C)

%%%%%%%%%%%%%
% unpacking %
%%%%%%%%%%%%%
ninc    = para.ninc;
sal     = para.sortie;
ns      = (sal.Ux + sal.Uy + sal.Uz)*sal.Ut;
nss   	= sal.sxx + sal.syy + sal.szz+sal.sxy+sal.sxz+sal.syz;

%posicion de las fuentes virtuales
xf      = coordf.x;
yf      = coordf.y;
zf      = coordf.z;

dr      = coordf.drxz;
dA      = coordf.dA;
phi     = coordf.phi;
nbeq    = coordf.nbeq;

%init el campo difractado
Udiff   = zeros(ns ,ninc);
Sdiff   = zeros(nss,ninc);


%posicion de las fv que hay que tomar en cuenta
j       = 1:coordf.nbpt;
jphi    = 1:coordf.nbeq;

%indice (logic y natural) de los puntos de colocacion perteneciendo a m
jjx     = coordf.indm(m).ind;
ii      = j(jjx);

%indice (logic y natural) de los phi que hay que contemplar (columnas
%de la matriz)
jjphi   = false(coordf.nbeq,1);
jjphi(phi(j(jjx),m)) = true(1);
jj      = jphi(jjphi);
njj     = length(jj);

%calculo de la suma de cada una de las contribuciones
xij     = xr-xf(ii);
yij     = yr-yf(ii);
zij     = zr-zf(ii);
drj     = dr(ii);
dAj     = dA(ii);

rij     = sqrt(xij.^2+yij.^2+zij.^2);
g       = zeros(3,njj);
g(1,:)  = xij./rij;
g(2,:)  = yij./rij;
g(3,:)  = zij./rij;

%integracion gaussiana donde se requiere
if para.dim == 4 % 3Dgeneral
jjtmp2  = find((rij<=2.0*drj).*(rij~=0)); % a menos de 1.5 radios  
else %3D axisimétrico
jjtmp2  = find(rij<para.npplo/2*drj);
end
if sal.Ut==1
    Gij             = Gij_3D(ksi,kpi,rij,g,C,njj);
    if para.dim == 4
    ex = false;
    Gij(:,:,jjtmp2) = Gij_3D_r_smallGEN(coordf,[xr,yr,zr],ii(jjtmp2),ksi,kpi,para,C,ex);
    jjtmp2              = find(rij==0);
    ex = true;
    Gij(:,:,jjtmp2) = Gij_3D_r_smallGEN(coordf,[xr,yr,zr],ii(jjtmp2),ksi,kpi,para,C,ex);
    else
    Gij(:,:,jjtmp2) = Gij_3D_r_small(coordf,xr,yr,zr,ii(jjtmp2),ksi,kpi,gaussian,C);
    end
end
if sal.sxx==1 || sal.syy==1 || sal.szz==1 || sal.sxy==1 || sal.sxz==1 || sal.syz==1 
    [S_fx,S_fy,S_fz]    = S_3D(rij,g,ksi,kpi,[1 1 1]);
    [S_fx(:,:,jjtmp2),S_fy(:,:,jjtmp2),S_fz(:,:,jjtmp2)] = ...
        S_3D_r_small(coordf,jjtmp2,xr,yr,zr,ksi,kpi,[1 1 1],gaussian);
end
    

for iinc=1:ninc
    k = 1;
    if sal.Ut==1
        if sal.Ux==1
            Udiff(k,iinc)=sum(dAj.'.*(...
                squeeze(Gij(1,1,:)).*phi_fv(jj       ,iinc)+...
                squeeze(Gij(1,2,:)).*phi_fv(jj+  nbeq,iinc)+...
                squeeze(Gij(1,3,:)).*phi_fv(jj+2*nbeq,iinc)));
            
            k=k+1;
        end
         if sal.Uy==1
            Udiff(k,iinc)=sum(dAj.'.*(...
                squeeze(Gij(2,1,:)).*phi_fv(jj       ,iinc)+...
                squeeze(Gij(2,2,:)).*phi_fv(jj+  nbeq,iinc)+...
                squeeze(Gij(2,3,:)).*phi_fv(jj+2*nbeq,iinc)));
            k=k+1;
        end
        if sal.Uz==1
            Udiff(k,iinc)=sum(dAj.'.*(...
                squeeze(Gij(3,1,:)).*phi_fv(jj       ,iinc)+...
                squeeze(Gij(3,2,:)).*phi_fv(jj+  nbeq,iinc)+...
                squeeze(Gij(3,3,:)).*phi_fv(jj+2*nbeq,iinc)));
        end
    end
    k = 1;
    if sal.sxx==1
        Sdiff(k,iinc)=sum(dAj.'.*(...
            squeeze(S_fx(1,1,:)).*phi_fv(jj       ,iinc)+...
            squeeze(S_fy(1,1,:)).*phi_fv(jj+  nbeq,iinc)+...
            squeeze(S_fz(1,1,:)).*phi_fv(jj+2*nbeq,iinc)));
        k=k+1;
    end
    if sal.syy==1
        Sdiff(k,iinc)=sum(dAj.'.*(...
            squeeze(S_fx(2,2,:)).*phi_fv(jj       ,iinc)+...
            squeeze(S_fy(2,2,:)).*phi_fv(jj+  nbeq,iinc)+...
            squeeze(S_fz(2,2,:)).*phi_fv(jj+2*nbeq,iinc)));
        k=k+1;
    end
    if sal.szz==1
        Sdiff(k,iinc)=sum(dAj.'.*(...
            squeeze(S_fx(3,3,:)).*phi_fv(jj       ,iinc)+...
            squeeze(S_fy(3,3,:)).*phi_fv(jj+  nbeq,iinc)+...
            squeeze(S_fz(3,3,:)).*phi_fv(jj+2*nbeq,iinc)));
        k=k+1;
    end
    if sal.sxy==1
        Sdiff(k,iinc)=sum(dAj.'.*(...
            squeeze(S_fx(1,2,:)).*phi_fv(jj       ,iinc)+...
            squeeze(S_fy(1,2,:)).*phi_fv(jj+  nbeq,iinc)+...
            squeeze(S_fz(1,2,:)).*phi_fv(jj+2*nbeq,iinc)));
        k=k+1;
    end
    if sal.sxz==1
        Sdiff(k,iinc)=sum(dAj.'.*(...
            squeeze(S_fx(1,3,:)).*phi_fv(jj       ,iinc)+...
            squeeze(S_fy(1,3,:)).*phi_fv(jj+  nbeq,iinc)+...
            squeeze(S_fz(1,3,:)).*phi_fv(jj+2*nbeq,iinc)));
        k=k+1;
    end
    if sal.syz==1
        Sdiff(k,iinc)=sum(dAj.'.*(...
            squeeze(S_fx(2,3,:)).*phi_fv(jj       ,iinc)+...
            squeeze(S_fy(2,3,:)).*phi_fv(jj+  nbeq,iinc)+...
            squeeze(S_fz(2,3,:)).*phi_fv(jj+2*nbeq,iinc)));
    end
end