function [uxzdiff,szxdiff]=int_fv_Gij_PSV(phi_fv,coordf,para,m,xr,zr,gaussian,kpi,ksi,C)
%integracion de la contribucion de todas las fuentes virtuales

%%%%%%%%%%%%%
% unpacking %
%%%%%%%%%%%%%
ninc    = para.ninc;
sal     = para.sortie;
ns      = (sal.Ux + sal.Uz)*(sal.Ut+sal.UPh+sal.USh+sal.UIh+sal.UPt+sal.USt);
nss     = (sal.sxx + sal.szz + sal.sxz);
%init el campo difractado
uxzdiff  = zeros(ns,ninc);
szxdiff  = zeros(nss,ninc);

%posicion de las fuentes virtuales
xf      = coordf.x;
zf      = coordf.z;
dr      = coordf.dr;
phi     = coordf.phi;
nbeq    = coordf.nbeq;

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
zij     = zr-zf(ii);
drj     = dr(ii);
rij     = sqrt(xij.^2+zij.^2);
g       = zeros(2,njj);
g(1,:)  = xij./rij;
g(2,:)  = zij./rij;
if (sal.sxx || sal.sxz || sal.szz); anyStress = true;else anyStress = false;end
if C(6,6)~=0 % rho bet^2 = mu ~= 0 means it is a solid
    Gij     = Gij_PSV(ksi,kpi,rij,g,C,njj);
    if(anyStress) 
      Sij = Sij_PSV(rij,g,ksi,kpi,C);
    end
    %integracion gaussiana donde se requiere
    jjtmp2  = find(rij<3*para.npplo*drj);
    if ~isempty(jjtmp2)
      Gij(:,:,jjtmp2)=Gij_PSV_r_small_FP(coordf,xr,zr,jjtmp2,ksi,kpi,gaussian,C);
      if(anyStress)
        Sij(:,:,jjtmp2) = Sij_PSV_r_small_FP(coordf,xr,zr,jjtmp2,ksi,kpi,gaussian,C);
      end
    end
else % es un fluido
    Gij         = zeros(2,2,njj);
    Gij0        = 1i/4*kpi.*besselh(1,2,kpi.*rij);
    Gij(1,1,:)  = Gij0.*g(1,:);
    Gij(2,1,:)  = Gij0.*g(2,:);
    %TODO esfuerzos
    if(anyStress); error('falan los esfuerzos');end
    %integracion gaussiana donde se requiere
    jjtmp2              = find((rij<=para.npplo*drj));
    jr                  = j(jjx);
    
    coordf.x(coordf.nbpt+1)     = xr;
    coordf.z(coordf.nbpt+1)     = zr;
    coordf.vnx(coordf.nbpt+1)   = 1;
    coordf.vnz(coordf.nbpt+1)   = 0;
    Gij(1,1,jjtmp2) = T22_SH_r_small(coordf,coordf.nbpt+1,jr(jjtmp2),kpi,gaussian);
    coordf.vnx(coordf.nbpt+1)   = 0;
    coordf.vnz(coordf.nbpt+1)   = 1;
    Gij(2,1,jjtmp2) = T22_SH_r_small(coordf,coordf.nbpt+1,jr(jjtmp2),kpi,gaussian);
end

if sal.UPh==1 || sal.USh==1 || sal.UIh==1
    % Separacion del campo de desplazamiento total en 3 contribuciones: P
    % homogeneo, SV homogeneo, y ondas inhomogeneas (k>kp y >ks).
    % Se obtiene a traves de la descomposicion de las funciones de Hankel
    % en ondas planas homogeneas e inhomogeneas
    % Les fonctions de green sont exprimees ds le ref normal de l IBEM grace
    % aux g alors que la decomposition des fonction de Hankel est ds le
    % repere associe aux elements, ce qui permet prendre en compte la
    % longueur de la source facilement par convolution en kx
    
    % angle entre le nouveau repere lie a l element et le repere IBEM
    vnx   	= coordf.vnx(ii);
    vnz 	= coordf.vnz(ii);
    th1     = atan2(vnz,vnx)+pi/2;%vnz est opposee a z-pi/2-pi=pi/2(2.pi)
    th0     = atan2(zij,xij);
    xij1    = rij.*cos(th0-th1);
    zij1    = rij.*sin(th0-th1);
    %     g1(1,:)  = xij1./rij;
    %     g1(2,:)  = zij1./rij;
    [Gpr,Gsr]	= Gij_P_SV_R(ksi,kpi,xij1,zij1,rij,g,C,drj);
    Gi      = Gij-(Gpr+Gsr);    
end

if sal.UPt==1 || sal.USt==1
    % Separacion de los campos de polarizacion P y S 
    Gpt             = Gij_Pt_SVt(kpi,rij,g,C,njj);
    if ~isempty(jjtmp2)
        Gpt(:,:,jjtmp2) = Gij_Pt_SVt_r_small(coordf,xr,zr,ii(jjtmp2),kpi,gaussian,C);
    end
    Gst             = Gij-Gpt;
end

jjz=coordf.indphiz0(jj);
jjz(jjz==0)=[];
for iinc=1:ninc
    k       = 1;
    ks      = 1;
    if sal.Ut==1
        if sal.Ux==1
            uxzdiff(k,iinc)=sum(squeeze(Gij(1,1,:)).*drj.'.*phi_fv(jj,iinc));
%             uxzdiff(k,iinc)=sum(squeeze(Gij(1,1,:)).*drj.'.*phi_fv(jj,iinc)+squeeze(Gij(1,2,:)).*drj.'.*phi_fv(jj+nbeq,iinc));
            if ~isempty(jjz)
                uxzdiff(k,iinc)=uxzdiff(k,iinc)+sum(squeeze(Gij(1,2,:)).*drj.'.*phi_fv(jjz+nbeq,iinc));
            end
            k=k+1;
        end
        if sal.Uz==1
            uxzdiff(k,iinc)=sum(squeeze(Gij(2,1,:)).*drj.'.*phi_fv(jj,iinc));
            if ~isempty(jjz)
                uxzdiff(k,iinc)=uxzdiff(k,iinc)+sum(squeeze(Gij(2,2,:)).*drj.'.*phi_fv(jjz+nbeq,iinc));
            end
            k=k+1;
        end
        if sal.sxx==1
            szxdiff(ks,iinc)=sum(squeeze(Sij(1,1,:)).*drj.'.*phi_fv(jj,iinc));
            if ~isempty(jjz)
                szxdiff(ks,iinc)=szxdiff(ks,iinc)+sum(squeeze(Sij(1,2,:)).*drj.'.*phi_fv(jjz+nbeq,iinc));
            end
            ks=ks+1;
        end
        if sal.sxz==1
            szxdiff(ks,iinc)=sum(squeeze(Sij(2,1,:)).*drj.'.*phi_fv(jj,iinc));
            if ~isempty(jjz)
                szxdiff(ks,iinc)=szxdiff(ks,iinc)+sum(squeeze(Sij(2,2,:)).*drj.'.*phi_fv(jjz+nbeq,iinc));
            end
            ks=ks+1;
        end
        if sal.szz==1
            szxdiff(ks,iinc)=sum(squeeze(Sij(3,1,:)).*drj.'.*phi_fv(jj,iinc));
            if ~isempty(jjz)
                szxdiff(ks,iinc)=szxdiff(ks,iinc)+sum(squeeze(Sij(3,2,:)).*drj.'.*phi_fv(jjz+nbeq,iinc));
            end
            ks=ks+1;
        end
    end
    if sal.UPh==1
        if sal.Ux==1
            uxzdiff(k,iinc)=sum(squeeze(Gpr(1,1,:)).*drj.'.*phi_fv(jj,iinc)+squeeze(Gpr(1,2,:)).*drj.'.*phi_fv(jj+nbeq,iinc));k=k+1;
        end
        if sal.Uz==1
            uxzdiff(k,iinc)=sum(squeeze(Gpr(2,1,:)).*drj.'.*phi_fv(jj,iinc)+squeeze(Gpr(2,2,:)).*drj.'.*phi_fv(jj+nbeq,iinc));k=k+1;
        end
    end
    if sal.USh==1
        if sal.Ux==1
            uxzdiff(k,iinc)=sum(squeeze(Gsr(1,1,:)).*drj.'.*phi_fv(jj,iinc)+squeeze(Gsr(1,2,:)).*drj.'.*phi_fv(jj+nbeq,iinc));k=k+1;
        end
        if sal.Uz==1
            uxzdiff(k,iinc)=sum(squeeze(Gsr(2,1,:)).*drj.'.*phi_fv(jj,iinc)+squeeze(Gsr(2,2,:)).*drj.'.*phi_fv(jj+nbeq,iinc));k=k+1;
        end
    end
    if sal.UIh==1
        if sal.Ux==1
            uxzdiff(k,iinc)=sum(squeeze(Gi(1,1,:)).*drj.'.*phi_fv(jj,iinc)+squeeze(Gi(1,2,:)).*drj.'.*phi_fv(jj+nbeq,iinc));k=k+1;
        end
        if sal.Uz==1
            uxzdiff(k,iinc)=sum(squeeze(Gi(2,1,:)).*drj.'.*phi_fv(jj,iinc)+squeeze(Gi(2,2,:)).*drj.'.*phi_fv(jj+nbeq,iinc));k=k+1;
        end
    end
    if sal.UPt==1
        if sal.Ux==1
            uxzdiff(k,iinc)=sum(squeeze(Gpt(1,1,:)).*drj.'.*phi_fv(jj,iinc)+squeeze(Gpt(1,2,:)).*drj.'.*phi_fv(jj+nbeq,iinc));k=k+1;
        end
        if sal.Uz==1
            uxzdiff(k,iinc)=sum(squeeze(Gpt(2,1,:)).*drj.'.*phi_fv(jj,iinc)+squeeze(Gpt(2,2,:)).*drj.'.*phi_fv(jj+nbeq,iinc));k=k+1;
        end
    end
    if sal.USt==1
        if sal.Ux==1
            uxzdiff(k,iinc)=sum(squeeze(Gst(1,1,:)).*drj.'.*phi_fv(jj,iinc)+squeeze(Gst(1,2,:)).*drj.'.*phi_fv(jj+nbeq,iinc));k=k+1;
        end
        if sal.Uz==1
            uxzdiff(k,iinc)=sum(squeeze(Gst(2,1,:)).*drj.'.*phi_fv(jj,iinc)+squeeze(Gst(2,2,:)).*drj.'.*phi_fv(jj+nbeq,iinc));k=k+1;
        end
    end
end