function [Aijx,Aijz]=Aij_Gn_PSV(i,coord,para)

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
    
    if para.reg(mmat).bet==0
        vnxi= coord.vnx(i);
        vnzi= coord.vnz(i);
        drdn= (xij*vnxi+zij*vnzi)./rij;

        Gij(1,1,jjphi)      = 1i/4*kpi.*besselh(1,2,kpi.*rij).*drdn;
        %integracion gaussiana donde se requiere
        jjtmp2              = find((rij<=para.npplo*drj));
        Gij(1,1,jj(jjtmp2)) = T22_SH_r_small(coord,i,jr(jjtmp2),kpi,gaussian);

        % cuando rij=0
        signo1= (coord.Xim(i,m)==1) - (coord.Xim(i,m)==2);
        Gij(1,1,jj(rij==0))  = 0*Gij(1,1,jj(rij==0)) +signo1*0.5/drj(rij==0);
        
     
        
    else
        Gij(:,:,jj)	= Gij_PSV(ksi,kpi,rij,g,Ci,njj);
        %integracion gaussiana donde se requiere
        jjtmp2      = find((rij<=para.npplo*drj).*(rij>0));
        Gij(:,:,jj(jjtmp2)) = Gij_PSV_r_small(coord,x(i),z(i),jr(jjtmp2),ksi,kpi,gaussian,Ci);
        
        %cuando rij=0
        jjtmp2=find(rij==0);
        %     vn(1)=vnx(jr(jjtmp2));
        %     vn(2)=vnz(jr(jjtmp2));
        %     Gij(:,:,jj(jjtmp2)) = Gij_PSV_r_eq_0(ksi,kpi,0,drj(jjtmp2),Ci,vn);
        
        gn(1)=-vnz(jr(jjtmp2));
        gn(2)= vnx(jr(jjtmp2));
        Gij(:,:,jj(jjtmp2)) = Greenex_PSV(ksi,kpi,gn,Ci,drj(jjtmp2)); 
        if sum(coord.fl(i,:))==1
            %el solido esta en contacto con un fluido
            %se guarda los indices para proyectar sobre las normales
            jjf0=jj;
        end
    end
    
    %segun la normal
    for i0=1:2
        for j0=1:2
            Gij(i0,j0,jj) = ((m==im0) - (m~=im0))*squeeze(Gij(i0,j0,jj)).*drj.';
        end
    end
end

%Se reduce los Phiz a los que realmente existen
Aijx=[squeeze(Gij(1,1,:));squeeze(Gij(1,2,coord.logicphi))];%u_i=G_ij.Phi_j
Aijz=[squeeze(Gij(2,1,:));squeeze(Gij(2,2,coord.logicphi))];%u_i=G_ij.Phi_j

if sum(coord.fl(i,:))==1
    %uno de los medios es un fluido
    %proyection sobre normal y tangente para los componentes del solido
    %unicamente
    %cuando son dos fluidos o dos solidos no hay problema
    if exist('jjf0','var')
        
        %indices de los solidos:
        jjs = [jjf0,coord.nbeq+(1:sum(coord.logicphi))];
        un  = Aijx(jjs)*vnx(i)+Aijz(jjs)*vnz(i);%t_n=t_x*nx+t_z*nz
        ut  = Aijx(jjs)*vnz(i)-Aijz(jjs)*vnx(i);%t_t=t_x*nz-t_z*nx
        Aijx(jjs)= un;
        Aijz(jjs)= ut;
        %else
        %caso de un solo fluido
    end
end
