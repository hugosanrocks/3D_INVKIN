function Ai=Aij_Gn_SH(i,coord,para)

% coeficientes de la matriz A correspondiente a la continuidad de los desplazamientos
% i : indice de l element (de surface) sur lequel se calcul la contribution de toutes
% les autres sources virtuelles j (ponctuelle)

%%%%%%%%%%%%%
% unpacking %
%%%%%%%%%%%%%

x 	= coord.x;
z 	= coord.z;
dr	= coord.dr;
phi = coord.phi;

j       = 1:coord.nbpt;
jphi    = 1:coord.nbeq;
gn22    = zeros(1,coord.nbeq);

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
    ksi = para.reg(mmat).ksi;
    mu	= para.reg(mmat).Ci(6,6);
    
    %indice (logic y natural) de los puntos de colocacion perteneciendo a m
    jjx	= coord.indm(m).ind;
    jr  = j(jjx);
    
    %indice (logic y natural) de los phi que hay que contemplar (columnas
    %de la matriz)
    jjphi   = false(coord.nbeq,1);
    jjphi(phi(j(jjx),m)) = true(1);
    jj      = jphi(jjphi);
    
    drj = dr(jjx);
    xij = x(i)-x(jjx);
    zij = z(i)-z(jjx);
    rij = sqrt(xij.^2+zij.^2);
    
    gn22(jj)= G22_SH(ksi,rij,mu);
    
    %integracion gaussiana donde se requiere
    jjtmp2=find((rij<=1*para.npplo*drj).*(rij>0));
    gn22(jj(jjtmp2)) = G22_SH_r_small(coord,x(i),z(i),jr(jjtmp2),ksi,gaussian,mu);
    
    %cuando rij=0
    gn22(jj(rij==0)) = G22_SH_r_eq_0(ksi,drj(rij==0),mu);
    
    %este signo solo es por el orden en las ecuaciones con respeto al
    %termino fuente, el medio de menor indice entre los en contactos se
    %lleva el signo +
    gn22(jj) = ((m==im0) - (m~=im0))*gn22(jj).*drj;
end
Ai=gn22;
