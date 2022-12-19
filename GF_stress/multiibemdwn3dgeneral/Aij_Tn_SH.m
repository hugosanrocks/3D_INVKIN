function Ai=Aij_Tn_SH(i,coord,para)

% coeficientes de la matriz A correspondiente a la continuidad de
% los esfuerzos normales
% i : indice de l element (de surface) sur lequel se calcul la contribution de toutes
% les autres sources virtuelles j (ponctuelle)

%%%%%%%%%%%%%
% unpacking %
%%%%%%%%%%%%%

x 	=coord.x;
z 	=coord.z;
vnxi=coord.vnx(i);
vnzi=coord.vnz(i);
dr	=coord.dr;
phi =coord.phi;

j       = 1:coord.nbpt;
jphi    = 1:coord.nbeq;
t22     = zeros(1,coord.nbeq);

gaussian   = para.gaussian;

%indice de los medios en contactos
im  = find(coord.Xim(i,:)~=0);

%min de los indice de los medios en contacto que se lleva el signo +
im0 = min(im);

for m=im
    if m==1 && para.geo(1)==3
        continue
    end
    mmat    = para.subm1(m);
    ksi     = para.reg(mmat).ksi;
    
    %indice (logic y natural) de los puntos de colocacion perteneciendo a m
    jjx     = coord.indm(m).ind;
    jr      = j(jjx);
    
    %indice (logic y natural) de los phi que hay que contemplar (columnas
    %de la matriz)
    jjphi   = false(coord.nbeq,1);
    jjphi(phi(j(jjx),m)) = true(1);
    jj      = jphi(jjphi);
    
    drj     = dr(jjx);
    xij     = x(i)-x(jjx);
    zij     = z(i)-z(jjx);
    rij     = sqrt(xij.^2+zij.^2);
    drdn    =(xij*vnxi+zij*vnzi)./rij;
    
    t22(jjphi)= 1i/4*ksi.*besselh(1,2,ksi.*rij).*drdn;
    
    %integracion gaussiana donde se requiere
    jjtmp2  = find((rij<=1*para.npplo*drj));
    t22(jj(jjtmp2)) = T22_SH_r_small(coord,i,jr(jjtmp2),ksi,gaussian);
    
    t22(jj) = t22(jj).*drj;
    
    % cuando rij=0 
    % en este punto, hay discontinuidad de las tractiones ti-ti'=phi
    % el signo1 es igual a +1 (resp. -1) si la normal apunta hacia el exterior
    % (resp. a interior) del medio a que pertenece el phi del punto de
    % colocacion en la posicion singular rij=0
    % como la convencion es que las normales esten siempre dirigidas hacia
    % abajo, hay que conocer si el phi pertenece a un contorno arriba o
    % abajo. Esta informacion esta en Xim
    signo1= (coord.Xim(i,m)==1) - (coord.Xim(i,m)==2);
    t22(jj(rij==0))  = 0*t22(jj(rij==0))+signo1*0.5;
    
    % en un punto de colocacion, hay siempre 2 medios en contacto,
    % aunque a veces este medio es el vacio
    % las ecuaciones de continuidad de traccion o de desplazamiento
    % involven entonces siempre 2 medios  : sigma(m)=sigma(m1) o u(m)=u(m1)
    % y como cada campo corresponde a (por ejemplo) u=u_diff+u_inc
    % se reorganisa entonces como :
    % sigma_diff(m)-sigma_diff(m1)=u_inc(m1)-u_inc(m)
    % los signos se fijan deacuerdo con el indice del medio,
    % + si el indice es el mas pequeno de los medios en contacto
    % - en lo contrario
    t22(jj)  = ((m==im0) - (m~=im0))*t22(jj);
end
Ai=t22;