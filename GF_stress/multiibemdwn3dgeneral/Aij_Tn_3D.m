function [Aijx,Aijy,Aijz]=Aij_Tn_3D(i,coord,para)

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

gaussian   = para.gaussian;

%indice de los medios en contactos
im  = find(coord.Xim(i,:)~=0);

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
    
    vn(1,1)	= vnx(i);
    vn(2,1)	= vny(i);
    vn(3,1)	= vnz(i);
    
    tij(:,:,jjphi)=Tij_3D(ksi,kpi,rij,g,Ci,vn);
    
    %integracion gaussiana donde se requiere
    if para.dim == 4 % 3Dgeneral
    jjtmp2  = find((rij<=1.5*drj).*(rij~=0)); % a menos de 1.5 radios
    if i == 4
      disp(i)
    end
    tij(:,:,jj(jjtmp2)) = Tij_3D_r_smallGEN(coord,i,jr(jjtmp2),ksi,kpi,para,Ci);
    else %3D axisimétrico
    jjtmp2  = find((rij<para.npplo/2*drj).*(rij~=0));
    tij(:,:,jj(jjtmp2)) = Tij_3D_r_small(coord,i,jr(jjtmp2),ksi,kpi,gaussian,Ci);
    end
    for i0=1:3
        for j0=1:3
            tij(i0,j0,jj) = squeeze(tij(i0,j0,jj)).*dAj.';
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
    % involven entonces siempre 2 medios  : sigma(m)=sigma(m1) o u(m)=u(m1)
    % y como cada campo corresponde a (por ejemplo) u=u_diff+u_inc
    % se reorganisa entonces como :
    % sigma_diff(m)-sigma_diff(m1)=u_inc(m1)-u_inc(m)
    % los signos se fijan deacuerdo con el indice del medio,
    % + si el indice es el mas pequeno de los medios en contacto
    % - en lo contrario
    tij(:,:,jj) = ((m==im0) - (m~=im0))*tij(:,:,jj);
end
Aijx=[squeeze(tij(1,1,:));squeeze(tij(1,2,:));squeeze(tij(1,3,:))];%t_x=t_xj.Phi_j
Aijy=[squeeze(tij(2,1,:));squeeze(tij(2,2,:));squeeze(tij(2,3,:))];%t_y=t_yj.Phi_j
Aijz=[squeeze(tij(3,1,:));squeeze(tij(3,2,:));squeeze(tij(3,3,:))];%t_z=t_zj.Phi_j