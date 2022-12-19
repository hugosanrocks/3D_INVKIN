function [Aijx,Aijz]=Aij_Tn_PSV(i,coord,para)

% coeficientes de la matriz A correspondiente a la continuidad de
% los esfuerzos normales
% i : indice de l element (de surface) sur lequel se calcul la contribution de toutes
% les autres sources virtuelles j (ponctuelle)


%%%%%%%%%%%%%
% unpacking %
%%%%%%%%%%%%%

x 	=coord.x;
z 	=coord.z;
vnx	=coord.vnx;
vnz	=coord.vnz;
dr	=coord.dr;
phi =coord.phi;

j       = 1:coord.nbpt;
jphi    = 1:coord.nbeq;
tij     = zeros(2,2,coord.nbeq);

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
    
    vn(1,1)	= vnx(i);
    vn(2,1)	= vnz(i);
    
    if para.reg(mmat).bet==0
        tij(1,1,jjphi)  =-para.reg(mmat).rho*(2*pi*para.fj)^2*G22_SH(kpi,rij,1);
        %integracion gaussiana donde se requiere
        jjtmp2          = find((rij<=1*para.npplo*drj).*(rij>0));
        tij(1,1,jj(jjtmp2)) =-para.reg(mmat).rho*(2*pi*para.fj)^2*...
            G22_SH_r_small(coord,x(i),z(i),jr(jjtmp2),kpi,gaussian,1);
        
        %cuando rij=0
        tij(1,1,jj(rij==0)) =-para.reg(mmat).rho*(2*pi*para.fj)^2*G22_SH_r_eq_0(kpi,drj(rij==0),1);
        
        tij(1,1,jj) = squeeze(tij(1,1,jj)).*drj.';
    else
        tij(:,:,jjphi)=Tij_PSV(ksi,kpi,rij,g,Ci,vn);
        %integracion gaussiana donde se requiere
        jjtmp2=find((rij<=para.npplo*drj));
        tij(:,:,jj(jjtmp2)) = Tij_PSV_r_small(coord,i,jr(jjtmp2),ksi,kpi,gaussian,Ci);
        for i0=1:2
            for j0=1:2
                tij(i0,j0,jj) = squeeze(tij(i0,j0,jj)).*drj.';
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
        tij(:,:,jj(rij==0))  = 0*tij(:,:,jj(rij==0)) +signo1*0.5*eye(2);
        
        if sum(coord.fl(i,:))==1
            %el solido esta en contacto con un fluido
            %se guarda los indices para proyectar sobre las normales
            jjf0=jj;
        end
        
    end
    
    
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


Aijx=[squeeze(tij(1,1,:));squeeze(tij(1,2,coord.logicphi))];%t_x=t_xj.Phi_j
Aijz=[squeeze(tij(2,1,:));squeeze(tij(2,2,coord.logicphi))];%t_z=t_zj.Phi_j

if sum(coord.fl(i,:))==1
    %uno de los medios es un fluido
    %proyection sobre normal y tangente para los componentes del solido
    %unicamente
    %cuando son dos fluidos o dos solidos no hay problema
    if exist('jjf0','var')
        %indices de los solidos:
        jjs = [jjf0,coord.nbeq+(1:sum(coord.logicphi))];
        tn  = Aijx(jjs)*vnx(i)+Aijz(jjs)*vnz(i);%t_n=t_x*nx+t_z*nz
        tt  = Aijx(jjs)*vnz(i)-Aijz(jjs)*vnx(i);%t_t=t_x*nz-t_z*nx
        Aijx(jjs)= tn;
        Aijz(jjs)= tt;
        %else
        %caso de un solo fluido
    end
end