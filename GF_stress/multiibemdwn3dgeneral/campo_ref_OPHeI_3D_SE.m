function [u,t]=campo_ref_OPHeI_3D_SE(xs,ys,zs,coord,kpi,ksi,kri,C,tipo_onda,kxk,kyk)
% Funcion que da los campos de desplazamientos y las tractiones
% para incidencias de Ondas Planas P,SV,SH, Homogenea e Inhomogenea
% Rayleigh en un semi-espacio. Las ondas de Love son para un medio de
% fondo estratificado. Asi se considera como campo incidente, la onda
% directa y las reflejadas por la frontera libre.
% La onda incidente es unitaria con un componente horizontal del numero de
% onda definido con respecto al numero de onda k por un factor kxk y kyk.
% xs,ys,zs coordenadades del origen de la onda (fase nula)
% importante a tomar en cuenta cuando se incluye la atenuacion
% Las posiciones de los receptores x,y,z estan en coord.
% El sistema de referencia es tal que x1=x y x2=z



n       = coord.nbpt;
x       = coord.x;%siempre se pide por lo menos las tractiones
y       = coord.y;
z       = coord.z;


if tipo_onda==1
    %incidencia de ondas P
    
    
    u_incP	= zeros(3,n);
    u_refP  = zeros(3,n);
    u_incS  = zeros(3,n);
    u_refS  = zeros(3,n);
    
    kx      = kxk*kpi;
    ky      = kyk*kpi;
    kr      = sqrt(kx^2+ky^2);
    %projection ds repere rz
    [RPP,RPS,upir,uprr,usrr,kzP,kzS]=RPinc(C,kr,kpi,ksi);
    %projection ds repere xyz
    upi(1)=upir(1)*kx/kr;
    upi(2)=upir(1)*ky/kr;
    upi(3)=upir(2);
    
    upr(1)=uprr(1)*kx/kr;
    upr(2)=uprr(1)*ky/kr;
    upr(3)=uprr(2);
    
    usr(1)=usrr(1)*kx/kr;
    usr(2)=usrr(1)*ky/kr;
    usr(3)=usrr(2);
    
    for i=1:3
        u_incP(i,:)  = upi(i)    *exp(-1i*(kx.*(x-xs)+ky.*(y-ys)-kzP.*(z-zs)));
        u_refP(i,:)  = upr(i)*RPP*exp(-1i*(kx.*(x-xs)+ky.*(y-ys)+kzP.*(z-zs)));
        u_refS(i,:)  = usr(i)*RPS*exp(-1i*(kx.*(x-xs)+ky.*(y-ys)+kzS.*(z-zs)));
    end
    u=u_incP+u_incS+u_refP+u_refS;
elseif tipo_onda==2
    %incidencia de ondas SV
    u_incS  = zeros(3,n);
    u_refS  = zeros(3,n);
    u_incP  = zeros(3,n);
    u_refP  = zeros(3,n);
    
    kx      = kxk*kpi;
    ky      = kyk*kpi;
    kr      = sqrt(kx^2+ky^2);
    %projection ds repere rz
    [RSP,RSS,usir,uprr,usrr,kzP,kzS]=RSinc(C,kr,kpi,ksi);
    %projection ds repere xyz
    usi(1)=usir(1)*kx/kr;
    usi(2)=usir(1)*ky/kr;
    usi(3)=usir(2);
    
    upr(1)=uprr(1)*kx/kr;
    upr(2)=uprr(1)*ky/kr;
    upr(3)=uprr(2);
    
    usr(1)=usrr(1)*kx/kr;
    usr(2)=usrr(1)*ky/kr;
    usr(3)=usrr(2);
    
    for i=1:3
        u_incS(i,:)  = usi(i)    *exp(-1i*(kx.*(x-xs)+ky.*(y-ys)-kzS.*(z-zs)));
        u_refP(i,:)  = upr(i)*RSP*exp(-1i*(kx.*(x-xs)+ky.*(y-ys)+kzP.*(z-zs)));
        u_refS(i,:)  = usr(i)*RSS*exp(-1i*(kx.*(x-xs)+ky.*(y-ys)+kzS.*(z-zs)));
    end
    u=u_incP+u_incS+u_refP+u_refS;
elseif tipo_onda==3
    %incidencia de ondas SH
    u_incS  = zeros(3,n);
    u_refS  = zeros(3,n);
    u_incP  = zeros(3,n);
    u_refP  = zeros(3,n);
    
    kx      = kxk*kpi;
    ky      = kyk*kpi;
    kr      = sqrt(kx^2+ky^2);
    kzS     = sqrt(ksi.^2-kr.^2);

    usi     = [ky -kx 0]/kr;
    usr     = usi;
    
    for i=1:3
        u_incS(i,:)  = usi(i)*exp(-1i*(kx.*(x-xs)+ky.*(y-ys)-kzS.*(z-zs)));
        u_refS(i,:)  = usr(i)*exp(-1i*(kx.*(x-xs)+ky.*(y-ys)+kzS.*(z-zs)));
    end
    u=u_incP+u_incS+u_refP+u_refS;
elseif tipo_onda==4
    %incidencia de ondas Rayleigh en SE
    
    chi1    = sqrt(1-(kpi/kri)^2);
    chi2    = sqrt(1-(ksi/kri)^2);
    chi1    = sign(real(chi1)).*chi1;
    chi2    = sign(real(chi2)).*chi2;
    
    sq12    = sqrt(chi1*chi2);
    sq21    = sqrt(chi1/chi2);
    
    krix    = sign(kxk)*kri;
    kriy    = sign(kyk)*kri;
    krri    = sqrt(krix^2+kryi^2);
    kriz    = kri;
    
    u       = zeros(3,n);
    
    u(1,:)  =         (       exp(-kriz*chi1*abs(z))-sq12.*exp(-kriz*chi2*abs(z))).*exp(-1i*(krix.*(x-xs)+kriy.*(y-ys)))*kxk/krri;
    u(2,:)  =         (       exp(-kriz*chi1*abs(z))-sq12.*exp(-kriz*chi2*abs(z))).*exp(-1i*(krix.*(x-xs)+kriy.*(y-ys)))*kyk/krri;
    u(3,:)  = 1i*sq21*(-sq12.*exp(-kriz*chi1*abs(z))+      exp(-kriz*chi2*abs(z))).*exp(-1i*(krix.*(x-xs)+kriy.*(y-ys)));
elseif tipo_onda==5
    errordlg('no permitido en un semi-espacio');
end

if tipo_onda==1 || tipo_onda==2
    du1dx1=-1i*kx *u(1,:); %du1/dx1
    du2dx1=-1i*kx *u(2,:);
    du3dx1=-1i*kx *u(3,:);
    du1dx2=-1i*ky *u(1,:); %du1/dx2
    du2dx2=-1i*ky *u(2,:);
    du3dx2=-1i*ky *u(3,:);
    du1dx3= 1i*kzP*u_incP(1,:) + 1i*kzS*u_incS(1,:) -1i*kzP*u_refP(1,:) -1i*kzS*u_refS(1,:);
    du2dx3= 1i*kzP*u_incP(2,:) + 1i*kzS*u_incS(2,:) -1i*kzP*u_refP(2,:) -1i*kzS*u_refS(2,:);
    du3dx3= 1i*kzP*u_incP(3,:) + 1i*kzS*u_incS(3,:) -1i*kzP*u_refP(3,:) -1i*kzS*u_refS(3,:);
else
    du1dx1=-1i*krix*u(1,:); %du1/dx1
    du2dx1=-1i*krix*u(2,:);
    du3dx1=-1i*krix*u(3,:);
    du1dx2=-1i*kriy*u(1,:); %du1/dx2
    du2dx2=-1i*kriy*u(2,:);
    du3dx2=-1i*kriy*u(3,:);
    du1dx3= kriz*          (   -chi1*exp(-kriz*chi1*abs(z))+chi2*sq12.*exp(-kriz*chi2*abs(z))).*exp(-1i*(krix.*(x-xs)+kriy.*(y-ys)))*kxk/krri; %du1/dx3
    du2dx3= kriz*          (   -chi1*exp(-kriz*chi1*abs(z))+chi2*sq12.*exp(-kriz*chi2*abs(z))).*exp(-1i*(krix.*(x-xs)+kriy.*(y-ys)))*kyk/krri;
    du3dx3= kriz*1i*sq21*(sq12.*chi1*exp(-kriz*chi1*abs(z))-chi2*      exp(-kriz*chi2*abs(z))).*exp(-1i*(krix.*(x-xs)+kriy.*(y-ys)));
end


if isfield(coord,'vnx')
    vn(1,:)	= coord.vnx;
    vn(2,:)	= coord.vny;
    vn(3,:)	= coord.vnz;
    
    s       = zeros(3,3,n);
    % s(1,1,:)= C(1,1)*du1dx1+C(1,2)*du2dx2+C(1,3)*du3dx3;
    % s(2,2,:)= C(2,1)*du1dx1+C(2,2)*du2dx2+C(2,3)*du3dx3;
    % s(3,3,:)= C(3,1)*du1dx1+C(3,2)*du2dx2+C(3,3)*du3dx3;
    % s(1,2,:)= C(6,6)*(du1dx2+du2dx1);
    % s(1,3,:)= C(5,5)*(du1dx3+du3dx1);
    % s(2,3,:)= C(4,4)*(du2dx3+du3dx2);
    s(1,1,:)= C(1,1)*du1dx1+C(1,2)*du2dx2+C(1,2)*du3dx3;
    s(2,2,:)= C(1,2)*du1dx1+C(1,1)*du2dx2+C(1,2)*du3dx3;
    s(3,3,:)= C(1,2)*du1dx1+C(1,2)*du2dx2+C(1,1)*du3dx3;
    s(1,2,:)= C(6,6)*(du1dx2+du2dx1);
    s(1,3,:)= C(6,6)*(du1dx3+du3dx1);
    s(2,3,:)= C(6,6)*(du2dx3+du3dx2);

    t       = zeros(3,n);
    t(1,:)  = squeeze(s(1,1,:)).'.*vn(1,:)+squeeze(s(1,2,:)).'.*vn(2,:)+squeeze(s(1,3,:)).'.*vn(3,:);
    t(2,:)  = squeeze(s(1,2,:)).'.*vn(1,:)+squeeze(s(2,2,:)).'.*vn(2,:)+squeeze(s(2,3,:)).'.*vn(3,:);
    t(3,:)  = squeeze(s(1,3,:)).'.*vn(1,:)+squeeze(s(2,3,:)).'.*vn(2,:)+squeeze(s(3,3,:)).'.*vn(3,:);
    
else
    t=[];
end
end

