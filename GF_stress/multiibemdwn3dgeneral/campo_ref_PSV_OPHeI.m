function [u,t,sout]=campo_ref_PSV_OPHeI(xs,zs,coord,kpi,ksi,kri,C,tipo_onda,kxk)
% funcion que da los campos de desplazamientos y las tractiones
% para Ondas Planas P-SV, homogenea e inhomogenea incidente en un espacio
% completo
% la onda incidente es unitaria con un componente horizontal de la
% polarisacion cual es definido con respecto a k (el numero de onda) por un
% factor kxk
% asi cuando abs(akx)<1, se reencuentra theta kx=akx*k=sin(theta)*k
% y   cuando abs(akx)>1, se puede incluir la ondas evanecentes
% xs,zs coordenadades del origen de la onda (fase nula)
% importante a tomar en cuenta cuando se incluye la atenuacion
% x,z posicion del receptor
% fv=1 : si el calculo se occupa por una incidencia en una fuente virtual se
% toma en cuenta la parte anticausal, si es por un receptor entonces la
% parte anticausal no se toma en cuenta
% x1=x
% x2=z

% Funcion que da los campos de desplazamientos y las tractiones
% para incidencias de Ondas Planas P,SV,SH, Homogenea e Inhomogenea
% Rayleigh en un espacio completo. Las ondas de Love son para un medio de
% fondo estratificado.
% La onda incidente es unitaria con un componente horizontal del numero de
% onda definido con respecto al numero de onda k por un factor kxk y kyk.
% xs,ys,zs coordenadades del origen de la onda (fase nula)
% importante a tomar en cuenta cuando se incluye la atenuacion
% Las posiciones de los receptores x,y,z estan en coord.
% El sistema de referencia es tal que x1=x y x2=z


n       = coord.nbpt;
x       = coord.x;%siempre se pide por lo menos las tractiones
z       = coord.z;


if tipo_onda==1
    %incidencia de ondas P
    
    u_incS  = zeros(2,n);
    u_incP  = zeros(2,n);
    
    kx      = kxk*kpi;
    [~,~,upi,~,~,kzP,kzS]=RPinc(C,kx,kpi,ksi);
    
    for i=1:2
        u_incP(i,:)  = upi(i)*    exp(-1i*(kx.*(x-xs)-kzP.*(z-zs)));
    end
    u=u_incP;
elseif tipo_onda==2
    %incidencia de ondas S
    u_incS  = zeros(2,n);
    u_incP  = zeros(2,n);
    
    kx      = kxk*ksi;
    [~,~,usi,~,~,kzP,kzS]=RSinc(C,kx,kpi,ksi);
    for i=1:2
        u_incS(i,:)  = usi(i)*    exp(-1i*(kx.*(x-xs)-kzS.*(z-zs)));
    end
    u=u_incS;
else
    %incidencia de ondas Rayleigh
    
    chi1    = sqrt(1-(kpi/kri)^2);
    chi2    = sqrt(1-(ksi/kri)^2);
    chi1    = sign(real(chi1)).*chi1;
    chi2    = sign(real(chi2)).*chi2;
    
    sq12    = sqrt(chi1*chi2);
    sq21    = sqrt(chi1/chi2);
    
    krix    = sign(kxk)*kri;
    kriz    = kri;
    
    u       = zeros(2,n);
    
    u(1,:)  =         (       exp(-kriz*chi1*abs(z))-sq12.*exp(-kriz*chi2*abs(z))).*exp(-1i*(krix.*(x-xs)))*sign(kxk);
    u(2,:)  = 1i*sq21*(-sq12.*exp(-kriz*chi1*abs(z))+      exp(-kriz*chi2*abs(z))).*exp(-1i*(krix.*(x-xs)));
end

if tipo_onda==1 || tipo_onda==2
    du1dx1= (-1i*kx )*u(1,:); %du1/dx1
    du2dx1= (-1i*kx )*u(2,:);
    du1dx2= ( 1i*kzP)*u_incP(1,:) + ( 1i*kzS)*u_incS(1,:);
    du2dx2= ( 1i*kzP)*u_incP(2,:) + ( 1i*kzS)*u_incS(2,:);
else
    du1dx1= (-1i*krix )*u(1,:); %du1/dx1
    du2dx1= (-1i*krix )*u(2,:);
    du1dx2= kriz*          (   -chi1*exp(-kriz*chi1*abs(z))+chi2*sq12.*exp(-kriz*chi2*abs(z))).*exp(-1i*(krix.*(x-xs)))*sign(kxk);
    du2dx2= kriz*1i*sq21*(sq12.*chi1*exp(-kriz*chi1*abs(z))-chi2*      exp(-kriz*chi2*abs(z))).*exp(-1i*(krix.*(x-xs)));
end
s = zeros(2,2,n);
sout = zeros(3,n);
s(1,1,:)= C(1,1)*du1dx1+C(1,2)*du2dx2;
s(2,2,:)= C(1,1)*du2dx2+C(1,2)*du1dx1;%C22*du2/dx2+C12*du1/dx1
s(1,2,:)= C(6,6)*(du1dx2+du2dx1);
sout(1,:) = s(1,1,:); %sxx
sout(2,:) = s(1,2,:); %sxz
sout(3,:) = s(2,2,:); %szz
    
if isfield(coord,'vnx')
    vn(1,:)	= coord.vnx;
    vn(2,:)	= coord.vnz;
    t       = zeros(2,n);
    t(1,:)  = squeeze(s(1,1,:)).'.*vn(1,:)+squeeze(s(1,2,:)).'.*vn(2,:);
    t(2,:)  = squeeze(s(1,2,:)).'.*vn(1,:)+squeeze(s(2,2,:)).'.*vn(2,:);
else
    t=[];
end

end