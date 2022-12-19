function [uy,tn,s]=campo_ref_SH_OPHeI_SE(xs,zs,coord,ks,kxk,kzsigno,fv,mu,para)
% funcion que da los campos de desplazamientos y las tractiones
% para una Onda Plana SH, homogenea e inhomogenea incidente en un semi espacio
% la onda incidente es unitaria con un componente horizontal cual es
% definido con respecto a ks (el numero de onda) por un factor kxk
% asi cuando abs(akx)<1, se reencuentra theta kx=akx*ks=sin(theta)*ks
% y   cuando abs(akx)>1, se puede incluir la ondas evanecentes
% la onda reflejada sera incluida
% xs,zs coordenadades del origen de la onda (fase nula)
% importante a tomar en cuenta cuando se incluye la atenuacion
% x,z posicion del receptor
% fv=1 : si el calculo se occupa por una incidencia en una fuente virtual se
% toma en cuenta la parte anticausal, si es por un receptor entonces la
% parte anticausal no se toma en cuenta

x  	= coord.x;
z  	= coord.z;

kx	= kxk*ks;
kz 	= kzsigno.*sqrt(ks^2-kx.^2); %kz = + o -, la 2 soluciones son correctas

kz	=-kz.*(sign(imag(kz))+(imag(kz)==0)).*((sign(zs-z)+(z==zs))*(1-fv)+1*fv).*(abs(kxk)>1)+kz.*(abs(kxk)<=1);

u_inc  = exp(-1i*(kx.*(x-xs)-kz.*(z-zs)));
if para.geo(1)==2 && para.cont(1).ruggeo==1 && para.fuenteimagen==1 || para.geo(1)==3 && para.nsubmed==1
    u_ref  = exp(-1i*(kx.*(x-xs)+kz.*(z+zs))).*(abs(kxk)<=1);
else
    u_ref  = 0;
end

% if fv==0
%     itmpinc= abs(u_inc)>1;
%     u_inc(itmpinc)=0;
%     itmpref= abs(u_ref)>1;
%     u_ref(itmpref)=0;
% end

uy  = u_inc+u_ref;

%sx	=mu*du/dx
sxy	=-mu*1i*kx.*uy;

%sz =mu*du/dz
szy	=-mu*1i*kz.*(-u_inc + u_ref);

nss = (para.sortie.sxy + para.sortie.syz);
s   = zeros(nss,coord.nbpt);

if isfield(coord,'vnx')
    vnx	= coord.vnx;
    vnz	= coord.vnz;
    tn  = sxy.*vnx+szy.*vnz;
else
    tn=[];
    if para.sortie.sxy==1
        s(1,:)= sxy;
    end
    if para.sortie.syz==1
        s(nss,:)= szy;
    end
end

end