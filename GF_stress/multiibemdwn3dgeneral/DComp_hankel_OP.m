function y=DComp_hankel_OP(kx,ord,type,k,x,z,dr)
% calcula cada componente de las funciones de hankel en la descomposicion
% de helmoltz
% kx  = componente horizontal del numero de onda
% ord = orden de la funcion de hankel H_ord^type (nu)
% type= type de la funcion de hankel             (K )
% k   = numero de onda
%(x,z)= coordenadas del receptor relativas a la fuente

kz      = sqrt(k^2-kx.^2);

kz(kz==0)=1e6;

thr    	= atan(z/x);
% if abs(mod(ord*thr,pi)-pi/2)<1e-3
%     x   = x*(1+1e-6);
%     z   = z*(1-1e-6);
%     thr	= atan(z/x);
% end

y       =(-1)^(ord+type)*.5/pi/cos(ord*thr)*(((kz-1i*kx)/k).^ord+((-kz-1i*kx)/k).^ord).* ...
    exp(1i*(((-1)^(type+1))*kx*x-kz*abs(z)))./kz.*sinc(-dr*kx/2/pi);
% toto=0;

% y=(-1)^(ord+type)*.5/pi/cos(ord*pi/2)*(((kz-1i*kx)/k).^ord+((-kz-1i*kx)/k).^ord).* ...
%     exp(-1i*kz*abs(z))./kz;

%H2(2)
% y       =1i/(k^2*z*pi)*
% (-1)^(ord+type)*.5/pi/cos(ord*thr)*(((kz-1i*kx)/k).^ord+((-kz-1i*kx)/k).^ord).* ...
%     exp(1i*(((-1)^(type+1))*kx*x-kz*abs(z)))./kz;