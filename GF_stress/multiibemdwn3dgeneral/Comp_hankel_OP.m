function y=Comp_hankel_OP(kx,ord,type,k,z,dr)
% calcula cada componente de las funciones de hankel en la descomposicion
% de helmoltz
% kx  = componente horizontal del numero de onda
% ord = orden de la funcion de hankel H_ord^type (nu)
% type= type de la funcion de hankel             (K )
% k   = numero de onda
% z   = distancia r porque x = 0

kz      = sqrt(k^2-kx.^2);
kz      =-kz.*(sign(imag(kz))-(sign(imag(kz))==0));
kz(kz==0)=1e6;

% y=(-1)^(ord+type)*.5/pi/cos(ord*thr)*(((kz-1i*kx)/k).^ord+((-kz-1i*kx)/k).^ord).* ...
%     exp(1i*(((-1)^(type+1))*kx*x-kz*abs(z)))./kz;

y=(-1)^(ord+type)*.5/pi/cos(ord*pi/2)*(((kz-1i*kx)/k).^ord+((-kz-1i*kx)/k).^ord).* ...
    exp(-1i*kz*abs(z))./kz.*sinc(-dr*kx/2/pi);