function [Ht,Hr,Hi]=besselh_R_Im(ord,type,k,x,z)

if sign(imag(k))==-1 && type==1
    k=real(k)-1i*imag(k);
elseif sign(imag(k))==1 && type==2
    k=real(k)-1i*imag(k);
end

z   = sqrt(x^2+z^2);
x   = 0;
thr	= pi/2;

% thr    	= atan(z/x);

%integration lineaire en theta reel puis transpose au inhomogene puis k -
nth     = 1e3;

th      = linspace(90,0,nth)*pi/180;

%limite de l integration a chercher
%ref
kx      = 1.01*real(k);
kz      = sqrt(k^2-kx.^2);
kz      =-kz.*sign(imag(kz));
cOPref  = (-1)^(ord+type)*.5/pi/cos(ord*thr)*(((kz-1i*kx)/k).^ord+((-kz-1i*kx)/k).^ord).* ...
    exp(1i*(((-1)^(type+1))*kx*x-kz*abs(z)))./kz;
cOP     = cOPref;
while abs(imag(cOP))>abs(1e-2*imag(cOPref))
    kx  = kx*2;
    kz  = sqrt(k^2-kx.^2);
    cOP = (-1)^(ord+type)*.5/pi/cos(ord*thr)*(((kz-1i*kx)/k).^ord+((-kz-1i*kx)/k).^ord).* ...
        exp(1i*(((-1)^(type+1))*kx*x-kz*abs(z)))./kz;
end
krf     = real(k);
kxf     = kx;

kxr     = krf*cos(th);

kx      = kxr;
dkx     = zeros(1,nth);
dkx(1)  = (kx(2)	- kx(1)     )/2;
dkx(end)= (kx(nth)  - kx(nth-1) )/2;
dkx(2:end-1) = diff(kx(1:(nth-1)))/2 + diff(kx(2:(nth)))/2;

kz      = sqrt(k^2-kx.^2);
kz      =-kz.*sign(imag(kz));
kz(kz==0)=1e6;

cOP     = (-1)^(ord+type)*.5/pi/cos(ord*thr)*(((kz-1i*kx)/k).^ord+((-kz-1i*kx)/k).^ord).* ...
    exp(1i*(((-1)^(type+1))*kx*x-kz*abs(z)))./kz;
Hr      = 2*sum(cOP.*dkx);

nth2    = min(round(kxf/krf*nth),1e6);
th      = linspace(90,0,nth2)*pi/180;

kxi     = kxf-sin(th)*(kxf-krf);

kx      = kxi;
dkx     = zeros(1,nth2);
dkx(1)  = (kx(2)	- kx(1)     )/2;
dkx(end)= (kx(nth2) - kx(nth2-1))/2;
dkx(2:end-1) = diff(kx(1:(nth2-1)))/2 + diff(kx(2:(nth2)))/2;

kz      = sqrt(k^2-kx.^2);
kz      =-kz.*sign(imag(kz));
kz(kz==0)=1e6;

cOP     = (-1)^(ord+type)*.5/pi/cos(ord*thr)*(((kz-1i*kx)/k).^ord+((-kz-1i*kx)/k).^ord).* ...
    exp(1i*(((-1)^(type+1))*kx*x-kz*abs(z)))./kz;
Hi      = 2*sum(cOP.*dkx);

Ht      = Hr+Hi;