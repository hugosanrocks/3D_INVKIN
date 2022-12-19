function [l,dr]=malla_contorno_fast(nseg,cont)
% funcion que permite posicionar los puntos approximadamente separados de
% dr0, los puntos se situan en medio de cada segmento

% xc  = zeros(nseg,1);
% zc  = zeros(nseg,1);
% vnx = zeros(nseg,1);
% vnz = zeros(nseg,1);
r   = zeros(nseg,1);
dr  = zeros(nseg,1);
l   = zeros(nseg,1);

dr0 = cont.long/nseg;
rc  = dr0*(2*(1:nseg)-1)/2; %posicion deseada

j=1;
for i=1:nseg
    [~,imin]= min(abs(rc(i)-cont.vec.r(j:end)));
    %     xc(i)   = cont.x(imin);
    %     zc(i)   = cont.z(imin);
    %     vnx(i)  = cont.vnx(imin);
    %     vnz(i)  = cont.vnz(imin);
    imin=j-1+imin;
    r(i)    = cont.vec.r(imin);
    l(i)    = imin;
    j=imin;
end

i=2:(nseg-1);
dr(i)=(r(i+1)-r(i-1))/2;
dr(1)=2*r(1);
dr(nseg)=2*(cont.long-r(nseg));