function G=Greenex_PSV(ks,kp,g,C,dr)
%retorna la expression analytica de la integracion de la funcion de Green 
%en el elemento centrado en la fuente

geu   =0.57721566;

i2pi=1i*2/pi;
fac1= 1-i2pi*(geu-1);
fac2= 1-i2pi*(geu-1-1/3);
fac3= 1-i2pi*(geu-3/4-1/3);

logp= log(kp*dr/4);
logs= log(ks*dr/4);

h0P = (fac1-i2pi*logp) - dr.^2*kp.^2/48*(fac2-i2pi*logp);
h0S = (fac1-i2pi*logs) - dr.^2*ks.^2/48*(fac2-i2pi*logs);

h2P = 1i/pi + dr.^2*kp.^2/96*(fac3-i2pi*logp);
h2S = 1i/pi + dr.^2*ks.^2/96*(fac3-i2pi*logs);

A   = h0P/C(1,1)+h0S/C(6,6);
B   = h2P/C(1,1)-h2S/C(6,6);

G   = zeros(2,2);
d   = eye(2);

i=1;
for j=1:2
    G(i,j)=1/(1i*8)*(d(i,j)*A-(2*g(i).*g(j)-d(i,j)).*B);
end

i=2;j=2;
G(i,j)=1/(1i*8)*(d(i,j)*A-(2*g(i).*g(j)-d(i,j)).*B);
G(2,1)=G(1,2);