delta   =0* 0.001;
w       = 2-1i*delta;
bet     = 1;
k       = 1e2*(7.8540e-04 - 3.9270e-06i);%w/bet*(1-1i/(2*1000));
x       = 0;
z       = sqrt(13);
thr      = atan(z/x);
r       = sqrt( x^2 + z^2 );
kr      = k*r;


ref(1)  = besselh(0,2,real(kr));
txt{1}  = 'H0              ';
ref(2)  = besselh(0,2,(kr));
txt{2}  = 'H0 w cmplx      ';

ref2(1)  = besselh(2,2,real(kr));
txt2{1}  = 'H2              ';
ref2(2)  = besselh(2,2,(kr));
txt2{2}  = 'H2 w cmplx      ';

%integration lineaire en dkx
nkx     = 4*1e5;
kx      = linspace(-10*real(k),10*real(k),nkx);
dkx     = kx(2)-kx(1);
kz      = sqrt(k^2-kx.^2);

tmp     = exp(-1i*(kx*x+kz*abs(z)))./kz;
ref(3)  = 1/pi*sum(tmp)*dkx;
txt{3}  = 'int OP dkx lin  ';
figure(10);hold on
plot(kx,abs(tmp),'r')

tmp2     = 1/cos(2*thr)*(kz.^2-kx.^2)./k.^2.*exp(-1i*(kx*x+kz*abs(z)))./kz;
ref2(3)  = 1/pi*sum(tmp2)*dkx;
txt2{3}  = 'int OP dkx lin  ';
plot(kx,imag(tmp2),'m')


%integration lineaire en theta reel puis transpose au inhomogene
nth     = 1000;
th      = linspace(90,0,nth+1)*pi/180;
th      = th(1:nth);
kx      = real(k)*cos(th);
kxf     = 30*real(k);
kx      = [kx,  kxf-sin(th)*(kxf-real(k));];
% th      = linspace(180,90,nth)*pi/180;
kx      = [-fliplr(kx), kx];

kz      = sqrt(k^2-kx.^2);
dkx     = zeros(1,4*nth);
dkx(1)  = kx(2)     - kx(1);
dkx(end)= kx(4*nth)   - kx(4*nth-1);
dkx(2:end-1) = diff(kx(1:(4*nth-1)))/2 + diff(kx(2:4*nth))/2;

kz      =-kz.*sign(imag(kz));

tmp     = exp(-1i*(kx*x+kz*abs(z)))./kz;
ref(4)  = 1/pi*sum(tmp.*dkx);
txt{4}  = 'int OP dkx theta';
plot(kx,abs(tmp),'k.')


tmp2    = 1/cos(2*thr)*(kz.^2-kx.^2)./k.^2.*exp(-1i*(kx*x+kz*abs(z)))./kz;
ord     = 2;
tmp2    = .5/cos(ord*thr)*(((kz-1i*kx)/k).^ord+((-kz-1i*kx)/k).^ord).*exp(-1i*(kx*x+kz*abs(z)))./kz;

ref2(4)  = 1/pi*sum(tmp2.*dkx);
txt2{4}  = 'int OP dkx theta';
plot(kx,abs(tmp2),'.')


for i=1:4
    txt0{i}=[num2str(ref(i),'%#5.3s\t'),'   ',txt{i},'   ',num2str(ref2(i),'%s\t'),'    ',txt2{i}];
end

txt0.'