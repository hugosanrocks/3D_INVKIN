function Hr=besselh_R_Im2(ord,type,k,x,z,dr)
% [Ht,Hr,Hi]=besselh_R_Im2(ord,type,k,x,z)
r       = sqrt(x^2+z^2);

if sign(imag(k))==-1 && type==1
    k=real(k)-1i*imag(k);
elseif sign(imag(k))==1 && type==2
    k=real(k)-1i*imag(k);
end

kr      = real(k);
Hr      = 2*integral(@(kx) Comp_hankel_OP(kx,ord,type,k,r,dr),0 ,kr,'RelTol',1e-3);

% %limite de l integration a chercher
% %ref
% kx      = 1.01*kr;
% cOPref  = Comp_hankel_OP(kx,ord,type,k,r);
% cOP     = cOPref;
% while abs((cOP))>abs(1e-3*(cOPref))
%     kx  = kx*5;
%     cOP = Comp_hankel_OP(kx,ord,type,k,r);
% end
% kxf     = kx;
% 
% Hi      = 2*integral(@(kx) Comp_hankel_OP(kx,ord,type,k,r),kr,kxf,'RelTol',1e-5);
% 
% Ht      = Hr+Hi;%
