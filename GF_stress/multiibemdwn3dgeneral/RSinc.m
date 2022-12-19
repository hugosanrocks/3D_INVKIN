function [RSP,RSS,usi,upr,usr,kzP,kzS]=RSinc(C,kx,kp,ks)
kzP = sqrt(kp.^2-kx.^2);
kzS = sqrt(ks.^2-kx.^2);
if abs(kx)>abs(kp)
    kzP=-kzP.*sign(imag(kzP));
end
usi =sign(kx).*[kzS  kx] /ks;
upr =          [kx   kzP]/kp;
usr =sign(kx).*[kzS -kx] /ks;

if abs(kx)==abs(ks)
    RSP=0;
    RSS=0;
elseif kx~=0
    a(1,1)  = C(1,2)*kx.*upr(1)+C(1,1)*kzP.*upr(2);
    a(1,2)  = C(1,2)*kx.*usr(1)+C(1,1)*kzS.*usr(2);
    a(2,1)  = C(6,6)*(kzP.*upr(1)+kx.*upr(2));
    a(2,2)  = C(6,6)*(kzS.*usr(1)+kx.*usr(2));
    
    b(1,1)  =-(C(1,2)*kx.*usi(1)-C(1,1)*kzS.*usi(2));
    b(2,1)  =- C(6,6)*(-kzS.*usi(1)+kx.*usi(2));
    c       = a\b;
    
    RSP     = c(1);
    RSS     = c(2);
else
    RSP     = 0;
    RSS     = 1;
    usi = [kzS  kx] /ks;
    usr = [kzS -kx] /ks;
end
end