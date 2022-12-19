function [RPP,RPS,upi,upr,usr,kzP,kzS]=RPinc(C,kx,kp,ks)
kzP = sqrt(kp.^2-kx.^2);
kzS = sqrt(ks.^2-kx.^2);

upi = [kx -kzP]/kp;%incidencia hacia z<0
upr = [kx  kzP]/kp;
usr =sign(kx).*[kzS -kx]/ks;
if abs(kx)==abs(kp)
    RPP=0;
    RPS=0;
elseif kx~=0
    a(1,1)  = C(1,2)*kx.*upr(1)+C(1,1)*kzP.*upr(2);
    a(1,2)  = C(1,2)*kx.*usr(1)+C(1,1)*kzS.*usr(2);
    a(2,1)  = C(6,6)*(kzP.*upr(1)+kx.*upr(2));
    a(2,2)  = C(6,6)*(kzS.*usr(1)+kx.*usr(2));
    
    b(1,1)  =-(C(1,2)*kx.*upi(1)-C(1,1)*kzP.*upi(2));
    b(2,1)  =- C(6,6)*(-kzP.*upi(1)+kx.*upi(2));
    c       = a\b;
    
    RPP     = c(1);
    RPS     = c(2);
else
    RPP     =-1;
    RPS     = 0;
    usr     = [kzS -kx]/ks;
    
end

end
