function [ F ] = Croissant( X,Y )
%Croissant (Moon valley) Sánchez-Semsa & Luzón 1995
a = 4.0;
H = 0.4/a;
XP= 0.7*a;											%JCMV
% YP=2.0;											%JCMV

R12=(X-a).^2+Y.*Y;
B2=XP*XP;
F=0.0;
if (R12 <= B2)
return
end
R22=X.*X+Y.*Y;
F=0.0;
if(R22>=a*a)
return
end
F=-(R12-B2).*(1.0-2.0*a*(a-X)./R12);
F = H*F;
end

