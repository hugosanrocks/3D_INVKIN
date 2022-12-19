function [G220,G22z0,G22h,G22zh]=G22_SH_DWN(htop,hbot,zs,k1f,MAT,ics,ncapas)

%%%%%%%%%%%%%%%%%
% interface sup %
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des deplacements %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
zt      = htop-zs;%x=z dans l epaisseur
exp2t   = exp(-1i.*k1f.*abs(zt));
fac     = 1/(4*pi*MAT(ics).Ci(6,6));
G220    =-1i*fac./k1f.*exp2t;%pair en k2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des derivees utiles au calcul des contraintes normales %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% derivee p/x1 en X (x1=z)
%pair
G22z0=-sign(zt)*fac*exp2t;

%%%%%%%%%%%%%%%%%
% interface inf %
%%%%%%%%%%%%%%%%%
if ics~=ncapas
    zb      = hbot-zs;%x=z dans l epaisseur
    exp2b   = exp(-1i.*k1f.*abs(zb));
    
    %pair en k2
    G22h    =-1i*fac./k1f.*exp2b;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcul des derivees utiles au calcul des contraintes normales %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % derivee p/x1 en X (x1=z)
    %pair
    G22zh   =-sign(zb)*fac*exp2b;
else
    G22h    = 0;
    G22zh   = 0;
end