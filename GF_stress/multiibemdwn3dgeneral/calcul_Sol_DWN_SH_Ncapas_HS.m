function Solution=calcul_Sol_DWN_SH_Ncapas_HS(para,zs,DWN)
% METHODE : @ recuperacion de k1, u1, u2 (sol homogenea)
%           @ calculo del vector fuente
%           @ calculo de la solucion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des conditions aux interface dues aux forces internes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% attention, pour calculer les differentes composantes
% il est nécessaire de subdiviser le calcul en fonction
% des parties paire et impaire (k2)
% cf G11, G22 paires et G12 impaire
% ce qui peut être fait en traitant séparement fint1 et fint2
ncapas      = para.nsubmed;
k2          = DWN.k2;
k1          = DWN.k1;
A_DWN       = DWN.A_DWN;
nk2       	= length(k2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inversion de la matrice %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nk2
    A_DWN(:,:,i)=inv(A_DWN(:,:,i));
end


MAT         = para.reg(1).sub;

Fsource     = zeros(2*(ncapas-1)+1,nk2);
Solution    = zeros(2*(ncapas-1)+1,nk2);

%identification de la couche dans laquelle se trouve la force
ics = 1;
hc  = MAT(1).h;
while zs>hc && ics<ncapas
    ics=ics+1;
    hc=hc+MAT(ics).h;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changement de convention %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1f =k1(:,ics).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul du vecteur source et resolution des coefficients %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% le calcul du vecteur source se fait en utilisant
% la decomposition en ondes planes pour l element j
if ics>1
    htop=0;
    for inch=1:ics-1
        htop=htop+MAT(inch).h;
    end
else
    htop=0;
end
hbot=0;
for inch=1:ics
    hbot=hbot+MAT(inch).h;
end

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
end

%indice de los terminos en el vector fuente
if ics==1
    is120=1;%SIGMA12 EN 0
    if ncapas>1
        is12h=2;%SIGMA12 EN h(ics)
        iu2h =3;%U2 EN h(ics)
    end
elseif ics>1
    is120=(ics-2)*2+1+1;%SIGMA12 EN 0
    iu20=(ics-2)*2+1+2;%U2 EN 0
    if ics<=(ncapas-1)
        is12h=(ics-2)*2+1+3;%SIGMA12 EN h(ics)
        iu2h =(ics-2)*2+1+4;%U2 EN h(ics)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des conditions aux interface dues aux forces internes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ci=MAT(ics).Ci;
Fsource(is120,:)=Ci(6,6)*G22z0;                           %sigma12 en 0
if ics>1
    Fsource(iu20,:)=G220;%U2 EN 0
end
if (ncapas>1 && ics==1) || (ics<ncapas)
    Fsource(is12h,:)=Ci(6,6)*G22zh;                           %sigma12 en h
    Fsource(iu2h,:) =G22h;%U2 EN h
end
mic     = -(-1)^ics;
Fsource = Fsource*mic;

%Resolution du systeme des conditions aux limites
% for i=1:nk2
%     Solution(:,i)   = squeeze(A_DWN(:,:,i))*squeeze(Fsource(:,i));
% end

for i=1:2*(para.nsubmed-1)+1
    Solution(i,:)   = sum(squeeze(A_DWN(i,:,:)).*Fsource,1);
end