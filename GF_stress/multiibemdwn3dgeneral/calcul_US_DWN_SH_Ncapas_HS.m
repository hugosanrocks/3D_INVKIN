function [UXW,SXW,UyKW1]=calcul_US_DWN_SH_Ncapas_HS(para,xr,zr0,izr0,salu,sals,coordf,DWN)
% METHODE : @ recuperacion de k1, u1, u2 (sol homogenea)
%           @ calculo del vector fuente
%           @ calculo de la solucion
%           zr=zr0(izr0)
%           SXW = [SxyKW,SzyKW](xr,zr,w)

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
nk2t        = 2*(nk2-1);%para.DWNnbptkx;
MAT         = para.reg(1).sub;
nrec        = length(xr);
nrecz       = length(zr0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inversion de la matrice %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nk2
    A_DWN(:,:,i)=inv(A_DWN(:,:,i));
end


Kpos2       = k2;
Kneg2       = k2(2:(nk2-1));
DK          = DWN.dk2;
DK          = [DK,DK(end-1:-1:2)];

UyKW        = zeros(1,nk2);
SzyKW       = zeros(1,nk2);
SxyKW       = zeros(1,nk2);

UyKW0       = zeros(1,nk2t);
SzyKW0      = zeros(1,nk2t);
SxyKW0      = zeros(1,nk2t);

xs          = coordf.xs;
nxs       	= length(xs);
zs          = coordf.zs;

UXW         = zeros(nrec,nxs);
SXW         = zeros(2,nrec,nxs);

Fsource     = zeros(2*(ncapas-1)+1,nk2);
Solution    = zeros(2*(ncapas-1)+1,nk2);

%identification de la couche dans laquelle se trouve la force
ics = 1;
hc  = MAT(1).h;
while zs>hc && ics<ncapas
    ics=ics+1;
    hc=hc+MAT(ics).h;
end

k1f =k1(:,ics).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des conditions aux interface dues aux forces internes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%la convolution se fait directement au niveau de la source
%il faut donc separer les ondes qui montent de celles qui descendent et
%celles qui vont à droite de celles qui vont a gauche
mxspecs=1;      %para las ondas que suben
mxspecb=1;      %para las ondas que bajan
nopair=false;   %bandeja indicando si el espectro en k2 es par o no
if nxs==1
    %si hay solo una fuente se hace la convolucion con su extension spacial
    mxspec=1;
    if isfield(coordf,'vnx')
        %caso de un elemento (fuente virtual) de extension proyectada lsegx lsegz
        lsegx   = abs(coordf.dr*coordf.vnz);
        lsegz   = abs(coordf.dr*coordf.vnx);
        mxspecs = sinc(-lsegx*Kpos2/2/pi-lsegz*k1f/2/pi);%(-k1f ,+k2)
        mxspecb = sinc(-lsegx*Kpos2/2/pi+lsegz*k1f/2/pi);%(+k1f ,+k2)
        mxspecsm= sinc(+lsegx*Kpos2/2/pi-lsegz*k1f/2/pi);%(-k1f ,-k2)
        mxspecbm= sinc(+lsegx*Kpos2/2/pi+lsegz*k1f/2/pi);%(+k1f ,-k2)
        nopair  = true;
        
        if (lsegx==0) || (lsegz==0)
            %ds ce cas le sinc est pair
            nopair=false;
        end
    end
end

%%%%%%%%%%%%%%%%%
% interface sup %
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des deplacements %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
zt      = htop-zs;%x=z dans l epaisseur
exp2t   = exp(-1i.*k1f.*abs(zt)).*mxspecs;
fac     = 1/(4*pi*MAT(ics).Ci(6,6));
G220    =-1i*fac./k1f.*exp2t;%pair en k2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des derivees utiles au calcul des contraintes normales %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% derivee p/x1 en X (x1=z)
%pair
G22z0   =-sign(zt)*fac*exp2t;

if nopair
    exp2tm 	= exp(-1i.*k1f.*abs(zt)).*mxspecsm;
    G220m 	=-1i*fac./k1f.*exp2tm;%pair en k2
    G22z0m	=-sign(zt)*fac*exp2tm;
end

%%%%%%%%%%%%%%%%%
% interface inf %
%%%%%%%%%%%%%%%%%
if ics~=ncapas
    zb      = hbot-zs;%x=z dans l epaisseur
    exp2b   = exp(-1i.*k1f.*abs(zb)).*mxspecb;
    
    %pair en k2
    G22h    =-1i*fac./k1f.*exp2b;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcul des derivees utiles au calcul des contraintes normales %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % derivee p/x1 en X (x1=z)
    %pair
    G22zh   =-sign(zb)*fac*exp2b;
    
    if nopair
        exp2bm   = exp(-1i.*k1f.*abs(zb)).*mxspecbm;
        G22hm    =-1i*fac./k1f.*exp2bm;
        G22zhm   =-sign(zb)*fac*exp2bm;
    end
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

if nopair
    Fsourcem         = 0*Fsource;
    Fsourcem(is120,:)= Ci(6,6)*G22z0m;                           %sigma12 en 0
    if ics>1
        Fsourcem(iu20,:)=G220m;%U2 EN 0
    end
    if (ncapas>1 && ics==1) || (ics<ncapas)
        Fsourcem(is12h,:)=Ci(6,6)*G22zhm;                           %sigma12 en h
        Fsourcem(iu2h,:) =G22hm;%U2 EN h
    end
    Fsourcem = Fsourcem*mic;
end

%Resolution du systeme des conditions aux limites
% for i=1:nk2
%     Solution(:,i)   = squeeze(A_DWN(:,:,i))*squeeze(Fsource(:,i));
% end
if ncapas==1
    Solution   = Fsource.*squeeze(A_DWN(1,:,:)).';
    if nopair
        Solutionm   = Fsourcem.*squeeze(A_DWN(1,:,:)).';
    end
else
    for i=1:(2*(ncapas-1)+1)
        Solution(i,:)   = sum(squeeze(A_DWN(i,:,:)).*Fsource,1);
    end
    if nopair
        Solutionm=0*Solution;
        for i=1:(2*(ncapas-1)+1)
            Solutionm(i,:)   = sum(squeeze(A_DWN(i,:,:)).*Fsourcem,1);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des champs au niveau du recepteur %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for irz=1:nrecz
    ixr     = find(izr0==irz);
    nrecx   = length(ixr);
    
    goU=sum(salu(ixr(1:nrecx)));
    goS=sum(sals(ixr(1:nrecx)));
    
    
    %Ecriture de la solution en fonction de la position du recepteur
    zricr   = zr0(irz);
    
    icr     = 1;
    while zricr>MAT(icr).h && icr<ncapas
        zricr   = zricr-MAT(icr).h;
        icr     = icr+1;
    end
    
    exp1=exp(-1i.*k1(:,icr).*zricr).';
    exp3=exp( 1i.*k1(:,icr).*(zricr-MAT(icr).h)).';
    
    if goU>0
        %calcul des deplacements diffractes
        if icr<ncapas
            UyKW = ...
                +Solution(1+(icr-1)*2,:).* exp1 ...
                +Solution(2+(icr-1)*2,:).* exp3;
        else
            UyKW = Solution(1+(icr-1)*2,:).* exp1;
        end
        UyKW(isnan(UyKW))=0;
        
        if nopair
            if icr<ncapas
                UyKWm = ...
                    +Solutionm(1+(icr-1)*2,:).* exp1 ...
                    +Solutionm(2+(icr-1)*2,:).* exp3;
            else
                UyKWm = Solutionm(1+(icr-1)*2,:).* exp1;
            end
            UyKWm(isnan(UyKWm))=0;
        end
    end
    
    if goS>0
        %calcul des contraintes diffractees
        %pair en k2
        C0= -1i*MAT(icr).Ci(6,6)*k1(:,icr).';
        if icr<ncapas
            SzyKW = ...
                +Solution(1+(icr-1)*2,:).*C0.* exp1 ...
                -Solution(2+(icr-1)*2,:).*C0.* exp3;
        else
            SzyKW =Solution(1+(icr-1)*2,:).*C0.* exp1;
        end
        SzyKW(isnan(SzyKW))=0;
        
        %impair en k2
        C0= -1i*MAT(icr).Ci(6,6)*k2;
        if icr<ncapas
            SxyKW = ...
                +Solution(1+(icr-1)*2,:).*C0.* exp1 ...
                +Solution(2+(icr-1)*2,:).*C0.* exp3;
        else
            SxyKW =Solution(1+(icr-1)*2,:).*C0.* exp1;
        end
        SxyKW(isnan(SxyKW))=0;
        
        if nopair
            C0= -1i*MAT(icr).Ci(6,6)*k1(:,icr).';
            if icr<ncapas
                SzyKWm = ...
                    +Solutionm(1+(icr-1)*2,:).*C0.* exp1 ...
                    -Solutionm(2+(icr-1)*2,:).*C0.* exp3;
            else
                SzyKWm =Solutionm(1+(icr-1)*2,:).*C0.* exp1;
            end
            SzyKWm(isnan(SzyKWm))=0;
            
            %impair en k2
            C0= 1i*MAT(icr).Ci(6,6)*k2;
            if icr<ncapas
                SxyKWm = ...
                    +Solutionm(1+(icr-1)*2,:).*C0.* exp1 ...
                    +Solutionm(2+(icr-1)*2,:).*C0.* exp3;
            else
                SxyKWm =Solutionm(1+(icr-1)*2,:).*C0.* exp1;
            end
            SxyKWm(isnan(SxyKWm))=0;
        end
    end
    
    %rajouts des champs incidents % pb difficile a integrer selon x et z
    %     if icr==ics
    %         zrs=zr0(irz)-zs;
    %         exp2t= exp(-1i.*k1f.*abs(zrs));
    %
    %         fac=1/(4*pi*MAT(icr).Ci(6,6));
    %
    %         if goU>0
    %             %pair en k2
    %             G22r=-1i*fac./k1f.*exp2t;
    %
    %             %calcul des deplacements incidents
    %             UyKW = UyKW + G22r;
    %         end
    %
    %         if goS>0
    %             %%%derivee p/x1 en X (x1=z)
    %             %pair
    %             G22z=-sign(zrs)*fac*exp2t;
    %
    %             % derivee p/x2 en X (x2=x)
    %             %impair
    %             G22x=-k2.*fac./k1f.*exp2t;
    %
    %             %calcul des tractions t1j
    %             SzyKW = SzyKW + MAT(icr).Ci(6,6)*G22z;
    %             SxyKW = SxyKW + MAT(icr).Ci(6,6)*G22x;
    %         end
    %     end
    
    %rajout de la phase sur x et inversion
    for ixs=1:nxs
        if nxs~=1
            if ~isfield(coordf,'vnx')
                mxspec=1;
            else
                if abs(coordf.vnz(ixs))>2*abs(coordf.vnx(ixs))
                    lsegx=coordf.dr(ixs)*coordf.vnz(ixs);
                    mxspec=sinc(lsegx*Kpos2/2/pi);
                elseif abs(coordf.vnx(ixs))>2*abs(coordf.vnz(ixs))
                    lsegz=coordf.dr(ixs)*coordf.vnx(ixs);
                    mxspec=sinc(lsegz*k1f/2/pi);
                else
                    mxspec=1;
                end
            end
        end
        
        for irx=1:nrecx
            ir = ixr(irx);
            if salu(ir)==1
                %correction spectre
                UyKW0(1:nk2)             = UyKW(1:nk2).*mxspec;
                if nopair && nxs==1
                    UyKW0(nk2t:-1:(nk2+1))   = UyKWm(2:(nk2t/2))  .*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                else
                    UyKW0(nk2t:-1:(nk2+1))   = UyKW0(2:(nk2t/2))  .*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                end
                UyKW0(1:nk2)             = UyKW0(1:nk2)       .*exp(-1i*(xr(ir)-xs(ixs))*Kpos2);
                %                 tmp=DK*nk2t*ifft(UyKW0,nk2t);
                %                 UXW(ir,ixs)=tmp(1);
                UXW(ir,ixs) = sum(DK.*UyKW0);
                if para.DWNnbptkx/2+1>para.nkmaxKW
                    UyKW1   = 0;
                else
                    UyKW1 	= UyKW0(1:nk2);
                end
                
                %rajout du champ incident, integration gaussienne
                if icr==ics
                    rij     = sqrt((xs(ixs)-xr(ir))^2+(zs-zr0(irz))^2);
                    ksi 	= para.reg(1).sub(ics).ksi;
                    mu      = MAT(ics).Ci(6,6);
                    if isfield(coordf,'vnx')
                        dr  = coordf.dr(ixs);
                        if rij==0 %a faire en dehors des cycles for por accelerer le calcul
                            UXW(ir,ixs) = UXW(ir,ixs) + G22_SH_r_eq_0(ksi,dr,mu);
                        elseif rij<=1*para.npplo*dr && isfield(coordf,'vnx') %a faire en dehors des cycles for por accelerer le calcul
                            coordf0.z   = zs;
                            coordf0.x   = xs(ixs);
                            coordf0.vnx = coordf.vnx(ixs);
                            coordf0.vnz = coordf.vnz(ixs);
                            coordf0.dr  = dr;
                            UXW(ir,ixs) = UXW(ir,ixs) + G22_SH_r_small(coordf0,xr(ir),zr0(irz),1,ksi,para.gaussian,mu);
                        else
                            UXW(ir,ixs) = UXW(ir,ixs) + G22_SH(ksi,rij,mu);
                        end
                    else
                        UXW(ir,ixs) = UXW(ir,ixs) + G22_SH(ksi,rij,mu);
                    end
                end
            end
            
            if sals(ir)==1
                %correction spectre
                SxyKW0(1:nk2)           = SxyKW(1:nk2).*mxspec;
                if nopair && nxs==1
                    SxyKW0(nk2t:-1:(nk2+1)) = SxyKWm(2:(nk2t/2)).*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                else
                    SxyKW0(nk2t:-1:(nk2+1)) =-SxyKW0(2:(nk2t/2)).*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                end
                SxyKW0(1:nk2)           = SxyKW0(1:nk2).*exp(-1i*(xr(ir)-xs(ixs))*Kpos2);
                %                 tmp=DK*nk2t*ifft(SxyKW0,nk2t);
                %                 SXW(1,ir,ixs)=tmp(1);
                SXW(1,ir,ixs)=sum(DK.*SxyKW0);
                
                %correction spectre
                SzyKW0(1:nk2)           = SzyKW(1:nk2).*mxspec;
                if nopair && nxs==1
                    SzyKW0(nk2t:-1:(nk2+1)) = SzyKWm(2:(nk2t/2)).*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                else
                    SzyKW0(nk2t:-1:(nk2+1)) = SzyKW0(2:(nk2t/2)).*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                end
                SzyKW0(1:nk2)           = SzyKW0(1:nk2).*exp(-1i*(xr(ir)-xs(ixs))*Kpos2);
                %                 tmp=DK*nk2t*ifft(SzyKW0,nk2t);
                %                 SXW(2,ir,ixs)=tmp(1);
                SXW(2,ir,ixs)=sum(DK.*SzyKW0);
                
                %rajout du champ incident, integration gaussienne
                if icr==ics
                    xij     = xr(ir)    - xs(ixs);
                    zij     = zr0(irz) 	- zs;
                    rij     = sqrt(xij.^2+zij.^2);
                    ksi 	= para.reg(1).sub(ics).ksi;
                    if isfield(coordf,'vnx')
                        dr      = coordf.dr(ixs);
                        if rij<=1*para.npplo*dr
                            coordf0.z   = zs;
                            coordf0.x   = xs(ixs);
                            coordf0.vnx = coordf.vnx(ixs);
                            coordf0.vnz = coordf.vnz(ixs);
                            coordf0.dr  = dr;
                            [Sxy,Szy]   = Sxzy_SH_r_small(coordf0,xr(ir),zr0(irz),ksi,para.gaussian);
                            SXW(1,ir,ixs) = SXW(1,ir,ixs) + Sxy;
                            SXW(2,ir,ixs) = SXW(2,ir,ixs) + Szy;
                        else
                            tmp = 1i/4*ksi.*besselh(1,2,ksi.*rij);
                            SXW(1,ir,ixs) = SXW(1,ir,ixs) + tmp.*xij./rij;
                            SXW(2,ir,ixs) = SXW(2,ir,ixs) + tmp.*zij./rij;
                        end
                    else
                        tmp = 1i/4*ksi.*besselh(1,2,ksi.*rij);
                        SXW(1,ir,ixs) = SXW(1,ir,ixs) + tmp.*xij./rij;
                        SXW(2,ir,ixs) = SXW(2,ir,ixs) + tmp.*zij./rij;
                    end
                end
            end
        end
    end
end