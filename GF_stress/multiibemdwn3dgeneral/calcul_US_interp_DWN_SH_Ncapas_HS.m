function [UXW,SXW]=calcul_US_interp_DWN_SH_Ncapas_HS(para,xr,zr0,izr0,salu,sals,coordf,zinterp,DWN,Solutionint)
% METHODE : @ recuperacion de k1, u1, u2 (sol homogenea)
%           @ calculo del vector fuente
%           @ calculo de la solucion
%           zr=zr0(izr0)
%           SXW = [SxyKW,SzyKW](xr,zr,w)

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

nk2       	= length(k2);
nk2t        = 2*(nk2-1);%para.DWNnbptkx;
MAT         = para.reg(1).sub;
nrec        = length(xr);
nrecz       = length(zr0);

Kpos2       = k2;
Kneg2       = k2(2:(nk2-1));
% DK          = k2(2)-k2(1);
DK          = DWN.dk2;
DK          = [DK,DK(end-1:-1:2)];
UyKW        = zeros(1,nk2);
SzyKW       = zeros(1,nk2);
SxyKW       = zeros(1,nk2);
UyKW0       = zeros(1,nk2t);
SzyKW0      = zeros(1,nk2t);
SxyKW0      = zeros(1,nk2t);

xs0         = coordf.xs;
nxs       	= length(xs0);
zs0         = coordf.zs;

UXW         = zeros(nrec,nxs);
SXW         = zeros(2,nrec,nxs);

% identification de la couche dans laquelle se trouve la force,
% (tout le segment se situe dans la meme couche)
ics = 1;
hc  = MAT(1).h;
while zs0>hc && ics<ncapas
    ics=ics+1;
    hc=hc+MAT(ics).h;
end

% recuperation de la composante verticale du nombre d onde
k1f     = k1(:,ics).';

% etape de decision si l integration gaussienne est requise ou non
ngau    = 1;
wgau    = 1;
xgau    = 0;
nxs0    = 1;
for ixs=1:nxs
    if abs(coordf.vnz(ixs))<=1*abs(coordf.vnx(ixs))
        ngau    = para.gDWN.ngau;
        wgau    = 0.5*para.gDWN.wgau;
        xgau    = para.gDWN.xgau;
        nxs0    = nxs;
        nxs     = 1;
        break
    end
end

% cycle sur chacun des segments, on rappelle que plusieurs segments sont
% centrés sur le meme zs (ils ont les mêmes coefficients solutions) mais
% qu'ils possedent different xs
% si l'integration gaussienne n est pas requise, le cycle est transparant
% si l'integration gaussienne est requise, il importe peu que les
% coefficients solutions soient identiques pour le centre du segment
% ### on peut ameliorer cette etape en considerant d eventuelle symetrie pour
% des vecteur normaux aux segments.
for ik = 1:nxs0
    % cycle sur les pts d integration, si celle-ci n'est pas requise, le
    % cycle est transparant
    for kg = 1:ngau
        if ngau ==1
            xs  = xs0;
            zs	= zs0;
        else
            xs 	= xs0(ik)-(coordf.vnz(ik).*coordf.dr(ik))*xgau(kg)*0.5;
            zs	= zs0    +(coordf.vnx(ik).*coordf.dr(ik))*xgau(kg)*0.5;
        end
        
        %identification des coefficients associes a la source
        % Solution =Solutionint(:,:,1);%initialisation
        % for j=1:(2*(para.nsubmed-1)+1)
        %     for i=1:nk2
        %         Solution(j,i)= interp1(zinterp,squeeze(Solutionint(j,i,:)),zs,'spline');
        %     end
        % end
        
        % %interpolation polynome degree 2
        % nz      = length(zinterp);
        % tmp     = find(zs<=zinterp+1e-6,1,'first');
        % if isempty(tmp)
        %     indz(3) = nz;
        %     indz(2) = nz-1;
        %     indz(1) = nz-2;
        % else
        %     indz(1) = tmp-1;
        %     indz(2) = tmp;
        %     indz(3) = tmp+1;
        %     if indz(3)>nz
        %         indz = indz-1;
        %     end
        % end
        %
        % Solution=Solutionint(:,:,1);%initialisation
        % for j=1:(2*(para.nsubmed-1)+1)
        %     a   = Solutionint(j,:,indz(1));
        %     z0  = zinterp(indz(1));
        %     dz1 = zinterp(indz(2))-zinterp(indz(1));
        %     dz2 = zinterp(indz(3))-zinterp(indz(1));
        %     c   = dz2*Solutionint(j,:,indz(2))-dz1*Solutionint(j,:,indz(3))-a*(dz2-dz1);
        %     b   = (Solutionint(j,:,indz(2))-a-c*dz1^2)/dz1;
        %     Solution(j,:)=a+b*(zs-z0)+c*(zs-z0)^2;
        % end
        
        %interpolation polynome degree 3
        nz      = length(zinterp);
        tmp     = find(zs<=zinterp+1e-6,1,'first');
        if isempty(tmp)
            indz(4) = nz;
            indz(3) = nz-1;
            indz(2) = nz-2;
            indz(1) = nz-3;
        else
            indz(1) = tmp-1;
            indz(2) = tmp;
            indz(3) = tmp+1;
            indz(4) = tmp+2;
            while indz(4)>nz
                indz = indz-1;
            end
            while indz(1)<1
                indz = indz+1;
            end
        end
        
        Solution=Solutionint(:,:,1);%initialisation
        dz  = zinterp(indz(2))-zinterp(indz(1));
        z0  = zinterp(indz(2));
        for j=1:(2*(para.nsubmed-1)+1)
            a   = Solutionint(j,:,indz(2));
            c   = (Solutionint(j,:,indz(3))+Solutionint(j,:,indz(1))-2*Solutionint(j,:,indz(2)))/(2*dz^2);
            d   = (Solutionint(j,:,indz(4))+2*Solutionint(j,:,indz(1))-3*Solutionint(j,:,indz(2))-6*dz^2*c)/(6*dz^3);
            b   = (Solutionint(j,:,indz(3))-a-c*dz^2-d*dz^3)/dz;
            a   = Solutionint(j,:,indz(1));
            Solution(j,:)=a+b*(zs-z0)+c*(zs-z0)^2+d*(zs-z0)^3;
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
                zricr  = zricr-MAT(icr).h;
                icr = icr+1;
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
            end
            
            %             %rajouts des champs incidents
            %             if icr==ics
            %                 zrs=zr0(irz)-zs;
            %                 exp2t= exp(-1i.*k1f.*abs(zrs));
            %
            %                 fac=1/(4*pi*MAT(icr).Ci(6,6));
            %
            %                 if goU>0
            %                     %pair en k2
            %                     G22r=-1i*fac./k1f.*exp2t;
            %
            %                     %calcul des deplacements incidents
            %                     UyKW = UyKW + G22r;
            %                 end
            %
            %                 if goS>0
            %                     %%%derivee p/x1 en X (x1=z)
            %                     %pair
            %                     G22z=-sign(zrs)*fac*exp2t;
            %
            %                     % derivee p/x2 en X (x2=x)
            %                     %impair
            %                     G22x=-k2.*fac./k1f.*exp2t;
            %
            %                     %calcul des tractions t1j
            %                     SzyKW = SzyKW + MAT(icr).Ci(6,6)*G22z;
            %                     SxyKW = SxyKW + MAT(icr).Ci(6,6)*G22x;
            %                 end
            %             end
            
            for ixs=1:nxs
                ixs0    = max(ixs,ik);
                if nxs==1
                    xs1     = xs;
                else
                    xs1     = xs(ixs0);
                end
                if ~isfield(coordf,'vnx')
                    mxspec=1;
                else
                    if abs(coordf.vnz(ixs0))>1*abs(coordf.vnx(ixs0))
                        lsegx=coordf.dr(ixs0)*coordf.vnz(ixs0);
%                         lsegz=coordf.dr(ixs0)*coordf.vnx(ixs0);
                        %                         mxspec=sinc(lsegx*Kpos2/2/pi+lsegz*k1f/2/pi);%ok slmt qd la source est ds le semi espace%+sinc(lsegx*Kpos2/2/pi-lsegz/2*k1f/2/pi);%
                        mxspec=sinc(lsegx*Kpos2/2/pi);
                        %                         lsegz=coordf.dr(ixs0)*coordf.vnx(ixs0);
                        %                         mxspec=mxspec.*sinc(lsegz*k1f/2/pi);%
                    else
                        mxspec=1;
                    end
                end
                
                %                 figure(205);hold on;plot(imag(UyKW(1:nk2)),'k')
                for irx=1:nrecx
                    ir = ixr(irx);
                    if salu(ir)==1
                        %correction spectre
                        UyKW0(1:nk2)             = UyKW(1:nk2).*mxspec;
                        UyKW0(nk2t:-1:(nk2+1))   = UyKW0(2:(nk2t/2))  .*exp( 1i*(xr(ir)-xs1)*Kneg2);
                        UyKW0(1:nk2)             = UyKW0(1:nk2)       .*exp(-1i*(xr(ir)-xs1)*Kpos2);
                        
                        UXW(ir,ixs0)=UXW(ir,ixs0)+wgau(kg)*sum(DK.*UyKW0);
                        if icr==ics && isfield(coordf,'vnx') && kg==1
                            rij     = sqrt((xs1-xr(ir))^2+(zs-zricr)^2);
                            ksi 	= para.reg(1).sub(ics).ksi;
                            mu      = MAT(ics).Ci(6,6);
                            dr      = coordf.dr(ixs0);
                            
                            if rij==0
                                UXW(ir,ixs0) = UXW(ir,ixs0) + G22_SH_r_eq_0(ksi,dr,mu);
                            elseif rij<=1*para.npplo*dr
                                coordf0.z   = zs;
                                coordf0.x   = xs1;
                                coordf0.vnx = coordf.vnx(ixs0);
                                coordf0.vnz = coordf.vnz(ixs0);
                                coordf0.dr  = dr;
                                UXW(ir,ixs0) = UXW(ir,ixs0) + G22_SH_r_small(coordf0,xr(ir),zricr,1,ksi,para.gaussian,mu);
                            else
                                UXW(ir,ixs0) = UXW(ir,ixs0) + G22_SH(ksi,rij,mu);
                            end
                        elseif icr==ics && kg==1
                            UXW(ir,ixs0) = UXW(ir,ixs0) + G22_SH(ksi,rij,mu);
                        end
                        
                    end
                    
                    if sals(ir)==1
                        %correction spectre
                        SxyKW0(1:nk2)           = SxyKW(1:nk2).*mxspec;
                        SxyKW0(nk2t:-1:(nk2+1)) =-SxyKW0(2:(nk2t/2)).*exp( 1i*(xr(ir)-xs1)*Kneg2);
                        SxyKW0(1:nk2)           = SxyKW0(1:nk2).*exp(-1i*(xr(ir)-xs1)*Kpos2);
                        SXW(1,ir,ixs0)=SXW(1,ir,ixs0)+wgau(kg)*sum(DK.*SxyKW0);
                        
                        %correction spectre
                        SzyKW0(1:nk2)           = SzyKW(1:nk2).*mxspec;
                        SzyKW0(nk2t:-1:(nk2+1)) = SzyKW0(2:(nk2t/2)).*exp( 1i*(xr(ir)-xs1)*Kneg2);
                        SzyKW0(1:nk2)           = SzyKW0(1:nk2).*exp(-1i*(xr(ir)-xs1)*Kpos2);
                        SXW(2,ir,ixs0)=SXW(2,ir,ixs0)+wgau(kg)*sum(DK.*SzyKW0);
                        
                        if icr==ics && isfield(coordf,'vnx') && kg==1
                            xij     = xr(ir)    - xs1;
                            zij     = zricr     - zs;
                            rij     = sqrt(xij.^2+zij.^2);
                            ksi 	= para.reg(1).sub(ics).ksi;
                            dr      = coordf.dr(ixs0);
                            
                            if rij<=1*para.npplo*dr
                                coordf0.z   = zs;
                                coordf0.x   = xs1;
                                coordf0.vnx = coordf.vnx(ixs0);
                                coordf0.vnz = coordf.vnz(ixs0);
                                coordf0.dr  = dr;
                                [Sxy,Szy]=Sxzy_SH_r_small(coordf0,xr(ir),zricr,ksi,para.gaussian);
                                SXW(1,ir,ixs0) = SXW(1,ir,ixs0) + Sxy;
                                SXW(2,ir,ixs0) = SXW(2,ir,ixs0) + Szy;
                            else
                                tmp = 1i/4*ksi.*besselh(1,2,ksi.*rij);
                                SXW(1,ir,ixs0) = SXW(1,ir,ixs0) + tmp.*xij./rij;
                                SXW(2,ir,ixs0) = SXW(2,ir,ixs0) + tmp.*zij./rij;
                            end
                        elseif icr==ics && kg==1
                            tmp = 1i/4*ksi.*besselh(1,2,ksi.*rij);
                            SXW(1,ir,ixs0) = SXW(1,ir,ixs0) + tmp.*xij./rij;
                            SXW(2,ir,ixs0) = SXW(2,ir,ixs0) + tmp.*zij./rij;
                        end
                    end
                end
            end
        end
    end
end