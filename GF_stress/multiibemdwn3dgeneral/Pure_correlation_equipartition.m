function [U,S] = Pure_correlation_equipartition(para,fj)

%para el calculo de las funciones de Green a partir de
%la equiparticion y de las correlaciones

DWN.omegac	= 2*pi*fj;
% para.DWNomei= 1e-1;

%posiciones original de los receptores
xr          = para.rec.xr;
zr          = para.rec.zr;
%-------------------------------------------------%
%reorganisacion de las fuentes y de los receptores%
%(por ahora 1 fuente y 1 receptor)                %
%-------------------------------------------------%
%la fuente se vuelve temporalmente receptor #1
xr              = [para.xs(1) xr];
zr              = [para.zs(1) zr];
nr              = length(xr);
%para acelerar el calculo, se busca las profundidades eguales
[zr0,izr0,~]    = pseudo_unique(zr,para);
salu            = ones(nr,1);
sals            = ones(nr,1);
%despues de la corelacion regresamos al receptor inicial
nr              = 1;
%el origen de las OP de volumen es cualquiera, pero para las ondas de
%superficie, la normalizacion se hace tomando en cuienta una referencia
%a la altura de la ultima interface con el semi espacio
para.xs         = mean(xr);
zs=0;
for jh=1:para.nsubmed-1
    zs=zs+para.reg(1).sub(jh).h;
end
para.zs         = zs;

%velocidades
beta            = para.reg(1).sub(para.nsubmed).bet;
if para.pol==2
    alpha     	= para.reg(1).sub(para.nsubmed).alpha;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculo contribucion de las ondas de volumen %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kz reel
% cmpt=0;
% for nincBW=para.nincBWOP:-2:2

nr0         = length(izr0);
sals        = zeros(nr0,1);


%no se puede agregar attenuacion por ahora en las ondas incidentes
% %ni tampoco la attenuacion numerica en w
% for ms=1:para.nsubmed
%     para.reg(1).sub(ms).tipoatts=0;
% end
% para        = attenuation(para,real(fj));


%     if nincBW==para.nincBWOP
%         uref=uBW;
%     elseif abs(uref-uBW)/abs(uref)>1e-2%
%         cmpt=cmpt+1;
%         if cmpt==3
%             break
%         end
%     else
%         cmpt=0;
%     end
% end
%     figure(108);hold on;plot(real(fj),ninc2,'k.')
%     uBW=uref;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculo contribucion de las ondas de superficie %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if para.nsubmed>1
    %calculo exacto de las posiciones kx y vg a partir de las curvas de dispersion
    w               = DWN.omegac;
    [k20ok,vgS]     = dispersion_curve_vg_fw(para,real(w));
    %numero de modos de superficie
    para.nincS      = length(vgS);
    para.kxk        = real(k20ok./para.reg(1).sub(para.nsubmed).ksi);
    % para.kzsigno    = ones(para.nincS,1);%*para.kzsigno(1);
    para.xs         = ones(para.nincS,1)*para.xs(1);
    para.zs         = ones(para.nincS,1)*para.zs(1);
end
%identificacion de la profundidad de la ultima interfase
%se normaliza a partir de esta profundidad

zs=0;
for jh=1:para.nsubmed-1
    zs=zs+para.reg(1).sub(jh).h;
end

if para.dim ==1
    sal  	= para.sortie;
    
    if  para.pol==1 %SH
        para     = attenuation(para,real(fj)-1i*1e-6);

        nss         = sal.sxy+sal.syz;
        S           = zeros(nr,para.ninc,nss);
        
        % calculo contribucion de las ondas de volumen %
        nincBW   	= 500;%para.nincBWOP;
        ninc2       = nincBW;
        angletot    = pi/2;
        dth         = angletot/ninc2;
        gam         = linspace(dth/2,angletot-dth/2,ninc2);
        kxk         = cos(gam);
        
        
        coordf.xs   = para.xs(1)*ones(1,ninc2);
        coordf.zs   = para.zs(1)*ones(1,ninc2);
        
        [uS,~]      = calcul_US_DWN_SH_Ncapas_HS_incOP(para,xr,zr0,izr0,salu,sals,coordf,kxk,1);
        uBW         = real(sum(uS(1,:).*conj(uS(2,:))))/ninc2;
        
        % calculo contribucion de las ondas de superficie %
        % init
        coordf.xs   = 0;
        coordf.zs   = zs;
        uL          = 0;
        if para.nsubmed>1
            % suma sobre los diferentes modos (o estados) de superficie
            for iincS=1:para.nincS %a vectoriser
                %kz imag
                pesoS    	= 2*beta^2/w/vgS(iincS);
                kxk         = para.kxk(iincS);
                [u,~]       = calcul_US_DWN_SH_Ncapas_HS_incOPI(para,xr,zr0,izr0,salu,sals,coordf,kxk,1);
                uL          = uL + 2*real(pesoS*u(1)*conj(u(2)));%factor 2 por simetria de la parte real y antisymetria parte imaginaria con "+ u(2)*conj(u(1))"
            end
        end
        % suma superficie y volumen y normalizacion E0
        U              = -1i*1/4/para.reg(1).sub(para.nsubmed).bet^2/2*(uL+uBW);
    elseif para.pol==2%P-SV
        nss   	= sal.sxx + sal.szz+sal.sxz;
        S       = zeros(nss,nr,para.ninc);
        
        % calculo contribucion de las ondas de volumen %
        nincBW   	= 10;%para.nincBWOP;

        angletot    = pi;
        dth         = angletot/nincBW;
        gam         = linspace(-angletot/2+dth/2,angletot/2-dth/2,nincBW);
        kxk         = sin(gam);
        coordf.xs   = para.xs(1)*ones(1,nincBW);
        coordf.zs   = zs*ones(1,nincBW);

        para        = attenuation(para,real(fj)-1i*1e-3);
        R           = para.reg(1).sub(para.nsubmed).alpha/para.reg(1).sub(para.nsubmed).bet;
        pesoP       = 1/(1+R^2);
        [uP,~]      = calcul_US_DWN_PSV_Ncapas_HS_incOP(para,xr,zr0,izr0,salu,sals,coordf,kxk,1);%P
        
        pesoS       = R^2/(1+R^2);
        [uS,~]      = calcul_US_DWN_PSV_Ncapas_HS_incOP(para,xr,zr0,izr0,salu,sals,coordf,kxk,2);%S
        
        uBW         = zeros(2,2);
        for i=1:2
            for j=1:2
                uBW(j,i)     = ...
                    pesoP*real(sum(uP(i,1,:).*conj(uP(j,2,:))))/nincBW + ...
                    pesoS*real(sum(uS(i,1,:).*conj(uS(j,2,:))))/nincBW;
            end
        end

        
        % calculo contribucion de las ondas de superficie %
        uS          = zeros(2,2);
        if para.nsubmed==1
            vP          = alpha;
            vS          = beta;
            R           = vP/vS;
            
            w           = 1;%en realidad indepediente de w
            x0          = 1/R^2;
            tmpR        = roots([-1,8,8*(2*x0-3),16*(-x0+1)]);
            vR          = vS*sqrt(tmpR(abs(tmpR)<1));
            
            chi1        = sqrt(1-(vR/vP)^2);
            chi2        = sqrt(1-(vR/vS)^2);
            Int_E_Ray_dz= vR/w*(chi2-chi1)*(1/(2*chi1*chi2)-2*sqrt(chi1/chi2)/(chi1+chi2));
            pesoR       = vP^2/w/vR/(1+R^2)/Int_E_Ray_dz/2;

            para.ninc   = 1; % 2 ondas de cada lado, pero se ocupa symetrias
            coordf.xs   = 0;
            coordf.zs   = 0;
            [uR,~]      = calcul_US_DWN_PSV_Ncapas_HS_incOP(para,xr,zr0,izr0,salu,sals,coordf,1,3);%S
            
            for i=1:2
                for j=1:2
                    uS(j,i)     = pesoR*real(uR(i,1,:).*conj(uR(j,2,:))+conj(uR(i,1,:)).*uR(j,2,:));
                end
            end
        else
            para        = attenuation(para,real(fj)-1i*1e-4);
            for iincS=1:para.nincS % suma sobre los diferentes modos (o estados) de superficie
                pesoR    	= 2*alpha^2/w/vgS(iincS)/(1+R^2);
                kxk         = para.kxk(iincS);
                [uR,~]      = calcul_US_DWN_PSV_Ncapas_HS_incOPI(para,xr,zr0,izr0,salu,sals,kxk);
                for i=1:2
                    for j=1:2
                        uS(j,i)     =  uS(j,i) + pesoR*real(uR(i,1).*conj(uR(j,2)));
                    end
                end
                kxk         =-para.kxk(iincS);
                [uR,~]      = calcul_US_DWN_PSV_Ncapas_HS_incOPI(para,xr,zr0,izr0,salu,sals,kxk);
                for i=1:2
                    for j=1:2
                        uS(j,i)     =  uS(j,i) + pesoR*real(uR(i,1).*conj(uR(j,2)));
                    end
                end
            end
        end
        U         	= zeros(2,nr,2);
        
        %suma superficie y volumen y normalizacion E0
        U(:,1,:) 	= -1i/pesoS/4/beta^2/2*(uS+uBW);
        U         	= permute(U,[2 3 1]);%(2,nrec,ninc)->(nrec,ninc,2)
    end
elseif para.dim>=3
    
    rec     = para.rec;
    kr   	= (0.5+(0:nk))*DK/2;
    DWN.kr	= kr;
    DWN.k2 	= kr;%copia para rebuildk2, sino inutil
    DWN.dkr	= kr*0+DK;
    
    
    
    sal  	= para.sortie;
    nss   	= sal.sxx + sal.syy + sal.szz+sal.sxy+sal.sxz+sal.syz;
    U       = zeros(3,nr,para.ninc);
    S       = zeros(nss,nr,para.ninc);
    
    if para.nsubmed==1 && para.reg(1).sub(1).h<0
        % campo incidente solamente
        coordf.xs       = para.xs; coordr.x = para.rec.xr;
        coordf.ys       = para.ys; coordr.y = para.rec.yr;
        coordf.zs       = para.zs; coordr.z = para.rec.zr;
        ic =ones(nr,1);
        
        
        [Gc,~,~,~]  = Gij_3D_DWN_polar(para,DWN,coordr,coordf,ic,1);
        for iinc=1:para.ninc
            fij       	= para.fij(iinc,:);
            U(1,:,iinc) = fij(1)*Gc(1,1,:,iinc)+fij(2)*Gc(1,2,:,iinc)+fij(3)*Gc(1,3,:,iinc);
            U(2,:,iinc) = fij(1)*Gc(2,1,:,iinc)+fij(2)*Gc(2,2,:,iinc)+fij(3)*Gc(2,3,:,iinc);
            U(3,:,iinc) = fij(1)*Gc(3,1,:,iinc)+fij(2)*Gc(3,2,:,iinc)+fij(3)*Gc(3,3,:,iinc);
        end
        
        
        %         Sp=[Srrx Srtx Szrx; Srtx Sttx Sztx; Szrx Sztx Szzx];
        %         T =[ ct    -st    0  ;   st    ct    0  ;   0     0     1  ];
        %         Sc= T*Sp*(T.');
        %         %Sc=[SxxxW SxyxW SxzxW; SyxxW SyyxW SyzxW; SzxxW SzyxW SzzxW];
        %
        %         SXW_fx(:,:,ir,ixs) = Sc;
        %
        %         %cf S_3D
        %         SXW_fy(1,1,ir,ixs) = SXW_fx(2,2,ir,ixs);
        %         SXW_fy(1,2,ir,ixs) =-SXW_fx(1,2,ir,ixs);
        %         SXW_fy(1,3,ir,ixs) =-SXW_fx(2,3,ir,ixs);
        %         SXW_fy(2,1,ir,ixs) = SXW_fy(1,2,ir,ixs);
        %         SXW_fy(2,2,ir,ixs) = SXW_fx(1,1,ir,ixs);
        %         SXW_fy(2,3,ir,ixs) = SXW_fx(1,3,ir,ixs);
        %         SXW_fy(3,1,ir,ixs) = SXW_fy(1,3,ir,ixs);
        %         SXW_fy(3,2,ir,ixs) = SXW_fy(2,3,ir,ixs);
        %         SXW_fy(3,3,ir,ixs) = SXW_fx(3,3,ir,ixs);
        %
        %         SXW_fx(:,:,ir,ixs) = fij(ixs,1)*SXW_fx(:,:,ir,ixs);
        %         SXW_fy(:,:,ir,ixs) = fij(ixs,2)*SXW_fy(:,:,ir,ixs);
        %
        
        %pour le tenseur des stress, il faut developer Gij_3D_DWN_polar
        %pour donner Sc avec toutes les composantes
        for iinc=1:para.ninc
            xij     = coordr.x - para.xs(iinc);
            yij     = coordr.y - para.ys(iinc);
            zij     = coordr.z - para.zs(iinc);
            rij     = sqrt(xij.^2+yij.^2+zij.^2);
            ksi 	= para.reg(1).sub(1).ksi;
            kpi 	= para.reg(1).sub(1).kpi;
            gam(1,:)  = xij./rij;
            gam(2,:)  = yij./rij;
            gam(3,:)  = zij./rij;
            fij             = para.fij(iinc,:);
            [S_fx,S_fy,S_fz]= S_3D(rij,gam,ksi,kpi,fij);
            S1              = zeros(3,3,nr,para.ninc);
            S1(:,:,:,iinc)   = S_fx*fij(1)+S_fy*fij(2)+S_fz*fij(3);
        end
        j   = 1;
        if sal.sxx; S(j,:,:) = S1(1,1,:,:);j=j+1; end;
        if sal.syy; S(j,:,:) = S1(2,2,:,:);j=j+1; end;
        if sal.szz; S(j,:,:) = S1(3,3,:,:);j=j+1; end;
        if sal.sxy; S(j,:,:) = S1(1,2,:,:);j=j+1; end;
        if sal.sxz; S(j,:,:) = S1(1,3,:,:);j=j+1; end;
        if sal.syz; S(j,:,:) = S1(2,3,:,:);       end;
    else
        DWN         = calcul_A_DWN_3D_polar_Ncapas_HS(para,DWN);
        DWN         = rebuildk2(DWN);
        DWN.kr      = DWN.k2;
        DWN.dkr     = DWN.dk2;
        DWN.dkr(1)  = DWN.kr(2)-DWN.dkr(2)/2;
        DWN         = calcul_A_DWN_3D_polar_Ncapas_HS(para,DWN);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % inversion de la matrice %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:length(DWN.kr)
            DWN.A_DWN(:,:,i)=inv(DWN.A_DWN(:,:,i));
            DWN.B_DWN(:,:,i)=inv(DWN.B_DWN(:,:,i));
        end
        
        for iinc=1:para.ninc
            coordf.xs       = para.xs(iinc);
            coordf.ys       = para.ys(iinc);
            coordf.zs       = para.zs(iinc);
            if para.fuente==1 %OP
                polOP       = para.tipo_onda(iinc);
                kxk         = para.kxk(iinc);
                kzsigno     = para.kzsigno(iinc);
                phi         = para.phi(iinc);
                [u,~]     	= calcul_US_DWN_3D_Ncapas_HS_incOP(para,xr,yr,zr0,izr0,salu,sals,coordf,kxk,polOP,kzsigno,phi);
            else %if para.fuente==2 %FP
                fij       	= para.fij(iinc,:);
                %                 [u_f1,s_f1,u_f2,s_f2,u_f3,s_f3] = calcul_US_DWN_3D_polar_Ncapas_HS(para,rec,salu,sals,coordf,fij,DWN);
                
                nzr0                = length(zr0);
                zricrall            = zeros(1,nzr0);
                icrall              = zeros(1,nzr0);
                for ir=1:nzr0
                    icrv     = 1;
                    zricrv   = zr0(ir); % profundidad relativa a la interface de la capa del receptor
                    while zricrv>para.reg(1).sub(icrv).h && icrv<para.nsubmed
                        zricrv   = zricrv-para.reg(1).sub(icrv).h;
                        icrv     = icrv+1;
                    end
                    zricrall(ir) = zricrv;
                    icrall(ir)   = icrv;
                end
                
                rec.zr0         = zr0;
                rec.izr0        = izr0;
                rec.zricrall    = zricrall;
                rec.icrall      = icrall;
                [u_f1,s_f1,u_f2,s_f2,u_f3,s_f3] = calcul_US_DWN_3D_polar_Ncapas_HS2(para,rec,salu,sals,coordf,fij,DWN);
                
                u           = u_f1+u_f2+u_f3;
                S1          = s_f1+s_f2+s_f3;
            end
            U(:,:,iinc) = u;
            
            j   = 1;
            if sal.sxx; S(j,:,iinc) = S1(1,1,:);j=j+1; end;
            if sal.syy; S(j,:,iinc) = S1(2,2,:);j=j+1; end;
            if sal.szz; S(j,:,iinc) = S1(3,3,:);j=j+1; end;
            if sal.sxy; S(j,:,iinc) = S1(1,2,:);j=j+1; end;
            if sal.sxz; S(j,:,iinc) = S1(1,3,:);j=j+1; end;
            if sal.syz; S(j,:,iinc) = S1(2,3,:);       end;
        end
    end
    U = permute(U,[2 3 1]);%(3,nrec,ninc)->(nrec,ninc,3)
    S = permute(S,[2 3 1]);%(6,nrec,ninc)->(nrec,ninc,6)
end