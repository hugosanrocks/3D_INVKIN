function [U,UKW0,S] = Pure_DWN(para,fj)

DWN.omegac	= 2*pi*fj;
nk          = para.DWNnbptkx;
nk0         = nk;
U = 0; UKW0    = 0; S = 0;
%differente strategie de definir la borne d integration kx
%1)selon une valeur de l'interface, un kx max pour toutes les fq, OK pour
%les graph kw
DK          = 2*para.DWNkmax/nk0;
  
% 
% if para.DWNnbptkx/2+1>para.nkmaxKW
%     %2) un kx max = w/beta min, ok pour SH et slmt partie imaginaire
%     DK          = 2*para.reg(1).sub(para.nsubmed).ksi/nk0;
%     %3) un kx max = w/beta min*fac, ok pour PSV et slmt partie imaginaire
%     DK          = 2*1.3*real(para.reg(1).sub(1).ksi)/nk0;
% end

%xmax doit etre > vmax*tmax=alpha/dfe
k2          = (0:(fix(nk0/2)))*DK;
k2(1)       = k2(2)/1000;
DWN.k2      = k2;
DWN.dk2     = k2*0+DK;

%posiciones de los receptores
xr          = para.rec.xr;
yr          = para.rec.yr; %?
zr          = para.rec.zr;
nr          = length(xr);
[zr0,izr0,~]= pseudo_unique(zr,para);
salu        = ones(nr,1);
sals        = ones(nr,1);


%caso de las ondas de superficie
%calculo exacto de las posiciones kx y vg a partir de las curvas de dispersion
%principalmente para el calculo de las funciones de Green a partir de
%la equiparticion y de las correlaciones
%o para la incidencia de ondas de superficie
if       para.dim==1 && ((para.pol==1 && max(para.tipo_onda==2)==1) || ...
        (para.pol==2 &&  max(para.tipo_onda==3)==1 && para.nsubmed>1)) || ...
        (para.dim>=3 && (max(para.tipo_onda>3)==1))
    w               = real(DWN.omegac);
    [k20ok,vgLove]	= dispersion_curve_vg_fw(para,w);
    %numero de modos de superficie
    para.ninc       = length(vgLove);
    para.kxk        = real(k20ok./para.reg(1).sub(para.nsubmed).ksi);
    para.kzsigno    = ones(para.nincL,1);%*para.kzsigno(1);
    para.xs         = ones(para.nincL,1)*para.xs(1);
    para.zs         = ones(para.nincL,1)*para.zs(1);
end



if para.dim ==1
    sal  	= para.sortie;
    
    if  para.pol==1
        nss = sal.sxy+sal.syz;
        S   = zeros(nss,nr,para.ninc);
        U   = zeros(nr,para.ninc);
        
        for iinc=1:para.ninc
            coordf.xs       = para.xs(iinc);
            coordf.zs       = para.zs(iinc);
            if para.fuente==1 %OP
                kzsigno     = para.kzsigno(iinc);
                kxk         = para.kxk(iinc);
                if para.tipo_onda==1
                    [u,~]     	= calcul_US_DWN_SH_Ncapas_HS_incOP(para,xr,zr0,izr0,salu,sals,coordf,kxk,kzsigno);
                else
                    kzsigno     = para.kzsigno(iinc);
                    kxk         = para.kxk(iinc);
                    [u,~]       = calcul_US_DWN_SH_Ncapas_HS_incOPI(para,xr,zr0,izr0,salu,sals,coordf,kxk,kzsigno);%*1.000001
                    vS          = para.reg(1).sub(para.nsubmed).bet;
                end
            else %if para.fuente==2 %FP
                DWN     = calcul_A_DWN_SH_Ncapas_HS(para,DWN);
                if para.DWNnbptkx/2+1>para.nkmaxKW
                    DWN     = rebuildk2(DWN);
                    DWN     = calcul_A_DWN_SH_Ncapas_HS(para,DWN);
                end
                [u,s,UKW0] = calcul_US_DWN_SH_Ncapas_HS(para,xr,zr0,izr0,salu,sals,coordf,DWN);
                if para.DWNnbptkx/2+1>para.nkmaxKW
                    UKW0=0;
                end
            end
            
            U(:,iinc)       = u;
            j=1;
            if sal.sxy==1
                S(j,:,iinc)= s(1,:);
                j=j+1;
            end
            if sal.syz==1
                S(j,:,iinc)= s(2,:);
            end
        end
        if ~isempty(S)
        S = permute(S,[2 3 1]);%(2,nrec,ninc)->(nrec,ninc,2)
        end
    elseif para.pol==2
        nss   	= sal.sxx + sal.szz+sal.sxz;
        S       = zeros(nss,nr,para.ninc);
        
        DWN     = calcul_A_DWN_PSV_Ncapas_HS(para,DWN);
        %         DWN     = rebuildk2(DWN);
        %         DWN     = calcul_A_DWN_PSV_Ncapas_HS(para,DWN);
        
        U       = zeros(2,nr,para.ninc);
        nk2     = length(DWN.k2);

        for iinc=1:para.ninc
            coordf.xs       = para.xs(iinc);
            coordf.zs       = para.zs(iinc);
            if para.fuente==1 %OP
                polOP       = para.tipo_onda(iinc);
                kxk         = para.kxk(iinc);
                [u,~]     	= calcul_US_DWN_PSV_Ncapas_HS_incOP(para,xr,zr0,izr0,salu,sals,coordf,kxk,polOP);
            else %if para.fuente==2 %FP
                fij       	= fliplr(para.fij(iinc,:));
                [u2,~,u1,~,UKW00] = calcul_US_DWN_PSV_Ncapas_HS(para,xr,zr0,izr0,salu,sals,coordf,fij,DWN);%1& 2 oppose ds IBEM & DWN
                u           = u1+u2;
            end
            U(:,:,iinc) = u;
            %para la graphica KW no mas
            if iinc==1
                UKW0 = UKW00;
            else
                UKW0 = UKW0 + UKW00;
            end
        end
        U = permute(U,[2 3 1]);%(2,nrec,ninc)->(nrec,ninc,2)
    end
elseif para.dim>=3
    
    rec     = para.rec;
    
    xl = para.DWNxl;
    DK = 2*pi/xl;
    nk = fix(para.DWNkmax/DK)+1;
    nk = max(500,nk);
    
    if ~para.usingparfor
      disp([char(8) 'nk= ' num2str(nk)])
    end
    
%     kr   	= (0.5+(0:nk))*DK/2;
    kr   	= (0.5+(0:nk))*DK;
    DWN.kr	= kr;
    DWN.k2 	= kr;%copia para rebuildk2, sino inutil
    DWN.dkr	= kr*0+DK;
    
    UKW0    = 0;

    
    sal  	= para.sortie;
    nss   	= sal.sxx + sal.syy + sal.szz+sal.sxy+sal.sxz+sal.syz;
    U       = zeros(3,nr,para.ninc);
    S       = zeros(nss,nr,para.ninc);
    
    if para.nsubmed==1 && para.reg(1).sub(1).h<0 %espacio completo DWN
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
        
        %aide a utiliser moins de kr
%         DWN         = rebuildk2(DWN);
%         DWN.kr      = DWN.k2;
%         DWN.dkr     = DWN.dk2;
%         DWN.dkr(1)  = DWN.kr(2)-DWN.dkr(2)/2;
%         DWN         = calcul_A_DWN_3D_polar_Ncapas_HS(para,DWN);
        
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
                kxk         = para.kxk(iinc); %cartesiana x de la normal de la OP
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