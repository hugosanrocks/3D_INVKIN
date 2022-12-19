function [uw,sw]=campo_inc(para)
%calculo del espectro debido a una OP o FP en un espacio infinito
coordr.x	= para.rec.xr.';
coordr.z	= para.rec.zr.';
nrec        = para.rec.nrec;
coordr.nbpt = nrec;

xs          = para.xs;
zs          = para.zs;
ninc        = para.ninc;

if para.fuente==1
    kzsigno = para.kzsigno;
    kxk     = para.kxk;
    polOP   = para.tipo_onda;
else
    fij     = para.fij;
end

Ci          = para.reg(1).Ci;
ksi         = para.reg(1).ksi;
uw          = 0;
sw          = 0;

sal         = para.sortie;

if para.dim==1
    if para.pol==1 %SH
        uw      = zeros(nrec,ninc,1);%sal.Ut+sal.USh+sal.UIh; no tiene sentido en espacio sentido
        nss     = (sal.sxy + sal.syz);
        sw      = zeros(nrec,ninc,nss);%3a dim: sxy, syz
        
        for iinc=1:ninc
            salu=logical(true(nrec,1).*para.sortie.Uy);
            sals=logical(true(nrec,1).*max(para.sortie.sxy,para.sortie.syz));
            if para.fuente==1 %OP
                [uw(:,iinc,1),~,sw(:,iinc,:)]     = campo_ref_SH_OPHeI_SE(xs(iinc),zs(iinc),coordr,ksi,kxk(iinc),kzsigno(iinc),0,Ci(6,6),para);
            elseif para.fuente==2 %FP
                [uw(:,iinc,1),~,sw(:,iinc,:)]     = campo_ref_SH_FP(xs(iinc),zs(iinc),coordr,ksi,Ci(6,6),salu,sals,para);
            end
        end
    elseif para.pol==2%PSV
        kpi     = para.reg(1).kpi;
        kri     = para.reg(1).kri;

        nsu     = (sal.Ux + sal.Uz)*(sal.Ut+sal.UPh+sal.USh+sal.UIh+sal.UPt+sal.USt);
        nss     = (sal.sxx + sal.szz + sal.sxz);

        uw      = zeros(nrec,ninc,nsu);
        sw      = zeros(nrec,ninc,nss);
        
        for iinc=1:ninc
            salu=logical(true(nrec,1).*max(nsu));
            sals=logical(true(nrec,1).*max(nss));
            if para.fuente==1 %OP
                [uw(:,iinc,:),~,~]   = campo_ref_PSV_OPHeI(xs(iinc),zs(iinc),coordr,kpi,ksi,kri,Ci,polOP(iinc),kxk(iinc));
            elseif para.fuente==2 %FP
                [uw(:,iinc,:),~,sw(:,iinc,:)]  = campo_ref_PSV_FP(xs(iinc),zs(iinc),coordr,kpi,ksi,Ci,fij(iinc,:),salu,sals,para);
            end
        end
%         if para.fuente==1
%             uxz0c=correcion_uxz0_OP(uxz0,ns,polOP,para);
%         end
    end
else%3D
    
    kpi = para.reg(1).kpi;
    coordr.y	= para.rec.yr.';
    ys          = para.ys;
    uw          = zeros(nrec,para.ninc,3);
    sw          = zeros(nrec,para.ninc,6);
    for iinc=1:ninc
        if para.fuente==1 %OP
            kri = para.reg(1).kri;
            %             [ uw(j+1,:,:,iinc),~] = campo_ref_OPHeI_3D(xs(iinc),ys(iinc),zs(iinc),coordr,kpi,ksi,kri,Ci,polOP(iinc),kxk(iinc),kyk(iinc));
        elseif para.fuente==2 %FP
            nsu     = sal.Ux + sal.Uy + sal.Uz;
            nss     = sal.sxx + sal.syy + sal.szz + sal.sxy + sal.sxz + sal.syz;

            salu=logical(true(nrec,1).*max(nsu));
            sals=logical(true(nrec,1).*max(nss));
            [ utmp,~,S]  = campo_ref_FP3D(xs(iinc),ys(iinc),zs(iinc),coordr,kpi,ksi,Ci,para.fij(iinc,:),salu,sals,para);
            uw(:,iinc,:)= permute(utmp,[2 1]);
            sw(:,iinc,:)= permute(S,[2 1]);
        end
    end
end