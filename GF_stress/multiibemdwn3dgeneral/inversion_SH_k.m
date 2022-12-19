function [uy,sw]=inversion_SH_k(para,coord,phi_fv,gaussian,DWN)

%%%%%%%%%%%%%
% unpacking %
%%%%%%%%%%%%%

xr      = para.rec.xr;
zr      = para.rec.zr;
nrec    = para.rec.nrec;

xs      = para.xs;
zs      = para.zs;

if para.fuente==1
    kzsigno = para.kzsigno;
    kxk     = para.kxk;
end

ninc    = para.ninc;
ninc0   = para.ninc;

sal     = para.sortie;
ns      = sal.Ut+sal.USh+sal.UIh;

uy      = zeros(nrec,ninc,ns);
nss     = (sal.sxy + sal.syz);
sw      = zeros(nrec,ninc,nss);%3a dim: sxy, syz
if para.geo(1)==3
    if para.pol==1
        if para.fuente==1
            if abs(para.kxk(1))>1
                % incidencia de ondas heterogeneas
                % el numero de incidencias depende de la frecuencia
                w           = DWN.omegac;
                k20ok       = dispersion_curve_k2_fw(para,DWN,real(w));
                para.ninc   = length(k20ok);
                para.kxk    = sign(para.kxk(1))*real(k20ok./para.reg(1).sub(para.nsubmed).ksi);
                para.kzsigno= ones(para.ninc,1)*para.kzsigno(1);
                para.xs     = ones(para.ninc,1)*para.xs(1);
                para.zs     = ones(para.ninc,1)*para.zs(1);
                para.xzs    = ones(para.ninc,1)*para.xzs(1);
                ninc        = para.ninc;
            end
        end
    end
end

% calculo de las soluciones a partir de los coeficientes
for m=1:para.nmed0
    %el receptor esta en el medio m
    irec        = para.rec.m(m).ind;
    nrecm       = length(irec);
    mmat        = para.subm1(m);
    
    if nrecm~=0
        if m==1 && para.geo(1)==3
            mu	= para.reg(1).sub(1).Ci(6,6);
            ksi = para.reg(1).sub(1).ksi;
        else
            mu	= para.reg(mmat).Ci(6,6);
            ksi = para.reg(mmat).ksi;
        end
        
        coordr.x    = xr(irec).';
        coordr.z    = zr(irec).';
        coordr.nbpt = nrecm;
        
        %campo incidente
        uy0	= zeros(nrecm,ninc0);
        s0 	= zeros(nss,nrecm,ninc0);
        for iinc=1:ninc
            salu=logical(true(nrecm,1)*(sal.Uy));
            sals=logical(true(nrecm,1).*(nss>0 ));
            if para.xzs(iinc)==m %fuente en m
                if m==1 && para.nsubmed>1 && para.geo(1)==3
                    %DWN con varias capas,
                    %calculo ya hecho anteriormente en vector_fuente_XXX
                    %se recupera el valor
                    uy0(:,iinc)         = DWN.uy0(:,iinc);
                    s0(:,:,iinc)        = DWN.s0(:,:,iinc);
                elseif para.fuente==1 %OP
                    [uy0(:,iinc),~,s0(:,:,iinc)]     = campo_ref_SH_OPHeI_SE(xs(iinc),zs(iinc),coordr,ksi,kxk(iinc),kzsigno(iinc),0,mu,para);
                elseif para.fuente==2 %FP
                    [uy0(:,iinc),~,s0(:,:,iinc)]     = campo_ref_SH_FP(xs(iinc), zs(iinc),coordr,ksi,mu,salu,sals,para);
                    %casos en los cuales se debe de considerar la contribucion de la fuente real imagen
                    if (para.geo(1)==2 && para.cont(1).ruggeo==1  && m==1 && para.fuenteimagen==1 || para.geo(1)==3 && m==1)
                        [uy1,~,s1]       	= campo_ref_SH_FP(xs(iinc),-zs(iinc),coordr,ksi,mu,salu,sals,para);
                        uy0(:,iinc)         = uy0(:,iinc)+uy1.';
                        s0(:,:,iinc)        = squeeze(s0(:,:,iinc))+s1;
                    end
                end
            else
                uy0(:,iinc) = 0;
            end
        end
        
        %adjuntar campo difractado a campo incidente
        if m == 1 && para.geo(1)==3 && para.nsubmed > 1 
            %IBEM-DWN
            dr      = coord.dr;
            phi     = coord.phi;
            
            %posicion de las fv que hay que tomar en cuenta
            j       = 1:coord.nbpt;
            jphi    = 1:coord.nbeq;
            
            %indice (logic y natural) de los puntos de colocacion perteneciendo a m
            jjx     = coord.indm(1).ind;
            ii      = j(jjx);
            
            %indice (logic y natural) de los phi que hay que contemplar (columnas
            %de la matriz)
            jjphi   = false(coord.nbeq,1);
            jjphi(phi(j(jjx),m)) = true(1);
            jj      = jphi(jjphi);
            
            drj = dr(ii);
            %integracion de las contribuciones
            for iinc = 1:ninc
                %campo difractado
                for i=1:nrecm
                    uydiff          = sum(DWN.uydiff(jj,i).'.*drj.*phi_fv(jj,iinc).');
                    uy(irec(i),iinc)= uy0(i,iinc)+uydiff;
                    
                    j=1;
                    if sal.sxy==1
                        sdiff          = sum(DWN.sdiff(1,jj,i).*drj.*phi_fv(jj,iinc).');
                        sw(irec(i),:,j)= s0(j,i,iinc)+sdiff;
                        j=j+1;
                    end
                    if sal.syz==1
                        sdiff          = sum(DWN.sdiff(j,jj,i).*drj.*phi_fv(jj,iinc).');
                        sw(irec(i),:,j)=s0(j,i,iinc)+sdiff;
                    end
                end
                %campo total
            end
        else
            %IBEM-multi 
            for i=1:nrecm
                %campo difractado
                [uydiff,sdiff]  = int_fv_G22_SH(phi_fv,coord,para,m,xr(irec(i)),zr(irec(i)),gaussian);
                %campo total
                j=1;
                if sal.Ut==1
                    uy(irec(i),:,j)=uy0(i,:)+uydiff(j,:);
                    j=j+1;
                end
                if sal.USh==1
                    uy(irec(i),:,j)=uy0(i,:)+uydiff(j,:);
                    j=j+1;
                end
                if sal.UIh==1
                    uy(irec(i),:,j)= uydiff(j,:);
                end
                
                j=1;
                if sal.sxy==1
                    sw(irec(i),:,j)=squeeze(s0(j,i,:))+sdiff(j,:).';
                    j=j+1;
                end
                 if sal.syz==1
                    sw(irec(i),:,j)=squeeze(s0(j,i,:))+sdiff(j,:).';
                end
            end
        end
    end
end