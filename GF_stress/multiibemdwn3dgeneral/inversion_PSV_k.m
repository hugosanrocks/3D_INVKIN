function [uxz,sxz]=inversion_PSV_k(para,coord,phi_fv,gaussian,DWN)

%%%%%%%%%%%%%
% unpacking %
%%%%%%%%%%%%%

xr      = para.rec.xr;
zr      = para.rec.zr;
nrec    = para.rec.nrec;

xs      = para.xs;
zs      = para.zs;
if para.fuente==1
    kxk     = para.kxk;
    polOP   = para.tipo_onda;
else
    fij     = para.fij;
end

ninc    = para.ninc;

sal     = para.sortie;
ns      = (sal.Ux + sal.Uz)*(sal.Ut+sal.UPh+sal.USh+sal.UIh+sal.UPt+sal.USt);
uxz     = zeros(nrec,ninc,ns);
nbeq    = coord.nbeq;
nss     = (sal.sxx + sal.szz + sal.sxz);
sxz     = zeros(nrec,ninc,nss);%3a dim: sxx, xzz, sxz

% calculo de las soluciones a partir de los coeficientes
for m=1:para.nmed0
    %el receptor esta en el medio m
    irec        = para.rec.m(m).ind;
    nrecm       = length(irec);
    mmat        = para.subm1(m);
    
    if nrecm~=0
        Ci  = para.reg(mmat).Ci;
        ksi = para.reg(mmat).ksi;
        kpi = para.reg(mmat).kpi;
        kri = para.reg(mmat).kri;
        
        coordr.x	= xr(irec).';
        coordr.z	= zr(irec).';
        coordr.nbpt = nrecm;
        
        %%%%%%%%%%%%%%%%%%%%%
        %% Campo incidente %%
        %%%%%%%%%%%%%%%%%%%%%
        
        uxz0 =zeros(2,nrecm,ninc);
        uxz0c=zeros(ns,nrecm,ninc);
        szx0 =zeros(3,nrecm,ninc);
        szx0c =zeros(nss,nrecm,ninc);
        for iinc=1:ninc
            if para.xzs(iinc)==m %fuente en m
                salu=true(nrecm,1);
                salt=true(nrecm,1);%salt=false(nrecm,1); %-- ????????????
                if m==1 && para.geo(1)==3
                    %DWN con varias capas,
                    %calculo ya hecho anteriormente en vector_fuente_XXX
                    %se recupera el valor
                    uxz0c(:,:,iinc)      = DWN.uxz0(:,:,iinc);%(2,para.rec.m(1).nr,para.ninc)
                    szx0c(:,:,iinc)      = DWN.szx0(:,:,iinc);
                elseif para.fuente==1 %OP
                    if para.geo(1)==2 %semiespacio
                        if para.fuenteimagen==1
                            [uxz0(:,:,iinc),~,szx0(:,:,iinc)] = campo_ref_PSV_OPHeI_SE(xs(iinc),zs(iinc),coordr,kpi,ksi,kri,Ci,polOP(iinc),kxk(iinc));
                        else
                            [uxz0(:,:,iinc),~,szx0(:,:,iinc)] = campo_ref_PSV_OPHeI(xs(iinc),zs(iinc),coordr,kpi,ksi,kri,Ci,polOP(iinc),kxk(iinc));
                        end
                    elseif para.geo(1)==1
                        [uxz0(:,:,iinc),~,szx0(:,:,iinc)] = campo_ref_PSV_OPHeI(xs(iinc),zs(iinc),coordr,kpi,ksi,kri,Ci,polOP(iinc),kxk(iinc));
                    end
                elseif para.fuente==2 %FP
                    [uxz0c(:,:,iinc),~,szx0c(:,:,iinc)]  = campo_ref_PSV_FP(xs(iinc),zs(iinc),coordr,kpi,ksi,Ci,fij(iinc,:),salu,salt,para);
                end
            end
        end
        if para.fuente==1 %Ordenar de acuerdo a las variables solicitadas
            uxz0c = correcion_uxz0_OP(uxz0,ns,polOP,para);
            szx0c = correcion_szx0_OP(szx0,ns,polOP,para);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%
        %% Campo difractado %%
        %%%%%%%%%%%%%%%%%%%%%%
        
        % sumar campo difractado DWN a campo incidente
        if m == 1 && para.geo(1)==3
            % medio de fondo estratificado
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
            uxzdiff=zeros(2,ninc);
            szxdiff=xeros(3,ninc);
            %campo difractado
            for i=1:nrecm
                j=1;
                if sal.Ut==1
                    %campo total
                    if sal.Ux==1
                        for iinc = 1:ninc
                            uxzdiff(1,iinc)         = sum(DWN.uxzdiff(1,jj,i).*drj.*phi_fv(jj,iinc).'+DWN.uxzdiff(1,jj+nbeq,i).*drj.*phi_fv(jj+nbeq,iinc).');%(2,nbeq2,para.rec.m(1).nr)
                        end
                        uxz(irec(i),:,j)=             squeeze(uxz0c(1,i,:)).'+uxzdiff(1,:);
                        j=j+1;
                    end
                    if sal.Uz==1
                        for iinc = 1:ninc
                            uxzdiff(2,iinc)         = sum(DWN.uxzdiff(2,jj,i).*drj.*phi_fv(jj,iinc).'+DWN.uxzdiff(2,jj+nbeq,i).*drj.*phi_fv(jj+nbeq,iinc).');%(2,nbeq2,para.rec.m(1).nr)
                        end
                        uxz(irec(i),:,j)=             squeeze(uxz0c(2,i,:)).'+uxzdiff(2,:);
                    end
                end
            end
            for i=1:nrecm %2015-10
                j=1;
                if sal.Ut==1
                    if sal.sxx==1
                        for iinc = 1:ninc
                            szxdiff(1,iinc) = sum(...
                                DWN.szxdiff(1,jj,i)     .*drj.*phi_fv(jj,iinc).'+...
                                DWN.szxdiff(1,jj+nbeq,i).*drj.*phi_fv(jj+nbeq,iinc).');
                        end
                        sxz(irec(i),:,j)=squeeze(szx0c(1,i,:)).'+szxdiff(1,:);
                        j=j+1;
                    end
                    if sal.sxz==1
                        for iinc = 1:ninc
                            szxdiff(2,iinc) = sum(...
                                DWN.szxdiff(2,jj,i)     .*drj.*phi_fv(jj,iinc).'+...
                                DWN.szxdiff(2,jj+nbeq,i).*drj.*phi_fv(jj+nbeq,iinc).');
                        end
                        sxz(irec(i),:,j)=squeeze(szx0c(2,i,:)).'+szxdiff(2,:);
                        j=j+1;
                    end
                    if sal.szz==1
                        for iinc = 1:ninc
                            szxdiff(3,iinc) = sum(...
                                DWN.szxdiff(3,jj,i)     .*drj.*phi_fv(jj,iinc).'+...
                                DWN.szxdiff(3,jj+nbeq,i).*drj.*phi_fv(jj+nbeq,iinc).');
                        end
                        sxz(irec(i),:,j)=squeeze(szx0c(3,i,:)).'+szxdiff(3,:);
                    end
                end
            end
        else
            % medio de fondo (m==1) semi-espace ou fondo ilimitado
            for i=1:nrecm
                %campo difractado
                [uxzdiff,szxdiff] = int_fv_Gij_PSV(phi_fv,coord,para,m,xr(irec(i)),zr(irec(i)),gaussian,kpi,ksi,Ci);
                %campo total
                j=1;
                js=1;
                %escritura compacta en funcion de los campos permitidos
                sal_list={'Ut','UPh','USh','UIh','UPt','USt'};
                for il=1:length(sal_list)
                    tmpsal = getfield(sal,sal_list{il});
                    if tmpsal
                        if sal.Ux==1
                            uxz(irec(i),:,j)= squeeze(uxz0c(j,i,:)).'+uxzdiff(j,:);
                            j=j+1;
                        end
                        if sal.Uz==1
                            uxz(irec(i),:,j)= squeeze(uxz0c(j,i,:)).'+uxzdiff(j,:);
                            j=j+1;
                        end
                        if sal.sxx==1
                            sxz(irec(i),:,js)= squeeze(szx0c(js,i,:)).'+szxdiff(js,:);
                            js=js+1;
                        end
                        if sal.sxz==1
                            sxz(irec(i),:,js)= squeeze(szx0c(js,i,:)).'+szxdiff(js,:);
                            js=js+1;
                        end
                        if sal.szz==1
                            sxz(irec(i),:,js)= squeeze(szx0c(js,i,:)).'+szxdiff(js,:);
                            js=js+1;
                        end
                    end
                end
            end
        end
    end
end