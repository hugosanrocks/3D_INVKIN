function [ut,st]=inversion_3D_k(para,coord,phi_fv,gaussian,DWN)

%%%%%%%%%%%%%
% unpacking %
%%%%%%%%%%%%%

% coordenadas de los receptores
xr      = para.rec.xr;
yr      = para.rec.yr;
zr      = para.rec.zr;
nrec    = para.rec.nrec;

% coordenadas de la fuente
xs      = para.xs;
ys      = para.ys;
zs      = para.zs;

if para.fuente==1 % onda plana
    kxk     = para.kxk;
    kyk     = para.kyk;
    polOP   = para.tipo_onda;
else              % fuente puntual
    fij     = para.fij;
end

% Num de inclusiones y de ecuaciones
ninc    = para.ninc;
nbeq    = coord.nbeq;

% variables de salida
sal     = para.sortie;
ns      =(sal.Ux + sal.Uy + sal.Uz);
ut      = zeros(nrec,ninc,ns);
nss     = sal.sxx + sal.syy + sal.szz+sal.sxy+sal.sxz+sal.syz;
st      = zeros(nrec,ninc,nss);

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
        coordr.y	= yr(irec).';
        coordr.z	= zr(irec).';
        coordr.nbpt = nrecm;
        
        %campo incidente
        U0=zeros(3  ,nrecm,ninc);
        S0=zeros(nss,nrecm,ninc);
        for iinc=1:ninc
            salu=logical(ones(nrecm,1)*(para.sortie.Ut));
            sals=logical(ones(nrecm,1)*max([para.sortie.sxx,para.sortie.syy,para.sortie.szz, ...
                para.sortie.sxy,para.sortie.sxz,para.sortie.syz]));
            
            if para.xzs(iinc)==m %fuente en m
                if m==1 && para.geo(1)==3
                    %DWN con varias capas,
                    %calculo ya hecho anteriormente en vector_fuente_3D
                    %se recupera el valor
                    U0(:,:,iinc)      = DWN.U0(:,:,iinc);%(3  ,para.rec.m(1).nr,para.ninc)
                    S0(:,:,iinc)      = DWN.S0(:,:,iinc);%(nss,para.rec.m(1).nr,para.ninc)
                elseif para.fuente==1 %OP
                    if para.geo(1)==2
                        if para.fuenteimagen==1
                            [U0(:,:,iinc),~] = campo_ref_OPHeI_3D_SE(xs(iinc),ys(iinc),zs(iinc),coordr,kpi,ksi,kri,Ci,polOP(iinc),kxk(iinc),kyk(iinc));
                        else
                            [U0(:,:,iinc),~] = campo_ref_OPHeI_3D(xs(iinc),ys(iinc),zs(iinc),coordr,kpi,ksi,kri,Ci,polOP(iinc),kxk(iinc),kyk(iinc));
                        end
                    elseif para.geo(1)==1
                        [U0(:,:,iinc),~] = campo_ref_OPHeI_3D(xs(iinc),ys(iinc),zs(iinc),coordr,kpi,ksi,kri,Ci,polOP(iinc),kxk(iinc),kyk(iinc));
                    end
                elseif para.fuente==2 %FP
                    [U0(:,:,iinc),~,S0(:,:,iinc)] = campo_ref_FP3D(xs(iinc),ys(iinc),zs(iinc),coordr,kpi,ksi,Ci,fij(iinc,:),salu,sals,para);
                end
            end
        end
        
        
        if m == 1 && para.geo(1)==3
            %adjuntar campo difractado a campo incidente DWN
            dA      = coord.dA;
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
            
            dAj = dA(ii);
            %integracion de las contribuciones
            Udiff=zeros(3,ninc);
            Sdiff=zeros(nss,ninc);
            for i=1:nrecm
                j=1;
                if sal.Ut==1
                    if sal.Ux==1
                        for iinc = 1:ninc
                            Udiff(1,iinc) = sum(dAj.*(...
                                DWN.Udiff(1,jj       ,i).*phi_fv(jj       ,iinc).'+...
                                DWN.Udiff(1,jj+  nbeq,i).*phi_fv(jj+  nbeq,iinc).'+...
                                DWN.Udiff(1,jj+2*nbeq,i).*phi_fv(jj+2*nbeq,iinc).'));
                        end
                        ut(irec(i),:,j) = squeeze(U0(1,i,:)).'+Udiff(1,:);
                        j=j+1;
                    end
                    if sal.Uy==1
                        for iinc = 1:ninc
                            Udiff(2,iinc)         = sum(dAj.*(...
                                DWN.Udiff(2,jj       ,i).*phi_fv(jj       ,iinc).'+...
                                DWN.Udiff(2,jj+  nbeq,i).*phi_fv(jj+  nbeq,iinc).'+...
                                DWN.Udiff(2,jj+2*nbeq,i).*phi_fv(jj+2*nbeq,iinc).'));
                        end
                        ut(irec(i),:,j) = squeeze(U0(2,i,:)).'+Udiff(2,:);
                        j=j+1;
                    end
                    if sal.Uz==1
                        for iinc = 1:ninc
                            Udiff(3,iinc)         = sum(dAj.*(...
                                DWN.Udiff(3,jj       ,i).*phi_fv(jj       ,iinc).'+...
                                DWN.Udiff(3,jj+  nbeq,i).*phi_fv(jj+  nbeq,iinc).'+...
                                DWN.Udiff(3,jj+2*nbeq,i).*phi_fv(jj+2*nbeq,iinc).'));
                        end
                        ut(irec(i),:,j) = squeeze(U0(3,i,:)).'+Udiff(3,:);
                    end
                end
                
                %campo total esfuerzos
                j=1;
                if sal.sxx==1
                    for iinc = 1:ninc
                        Sdiff(j,iinc)         = sum(dAj.*(...
                            DWN.Sdiff(j,jj       ,i).*phi_fv(jj       ,iinc).'+...
                            DWN.Sdiff(j,jj+  nbeq,i).*phi_fv(jj+  nbeq,iinc).'+...
                            DWN.Sdiff(j,jj+2*nbeq,i).*phi_fv(jj+2*nbeq,iinc).'));
                    end
                    st(irec(i),:,j) = squeeze(S0(j,i,:)).'+Sdiff(j,:);
                    j=j+1;
                end
                if sal.syy==1
                    for iinc = 1:ninc
                        Sdiff(j,iinc)         = sum(dAj.*(...
                            DWN.Sdiff(j,jj       ,i).*phi_fv(jj       ,iinc).'+...
                            DWN.Sdiff(j,jj+  nbeq,i).*phi_fv(jj+  nbeq,iinc).'+...
                            DWN.Sdiff(j,jj+2*nbeq,i).*phi_fv(jj+2*nbeq,iinc).'));
                    end
                    st(irec(i),:,j) = squeeze(S0(j,i,:)).'+Sdiff(j,:);
                    j=j+1;
                end
                if sal.szz==1
                    for iinc = 1:ninc
                        Sdiff(j,iinc)         = sum(dAj.*(...
                            DWN.Sdiff(j,jj       ,i).*phi_fv(jj       ,iinc).'+...
                            DWN.Sdiff(j,jj+  nbeq,i).*phi_fv(jj+  nbeq,iinc).'+...
                            DWN.Sdiff(j,jj+2*nbeq,i).*phi_fv(jj+2*nbeq,iinc).'));
                    end
                    st(irec(i),:,j) = squeeze(S0(j,i,:)).'+Sdiff(j,:);
                    j=j+1;
                end
                if sal.sxy==1
                    for iinc = 1:ninc
                        Sdiff(j,iinc)         = sum(dAj.*(...
                            DWN.Sdiff(j,jj       ,i).*phi_fv(jj       ,iinc).'+...
                            DWN.Sdiff(j,jj+  nbeq,i).*phi_fv(jj+  nbeq,iinc).'+...
                            DWN.Sdiff(j,jj+2*nbeq,i).*phi_fv(jj+2*nbeq,iinc).'));
                    end
                    st(irec(i),:,j) = squeeze(S0(j,i,:)).'+Sdiff(j,:);
                    j=j+1;
                end
                if sal.sxz==1
                    for iinc = 1:ninc
                        Sdiff(j,iinc)         = sum(dAj.*(...
                            DWN.Sdiff(j,jj       ,i).*phi_fv(jj       ,iinc).'+...
                            DWN.Sdiff(j,jj+  nbeq,i).*phi_fv(jj+  nbeq,iinc).'+...
                            DWN.Sdiff(j,jj+2*nbeq,i).*phi_fv(jj+2*nbeq,iinc).'));
                    end
                    st(irec(i),:,j) = squeeze(S0(j,i,:)).'+Sdiff(j,:);
                    j=j+1;
                end
                if sal.syz==1
                    for iinc = 1:ninc
                        Sdiff(j,iinc)         = sum(dAj.*(...
                            DWN.Sdiff(j,jj       ,i).*phi_fv(jj       ,iinc).'+...
                            DWN.Sdiff(j,jj+  nbeq,i).*phi_fv(jj+  nbeq,iinc).'+...
                            DWN.Sdiff(j,jj+2*nbeq,i).*phi_fv(jj+2*nbeq,iinc).'));
                    end
                    st(irec(i),:,j) = squeeze(S0(j,i,:)).'+Sdiff(j,:);
                end
            end
        else
            for i=1:nrecm
                %campo difractado
                [Udiff,Sdiff]  = int_fv_Gij_3D(phi_fv,coord,para,m,xr(irec(i)),yr(irec(i)),zr(irec(i)),gaussian,kpi,ksi,Ci);
                
                %campo total desplazamientos
                j=1;
                if sal.Ut==1
                    if sal.Ux==1
                        ut(irec(i),:,j) = squeeze(U0(1,i,:)).'+Udiff(j,:);
                        j=j+1;
                    end
                    if sal.Uy==1
                        ut(irec(i),:,j) = squeeze(U0(2,i,:)).'+Udiff(j,:);
                        j=j+1;
                    end
                    if sal.Uz==1
                        ut(irec(i),:,j) = squeeze(U0(3,i,:)).'+Udiff(j,:);
                    end
                end
                
                %campo total esfuerzos
                j=1;
                if sal.sxx==1
                    st(irec(i),:,j) = squeeze(S0(j,i,:)).'+Sdiff(j,:);
                    j=j+1;
                end
                if sal.syy==1
                    st(irec(i),:,j) = squeeze(S0(j,i,:)).'+Sdiff(j,:);
                    j=j+1;
                end
                if sal.szz==1
                    st(irec(i),:,j) = squeeze(S0(j,i,:)).'+Sdiff(j,:);
                    j=j+1;
                end
                if sal.sxy==1
                    st(irec(i),:,j) = squeeze(S0(j,i,:)).'+Sdiff(j,:);
                    j=j+1;
                end
                if sal.sxz==1
                    st(irec(i),:,j) = squeeze(S0(j,i,:)).'+Sdiff(j,:);
                    j=j+1;
                end
                if sal.syz==1
                    st(irec(i),:,j) = squeeze(S0(j,i,:)).'+Sdiff(j,:);
                end
            end
        end
    end
end