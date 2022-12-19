function [phi_fv,coord,DWN]=sol_dfv(para,fj,DWN)
% % calcula la densidad de fuerzas virtuales (dfv)
% % funcion que permite construir la malla en funcion de la frecuencia y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discretizacion de los contornos %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coord   = malla_geom(para,fj);
nbeq    = coord.nbeq;
nbpt    = coord.nbpt;
ju      = coord.ju;
sal     = para.sortie;

if para.pol==1
    nbeq2   = nbeq;
elseif para.pol==2
    %hay una ecuacion para phix y otra para phiz
    nbeq2   = 2*nbeq-coord.nfl;
end
coord.nbeq2 = nbeq2;

% se calcula la matriz del DWN
% y se busca los puntos en los cuales hay que calcular G^DWN
% acceleracion de los calculos del DWN para rellenar la matriz del IBEM
if para.geo(1)==3
    DWN.omegac      = 2*pi*fj;
    
    indpt           = 1:coord.nbpt;
    indi            = coord.indm(1,1).ind;%true-false de los puntos perteneciendo a m=1
    indir           = indpt(indi);%indice de estos puntos
    
    if para.pol==2 || (para.nsubmed>1 && para.pol==1  )
        %con nsubmed==1 se hace con una fuente imagen
        
        %se calcula el vector de valores de la componente horizontal de k
        %         xmax        = para.DX*para.DWNnbptkx/2;
        nk          = para.DWNnbptkx;
        nk0         = nk;
        DK          = 2*para.DWNkmax/(nk0);
        %xmax doit etre > vmax*tmax=alpha/dfe
        k2          = (0:(fix(nk0/2)))*DK;
        k2(1)       = k2(2)/1000;
        DWN.k2      = k2;
        DWN.dk2     = k2*0+DK;
        
        %         k2          = logspace(0,(fix(nk0/2))*DK,fix(nk0/2)+1);
        %         k2(1)       = k2(2)/1000;
        %         DWN.k2      = k2;
        %         DWN.dk2     = k2*0+DK;
        %
        %         k2=DWN.k2;          %  DWN.k2=k2;
        
        % calculo de la matriz del DWN
        % y inicialisacion del vector de los diffractados
        
        if  para.pol==1
            DWN = calcul_A_DWN_SH_Ncapas_HS(para,DWN);
            DWN = rebuildk2(DWN);
            DWN = calcul_A_DWN_SH_Ncapas_HS(para,DWN);
            %                         DWN = rebuildk2(DWN);
            
            %ns          = (sal.Ux + sal.Uy + sal.Uz)*sal.Ut;
            ns          = sal.Ut;
            DWN.uydiff  = zeros(nbeq,para.rec.m(1).nr);
            nss         = sal.sxy + sal.syz;
            DWN.sdiff   = zeros(nss,nbeq,para.rec.m(1).nr);
        elseif para.pol==2
            DWN = calcul_A_DWN_PSV_Ncapas_HS(para,DWN);
            DWN = rebuildk2(DWN);
            DWN = calcul_A_DWN_PSV_Ncapas_HS(para,DWN);
            
            %             nk2=length(DWN.k2);
            %             tmp=zeros(nk2,1);
            %             for i=1:nk2;
            %                 tmp(i)=abs(det(DWN.A_DWN(:,:,i)));
            %             end;
            %             figure(205);hold on
            %             plot(DWN.k2,tmp,'-o')
            %
            ns         = sal.Ut;
            DWN.uxzdiff=zeros(2,nbeq2,para.rec.m(1).nr);
            nss         = sal.sxx + sal.sxz + sal.szz;
        end
        
        %para el calculo de la matriz del IBEM:
        %posicion puntos de colocacion: receptores y fuentes (virtuales)
        xrv                 = coord.x(indi);
        zrv                 = coord.z(indi);
        nxrv                = length(xrv);
        vn(1,1:nxrv)        = coord.vnx(indi);
        vn(2,1:nxrv)        = coord.vnz(indi);
        salu                = ones(nxrv,1);
        sals                = ones(nxrv,1);
        
        %reduccion du numero de profundidades de las funtes virtuales
        [zrfv,izrfv,zrref]	= pseudo_unique(zrv,para);
        %         zrfv=zrv.';zrref=zrv;izrfv=(1:length(zrv)).';
        coord.z(indi)       = zrref;
        
        %se incluye de una ves los receptores reales para el calculo de los
        %campos difractados, el campo difractado es indepediente de las
        %incidiencias hasta tomar en cuenta los phi cf inversion _SH-PSV_k
        %posicion receptor reales
        xrr                 = para.rec.xr(para.rec.m(1).ind).';
        zrr                 = para.rec.zr(para.rec.m(1).ind).';
        nrr                 = length(xrr);
        [zrr,izr0,~]        = pseudo_unique(zrr,para);
        %         zrr=zrr.';izr0=(1:length(zrr)).';
        
        %concatenacion de los puntos y de las salidas
        xr                  = [xrv,xrr];
        zr0                 = [zrfv;zrr];
        izr0                = [izrfv;izr0+length(zrfv)];
        
        salu                = [salu;ones(nrr,1)*(ns>0)];
        sals                = [sals;ones(nrr,1)*(nss>0)];
        
        %reduccion du numero de profundidades de los receptores
        DWN.zr0             = zr0;
        DWN.izr0            = izr0;
        DWN.xr              = xr;
        DWN.salu            = salu;
        DWN.sals            = sals;
        DWN.nxrv            = nxrv;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculo de los terminos independientes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% termino fuente conocido %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% esfuerzos nulos en la superficie libre
% continuidad de los esfuerzos y desplazamientos normales
% en los contornos fuera de la superficie libre

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vectores de condiciones  a las fronteras %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if para.pol==1
    % incidencia de ondas heterogeneas
    % el numero de incidencias depende de la frecuencia
    if para.fuente==1
        if max(abs(para.kxk))>1
            w           = DWN.omegac;
            tipoatts=zeros(para.nsubmed,1);
            for ms=1:para.nsubmed
                tipoatts(ms)=para.reg(1).sub(ms).tipoatts;
                para.reg(1).sub(ms).tipoatts=0;
            end
            DWNomei     = para.DWNomei;
            para.DWNomei= 0;
            [k20ok,dkdw]= dispersion_curve_k2_fw(para,DWN,real(w));
            %         for ms=1:para.nsubmed
            %             para.reg(1).sub(ms).tipoatts=0;
            %         end
            k20ok(k20ok==0)=[];
            %         figure(205);plot(k20ok,w,'+');
            
            para.ninc   = length(k20ok);
            para.kxk    = real(k20ok./para.reg(1).sub(para.nsubmed).ksi);
            para.kzsigno= ones(para.ninc,1)*para.kzsigno(1);
            para.xs     = ones(para.ninc,1)*para.xs(1);
            para.zs     = ones(para.ninc,1)*para.zs(1);
            para.dkdw   = dkdw;
        end
    end
    
    [B,DWN]	= vector_fuente_SH(para,coord,DWN);
    if para.fuente==1
        if max(abs(para.kxk))>1
            for ms=1:para.nsubmed
                para.reg(1).sub(ms).tipoatts=tipoatts(ms);
            end
            para.DWNomei=DWNomei;
        end
    end
elseif para.pol==2
    [B,DWN]	= vector_fuente_PSV(para,coord,DWN);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% matriz de coeficientes %
%%%%%%%%%%%%%%%%%%%%%%%%%%

A   = zeros(nbeq2,nbeq2);
%en el caso del IBEM-DWN, se rellena solo una parte de la matrice, la que
%corresponde al interior de los contornos diferente del multi-estratos
for i=1:coord.nbpt
    %si Xi pertenece a un solo medio, entonces solo se considera la
    %ecuacion de esfuerzos en superficie libre
    %sino se considera la continuidad de los esfuerzos normales y de los
    %desplazamientos
    if para.pol==1
        A(i,:)=Aij_Tn_SH(i,coord,para);
        if sum(coord.Xim(i,:)~=0)==2
            A(ju(i)+nbpt,:)=Aij_Gn_SH(i,coord,para);
        end
    elseif para.pol==2
        [Aijx,Aijz]=Aij_Tn_PSV(i,coord,para);
        %if coord.fl(i)==0
        % 2 solidos en contactco
        %elseif coord.fl(i)==1
        % 1 solido y 1 liquido en contacto
        % tx es en realidad tn y tz es tt (normal y tangente)
        % ux es en realidad un (normal) y ut =0
        %elseif coord.fl(i)==2
        % 2 fluidos en contacto
        % tx es en realidad tn y tt esta nula
        % ux es en realidad un (normal) y ut =0
                
        %continuidad tx
        A(i,:)        = Aijx;
        %continuidad tz
        if ~(sum(coord.fl(i,:))==2) && ~(sum(coord.Xim(i,:)~=0)==1 && sum(coord.fl(i,:))==1)
            %no cuando es fluido-fluido o una frontera libre de fluido
            A(coord.iectz(i),:)   = Aijz;
        end
        if sum(coord.Xim(i,:)~=0)==2
            %es una interfase entre 2 medios
            [Aijx,Aijz]             = Aij_Gn_PSV(i,coord,para);
            %continuidad ux
            A(ju(i)+nbpt,:)         = Aijx;
            if sum(coord.fl(i,:))==0
                %continuidad uz
                A(coord.iecuz(i),:)= Aijz;
            end
        end
    end
end

if para.geo(1)==3
    %en el caso del IBEM-DWN, es mas ventajoso, rellenar la matrice
    %considerando los phi, es a decir, columna por columna
    %tambien por cada punto de colocacion perteneciendo al multi-estratos
    %hay que considerar una continuidad de traciones y desplazamientos
    
    %     cada fuente esta considerada por separada, la integracion se hace
    %     bien, se puede mejorar por fuentes de misma orientation, mismo z,
    %     y misma longitud
    for j=indir %ciclo sobre las fuentes virtuales
        if para.pol==1 %SH
            if para.nsubmed==1
                %en este caso se ocupa las fuentes virtuales y sus imagenes y no el DWN
                [tn,gn22]                            = Aij_Tn_Gn_SH_HS(j,coord,para);
                A(indir              ,coord.phi(j,1))= tn;
                A(ju(indir)+nbpt 	 ,coord.phi(j,1))= gn22;
            else
                %en este caso se ocupa el DWN
                coordf.xs   = coord.x(j);
                coordf.zs   = zrref(j==indir);
                coordf.vnx  = coord.vnx(j);
                coordf.vnz  = coord.vnz(j);
                coordf.dr   = coord.dr(j);
                
                [U,S,~]       = calcul_US_DWN_SH_Ncapas_HS(para,xr,zr0,izr0,salu,sals,coordf,DWN);
                
                if nrr>0
                    DWN.uydiff(coord.phi(j,1),:) = U(nxrv+1:nxrv+nrr);
                    if sals(nxrv+1)
                        S0  = squeeze(reshape(S(DWN.inds,nxrv+1:nxrv+nrr,:),nss,nrr,1));
                        DWN.sdiff(:,coord.phi(j,1),:) = S0;
                    end
                end
                U               = U(1:nxrv);
                S               = S(:,1:nxrv);
                tn              = S(1,:).*vn(1,:)+S(2,:).*vn(2,:);
                ind_superpos    = logical((xrv==coordf.xs).*(zrref==coordf.zs));
                
                signo1          = (coord.Xim(indir(ind_superpos),1)==1) - (coord.Xim(indir(ind_superpos),1)==2);
                tn(logical(ind_superpos))=signo1*0.5/coord.dr(j);
                
                A(indir              ,coord.phi(j,1))=tn*coord.dr(j);
                A(ju(indir)+nbpt 	 ,coord.phi(j,1))=U *coord.dr(j);
            end
        elseif para.pol==2
            %phi = [0,1] y [1,0]
            
            coordf.xs   = coord.x(j);
            coordf.zs   = zrref(j==indir);
            coordf.vnx  = coord.vnx(j);
            coordf.vnz  = coord.vnz(j);
            coordf.dr   = coord.dr(j);
            ind_superpos= logical((xrv==coordf.xs).*(zrref==coordf.zs));
            signo1      = (coord.Xim(indir(ind_superpos),1)==1) - (coord.Xim(indir(ind_superpos),1)==2);
            
            [U2,S2,U1,S1] = calcul_US_DWN_PSV_Ncapas_HS(para,xr,zr0,izr0,salu,sals,coordf,[1 1],DWN);%1& 2 oppose ds IBEM & DWN
            
            DWN.uxzdiff(:,coord.phi(j,1),:)         = U1(:,nxrv+1:nxrv+nrr);%(2,nbeq2,para.rec.m(1).nr)
            U1         	= squeeze(U1(:,1:nxrv,1));
            S1         	= squeeze(S1(:,:,1:nxrv,1));
            DWN.uxzdiff(:,coord.phi(j,1)+nbeq,:)    = U2(:,nxrv+1:nxrv+nrr);
            U2         	= squeeze(U2(:,1:nxrv,1));
            S2         	= squeeze(S2(:,:,1:nxrv,1));
            
            tn1         = squeeze(S1(1,1,:)).'.*vn(1,:)+squeeze(S1(1,2,:)).'.*vn(2,:);
            tn2         = squeeze(S1(1,2,:)).'.*vn(1,:)+squeeze(S1(2,2,:)).'.*vn(2,:);
            tn1(ind_superpos)=signo1*0.5/coord.dr(j);
            tn2(ind_superpos)=0;
            %continuidad tx
            A(indir              ,coord.phi(j,1))= tn1*coord.dr(j);
            %continuidad tz
            A(indir+nbeq         ,coord.phi(j,1))= tn2*coord.dr(j);
            %continuidad ux
            A(ju(indir)+nbpt 	 ,coord.phi(j,1))= U1(1,:)*coord.dr(j);
            %continuidad uz
            A(ju(indir)+nbpt+nbeq,coord.phi(j,1))= U1(2,:)*coord.dr(j);
            
            tn1  = squeeze(S2(1,1,:)).'.*vn(1,:)+squeeze(S2(1,2,:)).'.*vn(2,:);
            tn2  = squeeze(S2(1,2,:)).'.*vn(1,:)+squeeze(S2(2,2,:)).'.*vn(2,:);
            tn1(ind_superpos)=0;
            tn2(ind_superpos)=signo1*0.5/coord.dr(j);
            %continuidad tx
            A(indir              ,coord.phi(j,1)+nbeq)= tn1*coord.dr(j);
            %continuidad tz
            A(indir+nbeq         ,coord.phi(j,1)+nbeq)= tn2*coord.dr(j);
            %continuidad ux
            A(ju(indir)+nbpt 	 ,coord.phi(j,1)+nbeq)= U2(1,:)*coord.dr(j);
            %continuidad uz
            A(ju(indir)+nbpt+nbeq,coord.phi(j,1)+nbeq)= U2(2,:)*coord.dr(j);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% resolucion del sistema %
%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_fv =zeros(nbeq2,para.ninc);

if para.ninc>1
    A1=inv(A);
    for iinc=1:para.ninc
        phi_fv(:,iinc)=A1*B(:,iinc); %#ok<*MINV>
        %regularizacion de phi_fv
    end
else
    
    phi_fv=A\B;
    %         regularizacion_phi;
    
    %     if issparse(A), R = qr(A);
    %     else R = triu(qr(A)); end
    %     for iinc=1:para.ninc
    %         phi_fv(:,iinc) = R\(R'\(A'*B(:,iinc)));
    %         r = B(:,iinc) - A*phi_fv(:,iinc);
    %         e = R\(R'\(A'*r));
    %         phi_fv(:,iinc) = phi_fv(:,iinc) + e;
    %         regularizacion_phi;
    %     end
end