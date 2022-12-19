function [B,DWN]=vector_fuente3D(para,coord,DWN)

% unpack
xs      = para.xs;
ys      = para.ys;
zs      = para.zs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vectores de condiciones a las fronteras %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbeq    = coord.nbeq;
nbpt    = coord.nbpt;
ju      = coord.ju;
j       = 1:coord.nbpt;

B       = zeros(3*nbeq,para.ninc);

if para.fuente==1 %OP
    % incidencia de ondas planas homogeneas e inhomogeneas
    kxk     = para.kxk; %cartesiana x de la normal de la OP
    kyk     = para.kyk; %cartesiana y de la normal de la OP
    polOP   = para.tipo_onda;
    
    m       = para.xzs(1);   % medio al que pertence la primer fuente
    if isfield(para,'sumb1')
      mmat    = para.subm1(m); % indice correcto de m si hubo subdivición
    else
      mmat    = m; % en el caso 3D gen
    end
    
    Ci      = para.reg(mmat).Ci;
    ksi     = para.reg(mmat).ksi;
    kpi     = para.reg(mmat).kpi;
    kri     = para.reg(mmat).kri;
    
    %indice de los puntos perteneciendo al medio de la fuente, todos
    %estos cumplen con la continuidad de los esfuerzos
    jjx     = coord.indm(m).ind;
    %se busca si el medio de la fuente es el de indice mas pequeno de
    %los medios en contactos en el punto de colocacion para fijar el
    %signo de las ecuaciones
    %(cf vector fuente B y construcion matriz A)
    tmp = coord.Xim(jjx,:);
    im  = zeros(1,size(tmp,1));
    for i=1:size(tmp,1)
        im(i) = find(tmp(i,:)~=0, 1, 'first' );
    end
    
    %indice de los puntos perteneciendo al medio de la fuente y que
    %tienen que cumplir con continuidad de desplazamientos
    salu    = false(nbpt,1);
    salu(j(jjx))= sum(coord.Xim(j(jjx),:)~=0,2)==2;
    %se busca si el medio de la fuente es el de indice mas pequeno de
    %los medios en contactos en el punto de colocacion para fijar el
    %signo de las ecuaciones
    %(cf vector fuente B y construcion matriz A)
    tmp=coord.Xim(salu,:);
    imu=zeros(1,size(tmp,1));
    for i=1:size(tmp,1)
        imu(i) = find(tmp(i,:)~=0, 1, 'first' );
    end
    
    for iinc=1:para.ninc
        if m==1 && para.geo(1)==3 % Medio de fondo & Estratificado
            [u,t,DWN] = campo_ref_OP3D_DWN(iinc,coord,para,DWN,polOP(iinc),kxk(iinc));
        elseif para.geo(1)==2 % Semiespacio ...
            if para.fuenteimagen==1 % ... usando la fuente imagen
                [u,t] = campo_ref_OPHeI_3D_SE(xs(iinc),ys(iinc),zs(iinc),coord,kpi,ksi,kri,Ci,polOP(iinc),kxk(iinc),kyk(iinc));
            else % ... usando puntos de colocación
                [u,t] = campo_ref_OPHeI_3D(xs(iinc),ys(iinc),zs(iinc),coord,kpi,ksi,kri,Ci,polOP(iinc),kxk(iinc),kyk(iinc));
            end
        elseif para.geo(1)==1 % Espacio completo
            [u,t] = campo_ref_OPHeI_3D(xs(iinc),ys(iinc),zs(iinc),coord,kpi,ksi,kri,Ci,polOP(iinc),kxk(iinc),kyk(iinc));            
        end
        %esfuerzos
        B(j(jjx)       ,iinc)       = -((m==im ) - (m~=im )).*t(1,j(jjx));%tx
        B(j(jjx)+  nbeq,iinc)       = -((m==im ) - (m~=im )).*t(2,j(jjx));%ty
        B(j(jjx)+2*nbeq,iinc)    	= -((m==im ) - (m~=im )).*t(3,j(jjx));%tz
        %desplazamientos
        B(nbpt+ju(salu)       ,iinc)= -((m==imu) - (m~=imu)).*u(1,j(salu));%ux
        B(nbpt+ju(salu)+  nbeq,iinc)= -((m==imu) - (m~=imu)).*u(2,j(salu));%uy
        B(nbpt+ju(salu)+2*nbeq,iinc)= -((m==imu) - (m~=imu)).*u(3,j(salu));%uz
    end
else %if para.fuente==2 FP
    fij     = para.fij;
    % fuente puntual
    for iinc=1:para.ninc
        m       = para.xzs(iinc);
        mmat    = para.subm1(m);
        Ci      = para.reg(mmat).Ci;
        ksi     = para.reg(mmat).ksi;
        kpi     = para.reg(mmat).kpi;
        
        %indice de los puntos perteneciendo al medio de la fuente, todos
        %estos cumplen con la continuidad de los esfuerzos
        jjx     = coord.indm(m).ind;
        %indice de los puntos perteneciendo al medio de la fuente y que
        %tienen que cumplir con continuidad de desplazamientos
        salu    = false(nbpt,1);
        salu(j(jjx))= sum(coord.Xim(j(jjx),:)~=0,2)==2;
        
        %se busca si el medio de la fuente es el de indice mas pequeno de
        %los medios en contactos en el punto de colocacion para fijar el
        %signo de las ecuaciones
        %(cf vector fuente B y construcion matriz A)
        tmp=coord.Xim(jjx,:);
        im =zeros(1,size(tmp,1));
        for i=1:size(tmp,1)
            im(i) = find(tmp(i,:)~=0, 1, 'first' );
        end
        
        tmp=coord.Xim(salu,:);
        imu=zeros(1,size(tmp,1));
        for i=1:size(tmp,1)
            imu(i) = find(tmp(i,:)~=0, 1, 'first' );
        end
        
        
        %se calcula los terminos fuentes
        if m==1 && para.geo(1)==3 % [la primer fuente corresponde al medio 1 (fondo)] && [es fondo estratificado]
            [u,t,DWN]   = campo_ref_FP3D_DWN(iinc,coord,para,DWN);
        else
            [u,t]       = campo_ref_FP3D(xs(iinc),ys(iinc),zs(iinc),coord,kpi,ksi,Ci,fij(iinc,:),salu,jjx,para);
        end
        
        %esfuerzos
        B(j(jjx)       ,iinc)      	= -((m==im ) - (m~=im )).*t(1,j(jjx));%tx
        B(j(jjx)+  nbeq,iinc)       = -((m==im ) - (m~=im )).*t(2,j(jjx));%ty
        B(j(jjx)+2*nbeq,iinc)      	= -((m==im ) - (m~=im )).*t(3,j(jjx));%tz
        
        %desplazamientos
        B(nbpt+ju(salu)       ,iinc)= -((m==imu) - (m~=imu)).*u(1,j(salu));%ux
        B(nbpt+ju(salu)+  nbeq,iinc)= -((m==imu) - (m~=imu)).*u(2,j(salu));%uz
        B(nbpt+ju(salu)+2*nbeq,iinc)= -((m==imu) - (m~=imu)).*u(3,j(salu));%uz
        % en un punto de colocacion, hay siempre 2 medios en contacto,
        % aunque a veces este medio es el vacio
        % las ecuaciones de continuidad de traccion o de desplazamiento
        % involven entonces siempre 2 medios  : sigma(m)=sigma(m1) o u(m)=u(m1)
        % y como cada campo corresponde a (por ejemplo) u=u_diff+u_inc
        % se reorganisa entonces como :
        % sigma_diff(m)-sigma_diff(m1)=u_inc(m1)-u_inc(m)
        % los signos se fijan deacuerdo con el indice del medio,
        % + si el indice es el mas pequeno de los medios en contacto
        % - en lo contrario
    end
end