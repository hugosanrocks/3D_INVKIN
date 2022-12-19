function [B,DWN]=vector_fuente_PSV(para,coord,DWN)

% unpack
xs      = para.xs;
zs      = para.zs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vectores de condiciones a las fronteras %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbeq    = coord.nbeq;
nbeq2   = coord.nbeq2;
nbpt    = coord.nbpt;
ju      = coord.ju;
j       = 1:coord.nbpt;

B       = zeros(nbeq2,para.ninc);

if para.fuente==1 %OP
    % incidencia de ondas planas homogeneas e inhomogeneas
    kxk     = para.kxk;
    polOP   = para.tipo_onda;
    
    m       = para.xzs(1);
    mmat    = para.subm1(m);
    
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
    tmp=coord.Xim(jjx,:);
    im=zeros(1,size(tmp,1));
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
    
    %indice de los puntos perteneciendo al medio de la fuente, todos
    %estos cumplen con la continuidad de los esfuerzos tz
    % 1 solido por lo menos
    jjx2=jjx & coord.ectz;
    %se busca si el medio de la fuente es el de indice mas pequeno de
    %los medios en contactos en el punto de colocacion para fijar el
    %signo de las ecuaciones
    %(cf vector fuente B y construcion matriz A)
    tmp=coord.Xim(jjx2,:);
    imz=zeros(1,size(tmp,1));
    for i=1:size(tmp,1)
        imz(i) = find(tmp(i,:)~=0, 1, 'first' );
    end
    
    %identificacion des interfaces du milieu qui sont en contact avec
    %2 solides
    jjxu2=(sum(coord.fl,2)==0) & salu;
    tmp=coord.Xim(jjxu2,:);
    imuz=zeros(1,size(tmp,1));
    for i=1:size(tmp,1)
        imuz(i) = find(tmp(i,:)~=0, 1, 'first' );
    end
    
    
    for iinc=1:para.ninc
        if m==1 && para.geo(1)==3
            [uxz,txz,DWN] = campo_ref_PSV_OP_DWN(iinc,coord,para,DWN,polOP(iinc),kxk(iinc));
        elseif para.geo(1)==2
            if para.fuenteimagen==1
                [uxz,txz,~] = campo_ref_PSV_OPHeI_SE(xs(iinc),zs(iinc),coord,kpi,ksi,kri,Ci,polOP(iinc),kxk(iinc));
            else
                [uxz,txz,~] = campo_ref_PSV_OPHeI(xs(iinc),zs(iinc),coord,kpi,ksi,kri,Ci,polOP(iinc),kxk(iinc));
            end
        elseif para.geo(1)==1
            [uxz,txz,~] = campo_ref_PSV_OPHeI(xs(iinc),zs(iinc),coord,kpi,ksi,kri,Ci,polOP(iinc),kxk(iinc));
        end
        %         %esfuerzos
        %         B(j(jjx),iinc)             = -((m==im ) - (m~=im )).*txz(1,j(jjx));%tx
        %         B(j(jjx)+nbeq,iinc)        = -((m==im ) - (m~=im )).*txz(2,j(jjx));%tz
        %         %desplazamientos
        %         B(nbpt+ju(salu),iinc)      = -((m==imu) - (m~=imu)).*uxz(1,j(salu));%ux
        %         B(nbpt+ju(salu)+nbeq,iinc) = -((m==imu) - (m~=imu)).*uxz(2,j(salu));%uz
        
        %esfuerzos
        B(j(jjx),iinc)              = -((m==im ) - (m~=im )).*txz(1,j(jjx));%tx
        B(coord.iectz(j(jjx2)),iinc)= -((m==imz) - (m~=imz)).*txz(2,j(jjx2));%tz
        
        %desplazamientos
        B(nbpt+ju(salu),iinc)       = -((m==imu) - (m~=imu)).*uxz(1,j(salu));%ux
        B(coord.iecuz(jjxu2),iinc)  = -((m==imuz) - (m~=imuz)).*uxz(2,j(jjxu2));%uz
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
        %se busca si el medio de la fuente es el de indice mas pequeno de
        %los medios en contactos en el punto de colocacion para fijar el
        %signo de las ecuaciones
        %(cf vector fuente B y construcion matriz A)
        tmp=coord.Xim(jjx,:);
        im=zeros(1,size(tmp,1));
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
        
        %indice de los puntos perteneciendo al medio de la fuente, todos
        %estos cumplen con la continuidad de los esfuerzos tz
        % 1 solido por lo menos
        jjx2=jjx & coord.ectz;
        %se busca si el medio de la fuente es el de indice mas pequeno de
        %los medios en contactos en el punto de colocacion para fijar el
        %signo de las ecuaciones
        %(cf vector fuente B y construcion matriz A)
        tmp=coord.Xim(jjx2,:);
        imz=zeros(1,size(tmp,1));
        for i=1:size(tmp,1)
            imz(i) = find(tmp(i,:)~=0, 1, 'first' );
        end
        
        %identificacion des interfaces du milieu qui sont en contact avec
        %2 solides
        jjxu2=(sum(coord.fl,2)==0) & salu;
        tmp=coord.Xim(jjxu2,:);
        imuz=zeros(1,size(tmp,1));
        for i=1:size(tmp,1)
            imuz(i) = find(tmp(i,:)~=0, 1, 'first' );
        end
        
        
        if m==1 && para.geo(1)==3
            [uxz,txz,DWN] = campo_ref_PSV_FP_DWN(iinc,coord,para,DWN);
        else
            [uxz,txz,~] = campo_ref_PSV_FP(xs(iinc),zs(iinc),coord,kpi,ksi,Ci,fij(iinc,:),salu,jjx,para);
        end
        
        %esfuerzos
        B(j(jjx),iinc)              = -((m==im ) - (m~=im )).*txz(1,j(jjx));%tx
        B(coord.iectz(j(jjx2)),iinc)= -((m==imz) - (m~=imz)).*txz(2,j(jjx2));%tz
        
        %desplazamientos
        B(nbpt+ju(salu),iinc)       = -((m==imu) - (m~=imu)).*uxz(1,j(salu));%ux
        B(coord.iecuz(jjxu2),iinc)  = -((m==imuz) - (m~=imuz)).*uxz(2,j(jjxu2));%uz
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