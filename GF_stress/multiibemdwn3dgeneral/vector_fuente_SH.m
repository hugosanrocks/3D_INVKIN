function [B,DWN]=vector_fuente_SH(para,coord,DWN)

% unpack
xs      = para.xs;
zs      = para.zs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vectores de condiciones  a las fronteras %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbeq    = coord.nbeq;
nbpt    = coord.nbpt;
ju      = coord.ju;
j       = 1:coord.nbpt;

B=zeros(nbeq,para.ninc);

if para.fuente==1 %OP
    kzsigno = para.kzsigno;
    kxk     = para.kxk;
    % incidencia de ondas planas homogeneas e inhomogeneas
    
    m       = para.xzs(1);
    mmat    = para.subm1(m);
    
    if m==1 && para.geo(1)==3 %%&& para.nsubmed==1 sino se hace en el campo_ref_SH_FP_DWN
        mu      = para.reg(1).sub(para.nsubmed).Ci(6,6);
        ksi     = para.reg(1).sub(para.nsubmed).ksi;
    else
        mu      = para.reg(mmat).Ci(6,6);
        ksi     = para.reg(mmat).ksi;
    end
    
    jjx     = coord.indm(m).ind;
    salu    = false(nbpt,1);
    salu(j(jjx))= sum(coord.Xim(j(jjx),:)~=0,2)==2;
    
    tmp     = coord.Xim(jjx,:);
    im      = zeros(1,size(tmp,1));
    for i=1:size(tmp,1)
        im(i) = find(tmp(i,:)~=0, 1, 'first' );
    end
    
    tmp=coord.Xim(salu,:);
    imu=zeros(1,size(tmp,1));
    for i=1:size(tmp,1)
        imu(i) = find(tmp(i,:)~=0, 1, 'first' );
    end
    
    for iinc=1:para.ninc
        if para.geo(1)==3 && m==1
            if para.nsubmed==1
                %pour un semi espace slmt ou full-space
                %la decision de la reflexion se fait a l interieur
                [uy,tn] = campo_ref_SH_OPHeI_SE(xs(iinc),zs(iinc),coord,ksi,kxk(iinc),kzsigno(iinc),0,mu,para);
            else
                [uy,tn,DWN]	= campo_ref_SH_OP_DWN(iinc,coord,para,DWN,kxk(iinc),kzsigno(iinc));
            end
        else
            %pour un semi espace slmt
            %ds le cas IBEM-multi avec strates, le sous programme analyse campo_ref_SH_OPHeI_SE
            %si para.fuenteimagen==1 pour prendre en compte une eventuelle
            %reflexion
            [uy,tn] = campo_ref_SH_OPHeI_SE(xs(iinc),zs(iinc),coord,ksi,kxk(iinc),kzsigno(iinc),0,mu,para);
        end
        
        %esfuerzos
        B(j(jjx),iinc)             = -((m==im ) - (m~=im )).*tn(1,j(jjx));%tn
        %desplazamientos
        B(nbpt+ju(salu),iinc)      = -((m==imu) - (m~=imu)).*uy(1,j(salu));%uy
    end
else %if para.fuente==2 FP
    % fuente puntual 
    for iinc=1:para.ninc
        %medio de la fuente
        m       = para.xzs(iinc);
        mmat    = para.subm1(m);
        
        if m==1 && para.geo(1)==3 %%&& para.nsubmed==1 sino se hace en el campo_ref_SH_FP_DWN
            mu      = para.reg(1).sub(1).Ci(6,6);
            ksi     = para.reg(1).sub(1).ksi;
        else
            mu      = para.reg(mmat).Ci(6,6);
            ksi     = para.reg(mmat).ksi;
        end
        
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
        im=zeros(1,size(tmp,1));
        for i=1:size(tmp,1)
            im(i) = find(tmp(i,:)~=0, 1, 'first' );
        end
        
        tmp=coord.Xim(salu,:);
        imu=zeros(1,size(tmp,1));
        for i=1:size(tmp,1)
            imu(i) = find(tmp(i,:)~=0, 1, 'first' );
        end
        
        if para.geo(1)==3 && m==1
            if para.nsubmed==1
                %pour une semi espace slmt
                [uy,tn]     = campo_ref_SH_FP(xs(iinc), zs(iinc),coord,ksi,mu,salu,jjx,para);
                [uy1,tn1]  	= campo_ref_SH_FP(xs(iinc),-zs(iinc),coord,ksi,mu,salu,jjx,para); %fuente imagen
                uy          = uy+uy1;
                tn          = tn+tn1;
            else
                [uy,tn,DWN]	= campo_ref_SH_FP_DWN(iinc,coord,para,DWN);
            end
        else
            [uy,tn]     = campo_ref_SH_FP(xs(iinc), zs(iinc),coord,ksi,mu,salu,jjx,para);
            if para.geo(1)==2 && para.cont(1).ruggeo==1 && m==1 && para.fuenteimagen==1
                [uy1,tn1]  	= campo_ref_SH_FP(xs(iinc),-zs(iinc),coord,ksi,mu,salu,jjx,para); %fuente imagen
                uy          = uy+uy1;
                tn          = tn+tn1;
            end
        end
        
        %esfuerzos
        B(j(jjx),iinc)        	= -((m==im ) - (m~=im )).*tn(1,j(jjx));%tn
        %desplazamientos
        B(nbpt+ju(salu),iinc)	= -((m==imu) - (m~=imu)).*uy(1,j(salu));%uy
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