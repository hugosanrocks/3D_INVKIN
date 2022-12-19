function para=pos_rec(para)
if para.recpos<3 %En la superfice ó Malla constante
    para.rec.nrec   = para.rec.nrecx*para.rec.nrecy*para.rec.nrecz;
    nrec            = para.rec.nrec;
    
    xri     = para.rec.xri;
    dxr     = para.rec.dxr;
    nrecx   = para.rec.nrecx;
    yri     = para.rec.yri;
    dyr     = para.rec.dyr;
    nrecy   = para.rec.nrecy;
    zri     = para.rec.zri;
    dzr     = para.rec.dzr;
    nrecz   = para.rec.nrecz;
    
    
    xr1     = zeros(nrec,1);
    yr1     = zeros(nrec,1);
    zr1     = zeros(nrec,1);
    medio   = zeros(nrec,1);
    tipo    = zeros(nrec,1);
    
    for iz=1:nrecz
        for iy=1:nrecy
            yr      = yri+dyr*(iy-1);
            for ix=1:nrecx
                xr      = xri+dxr*(ix-1);
                
                if para.dim==1
                    rxy= xr;
                else
                    %axisimetria centrado en 0 solo (por ahora)
                    rxy= -sqrt(xr^2+yr^2);
                end
                
                if para.recpos==1
                    %posicion de los receptores puestos en las superficies
                    %prueba de los z de cada contorno
                    
                    %si ya esta los contornos calculados
                    if isfield(para,'nmed1')
                        zr=1e9*ones(para.nmed1,1);
                        if para.geo(1)==3
                            zr(para.nmed1+1)=0;
                        end
                        for i=para.nmed1:-1:1
                            x0  = para.cont1(i).vec.xc;
                            dx  = sign(diff(x0));
                            %tratar formas no biyectivas
                            indx= find(dx(1:(end-1)).*dx(2:end)<=0)+1;
                            indi= [1; indx];
                            indf= [indx; length(x0)];
                            
                            for j=1:length(indi)
                                if indf(j)==indi(j)+1
                                    zr(i)=para.cont1(i).vec.zc(indi(j));
                                else
                                    x1=x0(indi(j):indf(j));
                                    z1=para.cont1(i).vec.zc(indi(j):indf(j));
                                    if length(indi)==1
                                        if min(x1)<=rxy && max(x1)>=rxy
                                            zr(i) = interp1(x1,z1,rxy);
                                        else
                                            zr(i)=[];
                                            break
                                        end
                                    else
                                        if min(x1)<=rxy && max(x1)>=rxy
                                            if j==1
                                                zr(i) = interp1(x1,z1,rxy);
                                            else
                                                zr(i) = min(interp1(x1,z1,rxy),zr(i));
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        if ~isempty(para.cont1)
                            m1=para.cont1(i).mi;
                            if para.cont1(i).mv==1 && length(zr)>1 && para.reg(m1).rho==0
                                zr(para.nmed1+1)=[];
                            end
                        end
                        zr      = min(zr)+2e-6;
                    else
                        zr=zeros(para.nmed,2);
                        for m=2:para.nmed
                            for c=1:2
                                if para.geo(m)==1 %Dos contornos horizontales
                                    xm      = rxy-para.cont(m,1).xa-para.cont(m,1).a;
                                    if (xm>para.cont(m,1).a) || (xm<-para.cont(m,1).a)
                                        zr(m,c)=1e6;
                                    else
                                        zr(m,c) = eq_contour(xm,para.cont(m,c));
                                        zr(m,c) = zr(m,c)+para.cont(m,1).za;
                                    end
                                elseif para.geo(m)==2 %placa ilimitada
                                    xm      = rxy;
                                    if (xm>para.cont(m,1).a) || (xm<para.cont(m,1).xa)
                                        zr(m,c)=1e6;
                                    else
                                        zr(m,c) = eq_contour_P(xm,para.cont(m,c));
                                        zr(m,c) = zr(m,c)+para.cont(m,1).za;
                                    end
                                elseif para.geo(m)==3 %Semiplaca L
                                    if rxy<para.cont(m,1).xa
                                        zr(m,c) = 1e6;
                                    elseif rxy<(para.cont(m,1).xa+para.cont(m,1).ba)
                                        contL   =para.cont(m,1);
                                        contL.a =contL.ba;
                                        contL.ba=0;
                                        xm  	= rxy-para.cont(m,1).xa-para.cont(m,1).ba;
                                        zr(m,c) = eq_contour(xm,contL);
                                        zr(m,c) = zr(m,c)+para.cont(m,1).za;
                                    else
                                        xm      = rxy;
                                        zr(m,c) = eq_contour_P(xm,para.cont(m,c));
                                        zr(m,c) = zr(m,c)+para.cont(m,1).za;
                                    end
                                elseif para.geo(m)==4 %Semiplaca R
                                    if rxy>para.cont(m,1).xa+para.cont(m,1).a
                                        zr(m,c) = 1e6;
                                    elseif rxy>(para.cont(m,1).xa+para.cont(m,1).a-para.cont(m,1).ba)
                                        contL   =para.cont(m,1);
                                        contL.a =contL.ba;
                                        contL.ba=0;
                                        xm  	= rxy-para.cont(m,1).xa-para.cont(m,1).a+para.cont(m,1).ba;
                                        zr(m,c) = eq_contour(xm,contL);
                                        zr(m,c) = zr(m,c)+para.cont(m,1).za;
                                    else
                                        xm      = rxy;
                                        zr(m,c) = eq_contour_P(xm,para.cont(m,c));
                                        zr(m,c) = zr(m,c)+para.cont(m,1).za;
                                    end
                                end
                            end
                        end
                        % por ahora en 3D no se considera rugosidad
                        %background
                        %semi espacio
                        xm      = rxy;
                        zr(1,1) = eq_contour_m1(xm,para.cont(1,1),para.geo(1));
                        zr(1,2) = zr(1,1);
                        
                        zr      = min(min(zr));
                    end
                    %                     yr      = 0;
                    tipoi = 1;
                else
                    %posicion de los receptores en una malla
                    yr      = yri+dyr*(iy-1);
                    zr      = zri+dzr*(iz-1);
                    tipoi = 2;
                end
                
                ies     = ix+(iy-1)*nrecx+(iz-1)*nrecx*nrecy;
                
                xr1(ies)    = xr;
                yr1(ies)    = yr;
                zr1(ies)    = zr;
                tipo(ies)   = tipoi;
                
                if para.rafraichi==0
                  if para.dim < 4 % 2D,2.5D,3Daxisim
                    if isfield(para.cont,'vec')
                        medio(ies)  = inclusiontest(rxy,zr,para,1);
                    end
%                   elseif para.dim == 4 % 3Dgen
%                     medio(ies) = inclusiontest3G(xr,yr,zr,para);
                  end
                end
            end
        end
    end
    para.rec.xr     = xr1;
    para.rec.yr     = yr1;
    para.rec.zr     = zr1;
    
    if para.dim == 4 % 3Dgen
      % asignar medio a todos los puntos de una vez,
      medio = inclusiontest3G(xr1,yr1,zr1,para);
    end
else
    %entrada libre, dado por el usuario afuera del programa
    %deber proporcionar para.rec.nrec
    %                   para.rec.xr
    %                   para.rec.yr
    %                   para.rec.zr
    
    para.rec.nrec   = para.rec.nrecx;
    para.rec.nrecy  = 1;
    para.rec.nrecz  = 1;
    tipo    = zeros(para.rec.nrec,1);
    tipo(:) = 3;
    
    nrec            = para.rec.nrec;
    medio           = zeros(nrec,1);
    if para.dim>=3
        r = sqrt(para.rec.xr.^2+para.rec.yr.^2);
    else
        r = para.rec.xr;
    end
    
    if para.dim < 4 
    for ies=1:nrec
        if isfield(para.cont(1,1),'vec')
            medio(ies)  = inclusiontest(r(ies),para.rec.zr(ies),para,1);
        end
    end
    elseif para.dim == 4 % 3Dgen
      % asignar medio a todos los puntos de una vez,
      medio = inclusiontest3G(para.rec.xr,para.rec.yr,para.rec.zr,para);
    end   
end
% if para.smallscreen
% agregar receptores en los contornos 2D
para.rec.nresAtBoundary=0;
if para.rec.resatboundary
  if para.dim < 4 % 2D,2.5D,3Daxisim
    cmd_resatbou; %xctmp,zctmp,nrestmp
    para.rec.nresAtBoundary = nrestmp;
    para.rec.nrecx  = para.rec.nrecx + nrestmp;
    para.rec.nrec   = para.rec.nrec + nrestmp;
    nrecOld = nrec;
    nrec            = para.rec.nrec;
    medio           = [medio; zeros(nrestmp,1)];
    tipo            = [tipo;  zeros(nrestmp,1)];
    para.rec.xr = [para.rec.xr;xctmp];
    para.rec.zr = [para.rec.zr;zctmp];
    tipo(nrestmp+1:end) = 4;
    clear xctmp zctmp nrestmp
    % calcular el medio en los receptores en la frontera:
    for ies=nrecOld+1:nrec% 1:nrec %todos
        if isfield(para.cont(1,1),'vec')
            medio(ies)  = inclusiontest(...
                para.rec.xr(ies),...
                para.rec.zr(ies),para,1);
        end
    end
  elseif para.dim == 4 % 3Dgen
    error('Falta implementar receptores en automático en la frontera 3Dgen')
  end
end
% end
para.rec.medio  = medio;
para.rec.tipo   = tipo;
if isfield(para,'nmed0')
    for i=1:para.nmed0
        para.rec.m(i).ind     = find(medio==i);
    end
end
clear medio tipo

%en caso de un medio "fondo" estratificado, se identifica la capa de
%cada uno de los receptores perteneciendo a m=1
if para.geo(1)==3
    if isfield(para,'nmed0')
        indr=para.rec.m(1).ind;
    else
        indr=1:para.rec.nrec;
    end
    zr  = para.rec.zr(indr);
    nr  = length(zr);
    for ir=1:nr
        icr     = 1;
        zricr   = zr(ir); % profundidad relativa a la interface de la capa del receptor
        while zricr>para.reg(1).sub(icr).h && icr<para.nsubmed
            zricr   = zricr-para.reg(1).sub(icr).h;
            icr     = icr+1;
        end
        para.rec.zricr(indr(ir))   = zricr;
        para.rec.icr(indr(ir))     = icr;
    end
end
end
