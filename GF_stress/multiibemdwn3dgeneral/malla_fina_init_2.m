function para=malla_fina_init_2(para,h)
% discretization inicial fina de los contornos irregulares e originales
% sin tomar en cuenta los recortes
% se hace la discretizacion para una frecuencia 20 veces la mas grande para
% que no haya problema de ubicacion de las fuentes virtuales cuando se trata
% de frecuencias cercanas a la maxima
% eso permite despues evaluar rapidamente las mallas a mas baja frecuencia
% con solo buscar una fuente virtual cerca de los lugares requeridos

nf  = para.nf;
df  = para.df;
fmax= nf/2*df;

%como no se conoce en detalle los contornos, se identifica la velocidad mas
%baja y no nula par discretizar todos los contornos con respeto a esta
%velocidad
cmin    = zeros(para.nmed,1);
for i=1:para.nmed
    if para.reg(i).rho==0
        cmin(i) = 0;
    else
        cmin(i)	= para.reg(i).bet;
        if para.reg(i).bet==0
            cmin(i)	= para.reg(i).alpha;
        end
    end
end
cmin(cmin==0)   = [];
cmin            = squeeze(cmin);
cmin            = min(cmin);
para.cmin       = cmin;

%longitud de onda minima que hay que discretizar corectamente
lambda  = cmin/fmax;
m = 1; c = 1;
if para.geo(1)==2
    %discretizacion del background
    para.cont(m,c).long = long_contorno_m1(para.cont(m,c),para.geo(m));
    nbpt                = round(para.cont(m,c).long/lambda*para.npplo);
    para.cont(m,c).vec  = malla_contorno_dx_m1(nbpt,para.cont(m,c),para.geo(m));
    zmin=0;
    for m2=2:para.nmed
        if para.geo(m2)==2
            if para.cont(m2,2).za+para.cont(m2,2).h>zmin
                zmin = para.cont(m2,2).za+para.cont(m2,2).h;
            end
        end
    end
    if zmin~=0
        para.cont(m,c).vec.zc 	= para.cont(m,c).vec.zc+zmin;
    end
end

% ciclo sobre los contornos
for m=2:para.nmed
    for c=1:2
        if para.siDesktop
            waitbarSwitch(m/para.nmed,h,['Calculo del contorno ',num2str(c),' del medio ',num2str(m),' de manera fina']);
        else
            disp(['Calculo del contorno ',num2str(c),' del medio ',num2str(m),' de manera fina']);
        end
        if para.geo(m)==1
            para.cont(m,c).long = long_contorno(para.cont(m,c));
            nbpt                = max(round(para.cont(m,c).long/lambda*para.npplo),4);
            %             if m==2 && c==2
            %                 nbpt=nbpt*4;
            %             end
            para.cont(m,c).vec 	= malla_contorno_dx_new(nbpt,para.cont(m,c),para.geo(m));
            if para.cont(m,c).h==0
                para.cont(m,c).vec.zc=para.cont(m,c).vec.zc+(-(c==1)+(c==2))*1e-16;
            end
        elseif para.geo(m)==2
            para.cont(m,c).long     = long_contorno_m1(para.cont(m,c),para.geo(m));
            nbpt                    = max(round(para.cont(m,c).long/lambda*para.npplo),4);
            para.cont(m,c).vec      = malla_contorno_dx_m1(nbpt,para.cont(m,c),para.geo(m));
            para.cont(m,c).vec.zc 	= para.cont(m,c).vec.zc+para.cont(m,c).za+para.cont(m,c).h;
        elseif para.geo(m)==3 
            % parte no plana
            contL       = para.cont(m,c);
            contL.a     = contL.ba;
            contL.ba    = 0.25*contL.a;
            contL.long  = long_contorno(contL);
            nbpt        = max(round(contL.long/lambda*para.npplo),4);
            vecL        = malla_contorno_dx_new(nbpt,contL,1);
            if para.cont(m,c).h==0
                vecL.zc=vecL.zc+(-(c==1)+(c==2))*1e-16;
            end
            nbptL       = length(vecL.xc);

            %parte plana
            contR       = para.cont(m,c);
            contR.a     = contR.a -contR.ba;
            contR.xa    = contR.xa+contR.ba;
            contR.long  = long_contorno_m1(contR,2);
            nbpt        = max(round(contR.long/lambda*para.npplo),4);
            vecR        = malla_contorno_dx_m1(nbpt,contR,2);
            vecR.zc     = vecR.zc+contR.za+contR.h;
            
            para.cont(m,c).vec.xc     = [vecL.xc(1:nbptL/2);vecR.xc];
            para.cont(m,c).vec.zc     = [vecL.zc(1:nbptL/2);vecR.zc];
            para.cont(m,c).vec.r      = [vecL.r(1:nbptL/2) vecL.r(nbptL/2)+0.5*(vecL.r(nbptL/2)-vecL.r(nbptL/2-1))+vecR.r];
            para.cont(m,c).vec.vnx     = [vecL.vnx(1:nbptL/2);vecR.vnx];
            para.cont(m,c).vec.vnz     = [vecL.vnz(1:nbptL/2);vecR.vnz];
            
            para.cont(m,c).long     = contL.long/2+contR.long;
        elseif para.geo(m)==4
            % parte no plana
            contR       = para.cont(m,c);
            contR.xa    = contR.xa+contR.a-2*contR.ba;
            contR.a     = contR.ba;
            contR.ba    = abs(contR.h)/2*contR.a;
            contR.long  = long_contorno(contR);
            nbpt        = max(round(contR.long/lambda*para.npplo),4);
            vecR        = malla_contorno_dx_new(nbpt,contR,1);
            vecR.zc(end)= contR.za;
            if para.cont(m,c).h==0
                vecR.zc = vecR.zc+(-(c==1)+(c==2))*1e-16;
            end
            nbptR       = length(vecR.xc);
            
            %parte plana
            contL       = para.cont(m,c);
            contL.a     = contL.a -contL.ba;
            contL.xa    = contL.xa;
            contL.long  = long_contorno_m1(contL,2);
            nbpt        = max(round(contL.long/lambda*para.npplo),4);
            vecL        = malla_contorno_dx_m1(nbpt,contL,2);
            vecL.zc     = vecL.zc+contL.za+contL.h;
            
            para.cont(m,c).vec.xc     = [vecL.xc;vecR.xc(nbptR/2+1:end)];
            para.cont(m,c).vec.zc     = [vecL.zc;vecR.zc(nbptR/2+1:end)];
            para.cont(m,c).vec.r      = [vecL.r vecL.r(end)+0.5*(vecL.r(end)-vecL.r(end-1))+vecR.r(nbptR/2+1:end)-vecR.r(nbptR/2+1)];
            para.cont(m,c).vec.vnx     = [vecL.vnx;vecR.vnx(nbptR/2+1:end)];
            para.cont(m,c).vec.vnz     = [vecL.vnz;vecR.vnz(nbptR/2+1:end)];
            
            para.cont(m,c).long     = contL.long+contR.long/2;
        end
    end
end