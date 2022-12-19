function para=normalizacion(para)
% Normalización de variables y constantes para el cálculo
% OP Onda Plana
% FP Fuerza puntual *****************************************
%  ".gamr"   : OP :  gamma en radianes : ángulo desde eje z
%  ".phir"   : OP :  phi en radianes   : ángulo desde el eje x
%  ".kxk"    : OP :  cartesiana x de la normal de la OP
%  ".kyk"    : OP :  cartesiana y de la normal de la OP
%  ".kzsigno": OP :  - downward  + upward
%  ".fij"    : FP :  coord. vector dirección de la fuente 
% parémetros de geometría ***********************************
%  ".a" ".xa" ".za" ".thr" ".th"
% propiedades de los medios *********************************
%  ".C"      : (1,1) = rho alpha^2
%              (6,6) = rho beta^2 = mu
%              (1,2) = rho (alpha^2 - 2 beta^2)

if mod(para.nf,2)==1
    para.nf=para.nf-1;
end

%----------------------%
% conversion en radian %
%----------------------%
para.gamr    = pi/180*para.gam;
if para.fuente==1 %OP
    if para.dim>=3
        para.phir    = pi/180*para.phi;
        para.kxk   	=      sin(para.gamr).*cos(para.phir);
        para.kyk   	=      sin(para.gamr).*sin(para.phir);
        para.kzsigno= sign(cos(para.gamr));
    else
        para.kxk   	=      sin(para.gamr);
        para.kzsigno= sign(cos(para.gamr));
    end
else %FP
    if para.dim==1
        para.fij     = [sin(para.gamr).' -cos(para.gamr).'];
    else
        para.phir    = pi/180*para.phi;
        para.fij     = [(sin(para.gamr).*cos(para.phir)).' ...
            (sin(para.gamr).*sin(para.phir)).' ...
            -cos(para.gamr).'];
    end
end
%-----------------------%
% parametros geometrico %
%-----------------------%
% reescritura de parametros de geometria simetrica
% permite que se herede unas caracticas simetricas no puestas anterioremente
para.cont(1,1).za	= 0;
for m=1:para.nmed
    para.cont(m,2).a	= para.cont(m,1).a ;
    para.cont(m,2).xa	= para.cont(m,1).xa;
    para.cont(m,2).za	= para.cont(m,1).za;
    para.cont(m,1).thr  = para.cont(m,1).th *pi/180;
    para.cont(m,2).th	= para.cont(m,1).th;
end

%--------------------------%
% parametros de los medios %
%--------------------------%
if para.pol == 1 && para.dim == 1
    for i=para.nmed:-1:1
        if i==1 && para.geo(1) ==3
            for j=1:para.nsubmed
                para.reg(i).sub(j).C(6,6)  = para.reg(i).sub(j).rho*para.reg(i).sub(j).bet^2;
            end
        else
            para.reg(i).C(6,6)  = para.reg(i).rho*para.reg(i).bet^2;
        end
    end
else
    for i=para.nmed:-1:1;
         if i==1 && para.geo(1) ==3 
            for j=1:para.nsubmed
                para.reg(i).sub(j).C(1,1)  = para.reg(i).sub(j).rho*para.reg(i).sub(j).alpha^2;
                para.reg(i).sub(j).C(6,6)  = para.reg(i).sub(j).rho*para.reg(i).sub(j).bet^2;
                para.reg(i).sub(j).C(1,2)  = para.reg(i).sub(j).C(1,1) - 2*para.reg(i).sub(j).C(6,6);
            end
        else
            para.reg(i).C(1,1)  = para.reg(i).rho*para.reg(i).alpha^2;
            para.reg(i).C(6,6)  = para.reg(i).rho*para.reg(i).bet^2;
            para.reg(i).C(1,2)  = para.reg(i).C(1,1) - 2*para.reg(i).C(6,6);
        end
    end
end