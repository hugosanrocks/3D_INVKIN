function medio=medioencontacto(xs,zs,para)
%la primera prueba en el contexto actual es basada en el abscissa xs
%la secunda prueba es la de "Jordan curve theorem"
%el punto xs,zs pertenece al medio R (? E) si todas las semi recta saliendo de
%xs,zs cruza solamente una vez uno de los contornos
%si cruza dos o cero, no pertenece
%las semi rectas estan definidas por z=th(th).(x-xs)+zs con (x-xs)>0
%solamente se checa que los extremos estan de cada lado

medio=1;
m=para.nmed;
while m>1
    test=0;
    if m==2
        if abs(xs)>1 || zs>para.cont(2,2).h || zs<para.cont(2,1).h
            %de plano, no pertenece a m==2
            test=1;
        end
    else
        if      abs(xs-(para.cont(m,1).xa+para.cont(m,1).a))>para.cont(m,1).a || ...
                abs(zs-para.cont(m,1).za)>para.cont(m,2).h || ...
                abs(zs-para.cont(m,1).za)>abs(para.cont(m,1).h)
            %de plano, no pertenece a m
            test=1;
        end
    end
    if test==0
        %parte derecha
        %         x   = linspace(0,1.01*para.cont(m,1).a,1000);
        if m==2
            zc1 = eq_contour(xs,para.cont(m,1));
            zc2 = eq_contour(xs,para.cont(m,2));
        else
            zc1 = eq_contour(xs-(para.cont(m,1).xa+para.cont(m,1).a),para.cont(m,1))+para.cont(m,1).za;
            zc2 = eq_contour(xs-(para.cont(m,1).xa+para.cont(m,1).a),para.cont(m,2))+para.cont(m,1).za;
        end
        if (zs<=zc2) && (zs>=zc1)
            medio=m;
            return;
        end
    end
    m=m-1;
end
if medio==1 && zs<0
    medio=0;
end