function dz=dzdx_contorno_m1(x,cont,geo)
%funcion que calcula la derivada de las curvas del background en funcion de x

rh       = cont.rh;
ruggeo    = cont.ruggeo;
rba      = cont.rba;
x       = x-cont.xa;

if geo==2
    %semi espacio
    if ruggeo==1 %plano
        dz = 0.*x;
    elseif ruggeo==2 %sinus
        dz =  rh/2*2*pi/rba*cos(2*pi*x/rba);
    elseif ruggeo==3 %triangle
        dz = 2*rh/rba+0.*x;
    end
else
    dz = [];
end