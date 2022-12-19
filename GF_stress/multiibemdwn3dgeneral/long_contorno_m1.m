function long=long_contorno_m1(cont,geo)
% funcion que calcula la longitud total del background

rh     	= cont.rh;
ruggeo	= cont.ruggeo;
a       = cont.a;
rba 	= cont.rba;

if geo==2
    %semi espacio
    if ruggeo==1 %plano
        % dzdx=0;
        long = a;
    elseif ruggeo==2 %sinus
        dl =  @(x) sqrt((rh/2*2*pi/rba*cos(2*pi*x/rba)).^2+1);
        long = integral(dl,0,a,'RelTol',1e-6);
    elseif ruggeo==3 %triangle
        % dzdx = 2*rh/rba;
        long = sqrt((2*rh/rba)^2+1)*a;
    end
else
    long = 0;
end