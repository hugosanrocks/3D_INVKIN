function long=long_contorno(cont)
% funcion que calcula la longitud total de los contornos
% cual corresponden a las fronteras comun entre los medios
% = perimetro

%1 :tipo gaussiana
%2 :parabola
%3 :triangulo
%4 :coseno
%5 :eliptico
%6 :eliptico asimetrico
%7 :trapecio
%8 :base puesto en los bordes

h       =cont.h;
geom    =cont.geom;
a       =cont.a;

% if geom==2
%     sq  =sqrt(1.0+4.0*h^2);
%     if h~=0
%         long=sq+0.5/abs(h)*log(2.0*abs(h)+sq);
%     else
%         long=2*sq;
%     end
% else
if geom==3
    long=2.0*sqrt(h^2+a^2);
% elseif geom==5
%     long=pi*0.5*(1.5*(1.0+abs(h))-sqrt(abs(h)));
% elseif geom==8
%     long=pi*ba+2.*(1.0-ba);
else
    %caso general
    %perimetro calculado con el integral de la valor absoluta de la derivada
    long=quad(@(x) local_long_contorno(x,cont),-a,a,1e-6);
end