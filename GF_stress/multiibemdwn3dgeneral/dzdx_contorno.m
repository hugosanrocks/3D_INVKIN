function dz=dzdx_contorno(x,cont)
%funcion que calcula la derivada de las curvas de los contornos respecto a x

h 	=cont.h;
geom=cont.geom;
ba 	=cont.ba;
a 	=cont.a;
x   =x/a;
x2  =x.^2;

if geom==1
    dz=-2*h*x.*(4-3*x2)./exp(3*x2);
elseif geom==2
    dz=-2*h*x;
elseif geom==3
    dz=-h*x./abs(x);
elseif geom==4
    dz=-h*cos(pi*x/2).*sin(pi*x/2)*pi;
elseif geom==5
    dz=-h*x./sqrt(1-x2);
elseif geom==6
%     if x<=ba
%         dz=-2*h*(x-ba)/(1+ba)^2;
%     else %if x>ba
%         dz=-2*h*(x-ba)/(1-ba)^2;
%     end
    signo=sign(ba-x)+(x==ba);
    dz=-2*h*(x-ba)./(1+signo*ba).^2;
elseif geom==7
    d=0.05;
    p=-h/(2*d*(1-ba));
    aux1=abs(x)-(ba-d);
    aux2=abs(x)-(ba+d);
    dz=p*sign(x).*(aux1.*(aux1>0)-aux2.*(aux2>0));
elseif geom==8
%     if x<=(ba-1)
%         dz=-(x+1-ba)/sqrt(ba*ba-(x+1-ba)^2);
%     elseif x>(ba-1) && x<=(1-ba)
%         dz=0;
%     else %if x>(1-ba)
%         dz=-(x-1+ba)/sqrt(ba*ba-(x-1+ba)^2);
%     end
    dz= -(x+1-ba)./sqrt(ba^2-(x+1-ba).^2).*(x<=(ba-1)) +...
        -(x-1+ba)./sqrt(ba^2-(x-1+ba).^2).*(x>(1-ba));
    dz=dz/ba*h;   
end
dz=dz/a;
dz(isinf(dz))=10;
dz(isnan(dz))=0;