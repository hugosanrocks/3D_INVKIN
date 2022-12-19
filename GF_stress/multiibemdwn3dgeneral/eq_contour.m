function z=eq_contour(x,cont)
% fonction donnant z(x) pour la geometrie consideree

h	=cont.h;
geom=cont.geom;
a   =cont.a;

ba 	=cont.ba/a;
x   =x/a;

x2  =x.^2;

%1 :LS
%2 :parabola
%3 :triangulo
%4 :coseno
%5 :eliptico
%6 :eliptico asimetrico
%7 :trapeze
%8 :LG

if geom==1
    z=h*(1-x2)./exp(3*x2);
elseif geom==2
    z=h*(1-x2);
elseif geom==3
    z=h*(1-abs(x));
elseif geom==4
    z=h*cos(pi*x/2).^2;
elseif geom==5
    z=h*sqrt(1-x2);
elseif geom==6
    %     if x<=ba
    %         z=h*(1-((x-ba)/(1.+ba)).^2);
    %     else %if x>ba
    %         z=h*(1-((x-ba)/(1.-ba)).^2);
    %     end
    signo=sign(ba-x)+(x==ba);
    z=h*(1-((x-ba)./(1+signo*ba)).^2);
elseif geom==7
    d=0.05/a;%courbure trapeze
    p=-h/(4*d*(1-ba));
    if ba==1
        p=0;
    end
    aux1=abs(x)-(ba-d);
    aux2=abs(x)-(ba+d);
    z=h+p*(aux1.^2.*(aux1>0)-aux2.^2.*(aux2>0));
elseif geom==8
    %     if x<=(ba-1)
    %         z=sqrt(ba^2-(x+1-ba).^2);
    %     elseif x>(ba-1) && x<=(1-ba)
    %         z=ba;
    %     else %if x>(1-ba)
    %         z=sqrt(ba*ba-(x-1.+ba).^2);
    %     end
    z=sqrt(ba^2-(x+1-ba).^2)    .*(x<=(ba-1)) +...
        ba                      .*(x>(ba-1)).*(x<=(1-ba)) +...
        sqrt(ba^2-(x-1+ba).^2)  .*(x>(1-ba));
    z=z/ba*h;
end
z=real(z);
z(x==-1)=1e-16*sign(h);
z(x== 1)=1e-16*sign(h);
