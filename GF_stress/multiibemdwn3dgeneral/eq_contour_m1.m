function z=eq_contour_m1(x,cont,geo)
% fonction donnant z(x) pour la geometrie consideree

rh      = cont.rh;
ruggeo  = cont.ruggeo;
rba     = cont.rba;
x       = x-cont.xa;

if geo==2 
    if ruggeo==1 %plano
        z=0*x;
    elseif ruggeo==2 %sinus
        z=rh/2*sin(2*pi*x/rba);
    elseif ruggeo==3 %triangle
        z=2*rh/rba*abs(mod(x,rba)-rba/2)-rh/2;
    end
elseif geo==3
    z=0*x;
else
    z=[];
end