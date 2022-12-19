function z=eq_contour_P(x,cont)
% fonction donnant z(x) pour la geometrie consideree

rh      = cont.rh;
th      = cont.th;
ruggeo  = cont.ruggeo;
rba     = cont.rba;
h       = cont.h;
xa      = cont.xa;

if rh>=abs(h)/2
    rh=abs(h)/3;
end
if ruggeo==1 %plano
    z=0*x+sin(th)*(x-xa);
elseif ruggeo==2 %sinus
    z=rh/2*sin(2*pi*x/rba)+sin(th)*(x-xa);
elseif ruggeo==3 %triangle
    z=2*rh/rba*abs(mod(x,rba)-rba/2)-rh/2+sin(th)*(x-xa);
end
z=z+h;