function [x,z]=visu_courbe(geo,cont,npt)
%para visualizar los contornos

if geo==1
    x   = linspace(-cont.a,cont.a,npt);
    z   = eq_contour(x,cont);
    x   = x+cont.xa+cont.a;
    z   = z+cont.za;
elseif geo==2
    x   = linspace(0,cont.a,npt)+cont.xa;
    z   = eq_contour_P(x,cont);
    z   = z+cont.za;
elseif geo==3
    contL=cont;
    contL.a =contL.ba;
    contL.ba=0;
    if contL.a==0 || contL.h==0
        x   = linspace(0,cont.a,npt)+cont.xa;
        z   = eq_contour_P(x,cont);
        z   = z+cont.za-cont.h;
    else
        nptL= floor(npt*contL.a*2/cont.a);
        if mod(nptL,2)==1
            nptL=nptL+1;
        end
        x   = linspace(-contL.a,contL.a,nptL);
        z   = eq_contour(x,contL);
        x   = x(1:nptL/2)+contL.xa+contL.a;
        z   = z(1:nptL/2)+contL.za;
        
        nptR=npt-nptL/2;
        x1   = linspace(0,cont.a-cont.ba,nptR)+cont.xa+cont.ba;
        z1   = eq_contour_P(x1,cont);
        z1   = z1+cont.za;
        x=[x x1];
        z=[z z1];
    end
elseif geo==4
    contR=cont;
    contR.a =contR.ba;
    contR.ba=0;
    if contR.a==0|| contR.h==0
        x   = linspace(0,cont.a,npt)+cont.xa;
        z   = eq_contour_P(x,cont);
        z   = z+cont.za-cont.h;
    else
        nptR= floor(npt*contR.a*2/cont.a);
        if mod(nptR,2)==1
            nptR=nptR+1;
        end
        x   = linspace(-contR.a,contR.a,nptR);
        z   = eq_contour(x,contR);
        x   = x(nptR/2+1:nptR)+contR.xa+cont.a-contR.a;
        z   = z(nptR/2+1:nptR)+contR.za;
        
        nptL=npt-nptR/2;
        x1   = linspace(0,cont.a-cont.ba,nptL)+cont.xa;
        z1   = eq_contour_P(x1,cont);
        z1   = z1+cont.za;
        x=[x1 x];
        z=[z1 z];
    end
end

% dz  =dzdx_contorno(x,cont);

% figure(1);
% hold on
% plot(x,z)
% plot(x,dz,'r')