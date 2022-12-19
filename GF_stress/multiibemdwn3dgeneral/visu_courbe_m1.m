function [x,z]=visu_courbe_m1(geo,cont,npt)
%para visualizar los contornos

%background
if geo==2
    %semi espacio
    x = linspace(0,cont.a,npt)+cont.xa;
    z = eq_contour_m1(x,cont,geo);
else %espacio completo
    x = []; z = [];
end