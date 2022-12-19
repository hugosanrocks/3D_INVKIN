function para=Vp2Cij(para)
%conversion de las velocidades de fase en componentes de la matriz de
%Christoffel
%--------------------------%
% parametros de los medios %
%--------------------------%
if para.pol == 1
    for j=1:para.nsubmed
        para.reg(1).sub(j).C(6,6)  = para.reg(1).sub(j).rho*para.reg(1).sub(j).bet^2;
    end
else
    for j=1:para.nsubmed
        para.reg(1).sub(j).C(1,1)  = para.reg(1).sub(j).rho*para.reg(1).sub(j).alpha^2;
        para.reg(1).sub(j).C(6,6)  = para.reg(1).sub(j).rho*para.reg(1).sub(j).bet^2;
        para.reg(1).sub(j).C(1,2)  = para.reg(1).sub(j).C(1,1) -2*para.reg(1).sub(j).C(6,6);
    end
end