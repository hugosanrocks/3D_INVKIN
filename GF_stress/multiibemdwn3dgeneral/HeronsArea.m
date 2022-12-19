function [area] = HeronsArea(t)
% Area de un triángulo dados tres puntos usando fórmula de Hero de Alexandria
% Entra t de [3x3] con 3 vértices (columnas) con coordenadas x,y,z (renglones). 

% lados del tríangulo al cuadrado
% a2 = (sum((t(1:3,1) - t(1:3,2)).^2)); % v1 a v2
% b2 = (sum((t(1:3,2) - t(1:3,3)).^2)); % v2 a v3
% c2 = (sum((t(1:3,3) - t(1:3,1)).^2)); % v3 a v1
% area = 0.25 * sqrt(4*a2*b2 - (a2+b2-c2)^2);

% o con el producto vectorial
% A = 1/2 |(v1 - v3) x (v2 - v3)|
area = 0.5 * norm(cross(t(1:3,1) - t(1:3,3),...
                        t(1:3,2) - t(1:3,3)));
% error(['test Heron=',num2str(area),' and cross=',num2str(A)])
end
