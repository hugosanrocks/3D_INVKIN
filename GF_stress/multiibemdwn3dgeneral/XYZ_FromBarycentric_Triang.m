function [r] = XYZ_FromBarycentric_Triang(V,ABW)
% Cartesian coordinates from Triangle Vertex and Barycentric Coords

% Obtain r=([x,y,z],[1:Ncub],[1:Npuntos]) the cartesian global 
% coordiantes for the barycentric coordinate points in ABW also 
% called area or simplex coordinates.
%
% ABW first two columns are the simplex coordinates alpha and beta
% third comlum is W the cubature weight that get passed along.
% The cartesian coordinates are:
%   x = (1 - a - b)x1 + a x2 + b x3 
%   y = (1 - a - b)y1 + a y2 + b y3
%   z = (1 - a - b)z1 + a z2 + b z3
% where x_i y_i z_i are the i_th vertex coordinates from 
% arbitrary triangle "V" of the form:  
%               [3x3xip]
%                | |  '------ índice de triángulo
%                | '--------- índice de vértice
%                '----------- índice de coordenada
% REF : Quadrature in Two Dimensions.PDF

% simplex coordinates
a = ABW(:,1);   b = ABW(:,2);  g = 1 - a - b;%  W = ABW(:,3);

% La cantidad de puntos
Ncub = size(a,1);
Ntriang = size(V,3);  r = zeros(3,Ncub,Ntriang);
for it = 1:Ntriang
  i = 1;   r(i,:,it) = g*V(i,1,it) + a*V(i,2,it) + b*V(i,3,it); %x
  i = 2;   r(i,:,it) = g*V(i,1,it) + a*V(i,2,it) + b*V(i,3,it); %y
  i = 3;   r(i,:,it) = g*V(i,1,it) + a*V(i,2,it) + b*V(i,3,it); %z
end
end