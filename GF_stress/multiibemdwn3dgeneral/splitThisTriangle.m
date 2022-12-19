function tNuevos = splitThisTriangle(sdD,t)
% Partir un tríangulo desde el lado más largo al vértice opuesto
%  sdD contiene el para.cont(m,1).piece{p}.subdibData
% .triangles [3x3xn]   n : el número de triángulos
%             | '--------- índice de vértice
%             '----------- índice de coordenada

% Encontrar lado más largo, su centro y el índice de vértice opuesto
l(3,1:6) = lado(sdD.triangles(:,1,t),sdD.triangles(:,2,t),1,2);
l(1,1:6) = lado(sdD.triangles(:,2,t),sdD.triangles(:,3,t),2,3);
l(2,1:6) = lado(sdD.triangles(:,1,t),sdD.triangles(:,3,t),1,3);

lM = find(l == max(l(:,1)),1); %indice del lado más grande

% primer triángulo nuevo
tNuevos.triangles(:,1,1) = sdD.triangles(:,lM,t); % vértice opuesto
tNuevos.triangles(:,2,1) = sdD.triangles(:,l(lM,5),t); %primer vértice lado partido
tNuevos.triangles(:,3,1) = l(lM,2:4); %nuevo vértice (mitad de lado partido)

%segundo triángulo neuvo
tNuevos.triangles(:,1,2) = sdD.triangles(:,lM,t); % vértice opuesto
tNuevos.triangles(:,2,2) = l(lM,2:4); %nuevo vértice (mitad de lado partido)
tNuevos.triangles(:,3,2) = sdD.triangles(:,l(lM,6),t); %segundo vértice anterior

% Las demás varialbes .centers .areas .N  ;

tNuevos.centers(:,1) = mean(tNuevos.triangles(:,:,1),2); %primer nuevo centro
tNuevos.centers(:,2) = mean(tNuevos.triangles(:,:,2),2); %segundo neuvo centro

tNuevos.areas(1) = HeronsArea(tNuevos.triangles(:,:,1)); %primer área
tNuevos.areas(2) = HeronsArea(tNuevos.triangles(:,:,2)); %sgunda área

tNuevos.N(1,:) = sdD.N(t,:); % la misma normal para ambos
tNuevos.N(2,:) = sdD.N(t,:); % la misma normal para ambos

tNuevos.cv(1,:) = sdD.cv(t); % 1 o 2 si es contorno de arriba o abajo
tNuevos.cv(2,:) = sdD.cv(t);

tNuevos.radios(1) = sqrt(tNuevos.areas(1)/pi);
tNuevos.radios(2) = sqrt(tNuevos.areas(2)/pi);

tNuevos.minVel = sdD.minVel;

% %% Test me:
% figure; hold on; plot3(sdD.triangles(1,[1 2 3 1],t),...
%                        sdD.triangles(2,[1 2 3 1],t),...
%                        sdD.triangles(3,[1 2 3 1],t),'k-','LineWidth',0.5)
%                  plot3(tNuevos.triangles(1,[1 2 3 1],1),...
%                        tNuevos.triangles(2,[1 2 3 1],1),...
%                        tNuevos.triangles(3,[1 2 3 1],1),'b--','LineWidth',3)
%                  plot3(tNuevos.triangles(1,[1 2 3 1],2),...
%                        tNuevos.triangles(2,[1 2 3 1],2),...
%                        tNuevos.triangles(3,[1 2 3 1],2),'r--','LineWidth',3)
%                  plot3(tNuevos.centers(1,1),...
%                        tNuevos.centers(2,1),...
%                        tNuevos.centers(3,1),'b.')
%                  plot3(tNuevos.centers(1,2),...
%                        tNuevos.centers(2,2),...
%                        tNuevos.centers(3,2),'r.')
%                  view([-44 64]); grid on
end

function out = lado(a,b,l1,l2)
out = zeros(1,6);
out(1) = ( (a(1)-b(1))^2 + ...
           (a(2)-b(2))^2 + ...
           (a(3)-b(3))^2 ); %tamaño de la cara al cuadrado
 % punto medio:
 out(2) = b(1) + (a(1)-b(1))/2;
 out(3) = b(2) + (a(2)-b(2))/2;
 out(4) = b(3) + (a(3)-b(3))/2;
 out(5) = l1;
 out(6) = l2;
end