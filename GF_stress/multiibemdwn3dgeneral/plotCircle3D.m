function plotCircle3D(center,normal,radius,estilo)
% Graficar un disco en el espacio
% se indican centro y normal:  [3x1] 
% además el radio y el estilo: [1x1]

theta=0:0.01:2*pi;
v=null(normal);
% puntos todo alrededor del disco:
points=real(repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta)));
if(estilo == 1) %contornos
    plot3(points(1,:),points(2,:),points(3,:),'k-');
elseif(estilo == 2) %geometria
    h=fill3(points(1,:),points(2,:),points(3,:),[0.5 1.0 0.333]);
    alpha(h,0.9);
    vn(:,1) = center;
    vn(:,2) = center+0.5*radius*normal;
    plot3(vn(1,:),vn(2,:),vn(3,:),'b-');
else %phi IBEM
    h=fill3(points(1,:),points(2,:),points(3,:),[1 1 1]*(1-abs(len(estilo))));
    alpha(h,0.9);
    vn(:,1) = real(center);
    vn(:,2) = center'+real(estilo)*(radius/len(real(estilo)));
    plot3(vn(1,:),vn(2,:),vn(3,:),'b-','LineWidth',2);
    vn(:,1) = center;
    vn(:,2) = center'+imag(estilo)*(radius/len(imag(estilo)));
    plot3(vn(1,:),vn(2,:),vn(3,:),'r-','LineWidth',2);
end
end

function [l] = len(vec)
l = sqrt(vec(1)^2+vec(2)^2+vec(3)^2);
end