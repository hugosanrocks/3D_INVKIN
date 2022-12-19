function plotPatch(vert,alfa)
% Graficar pol�gono en el espacio 
% Se obtiene mejores resultados con tri�ngulos (3 v�rtices)
h=fill3(vert(:,1),vert(:,2),vert(:,3),[0.5 1.0 0.333]);
alpha(h,alfa);
if(alfa < 0.5),set(h,'EdgeColor','none');end
end
