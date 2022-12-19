% copie el contenido de los ejes h a una nueva figura
f=figure('Units','normalized','name',...
  'Geometria y contornos','numberTitle','off');
copyobj(h,gca);

% darle formato
if(para.dim>=3)
  set(gca,'zdir','reverse');
  view(3)
  xlabel('');ylabel('');zlabel('')
  grid on
  set(gca,'dataaspectratio',[1 1 1],'Projection','perspective','Box','off');
  light('Position',[2*min(get(gca,'xlim')) ...
                    2*max(get(gca,'ylim')) ...
                    5*min(get(gca,'zlim'))],'Style','infinite');
  axis tight
  set(gca, 'XColor', 'r');set(gca, 'YColor', 'g');set(gca, 'ZColor', 'b')
else
  set(gca,'ydir','reverse');
  view(2)
  xlabel('X');ylabel('Z');zlabel('');
  grid off
  set(gca,'dataaspectratio',[1 1 1],'Projection','orthographic','Box','on');
  axis image
  set(gca, 'XColor', 'k');set(gca, 'YColor', 'k');set(gca, 'ZColor', 'k')
end