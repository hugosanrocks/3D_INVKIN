function plotProfile(ParaRegSub,ShouldUseAlf)
% Dibujar el perfil de estratificación de un semiespacio
%   N    :  El número de capas sin contar el semiespacio
N = size(ParaRegSub,1)-1;
z = zeros(N+1,1); alf = z; bet = z;
%   ParaRegSub : estructura con los siguientes campos
%   z    :  N+1 profundidades de interfaz [0 z1 z2 .. zHalfSpace]
z(1) = 0;
for i=2:N+1
  z(i) = z(i-1) + ParaRegSub(i-1).h;
end
offset = 0.003 * z(min(2,length(z)));
%   alf  :  N+1 velocidades de P
%   bet  :  N+1 velocidades de S
for i=1:N+1
  alf(i) = ParaRegSub(i).alpha;
  bet(i) = ParaRegSub(i).bet;
end
if (~ShouldUseAlf); alf = []; end
maxcol = 0.4; % el más obscuro
mincol = 0.6; % el más claro

col = zeros(N+1,1);
maxVel = max(bet); minVel = min(bet);
if (maxVel == minVel)
mvc = 0;    
else
mvc = (maxcol-mincol)/(maxVel-minVel);
end
for i=1:N+1
    col(i) = mincol + bet(i)*mvc;
end
if ShouldUseAlf
  maxVel = max(maxVel,max(alf));
end
%figure(hfig); axes; set(hfig,'Visible', 'off');
% izq = subplot(1,2,1);
% set(gca,'YDir','reverse');
% cla
% hold on
% for i=1:N
%   h = z(i+1)-z(i);
%   if h > 0
%     rectangle('Position',[0,z(i),1,(z(i+1)-z(i))],'FaceColor',col(i).*[1 1 1])
%   end
% end
% i=N+1;
% rectangle('Position',[0,z(i),1,(1.1*z(i)-z(i))],'FaceColor',col(i).*[1 1 1])
% ylim([z(1) 1.1*z(N+1)])
% set(gca,'XTick',[])
% ylabel('Depth [meters]')
% 
% der = subplot(1,2,2);
set(gca,'YDir','reverse');
cla
hold on
view(2)
axis normal
for i=1:N
    line([bet(i) bet(i)],   [z(i) z(i+1)],  'Color','r','LineWidth',2)
    line([bet(i) bet(i+1)], [z(i+1) z(i+1)],'Color','r','LineWidth',2)
    if ShouldUseAlf
    line([alf(i) alf(i)],   [z(i) z(i+1)]+offset,  'Color','b','LineWidth',2)
    line([alf(i) alf(i+1)], [z(i+1) z(i+1)]+offset,'Color','b','LineWidth',2)
    end
end
i = N;
zend = 1.1*z(i+1); if(zend==0); zend = 1; end
line([bet(i+1) bet(i+1)], [z(i+1) zend],'Color','r','LineWidth',2)
if ShouldUseAlf
line([alf(i+1) alf(i+1)], [z(i+1) zend]+offset,'Color','b','LineWidth',2)
end
xlabel('Wave velocity [m/s]')
ylabel('Depth [meters]')
% set(gca,'YTick',[])
xlim([0.75*minVel 1.1*maxVel])
% ligar los ejes
% linkaxes([izq der],'y');

if (ShouldUseAlf)
p = plot(0,0,'r',0,0,'b');
legend(p,'beta', 'alpha');
else
p = plot(0,0,'r');
legend(p,'beta');
end
ylim([z(1) zend])
hold off
% set(hfig,'Visible', 'on');
