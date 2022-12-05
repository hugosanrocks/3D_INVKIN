clear all
close all
clc

jetn = load('jetneg.ascii');
jetn = jetn./255;

minx = 0;
maxx = 50;
miny = 0;
maxy = 20;
minz = 0;
maxz = 30;

angdip = 75;
raddip = angdip*pi/180;
ldip = 20;
depth = 10;

z = ldip*sin(raddip);
y = -ldip*cos(raddip);



corns = [0, 0, depth; 0, y, depth+z];

dline = (corns(2,:)-corns(1,:))./10;

plot3(corns(:,1), corns(:,2), corns(:,3),'color',jetn(17,1:3)),hold on
set(gca,'ZDir','reverse');
%set(gca,'YDir','reverse');

lines=15;
sumx = [2, 0, 0];
l(1,:) = corns(1,:) + sumx;
l(2,:) = corns(2,:) + sumx;
for i=1:lines
 plot3(l(:,1), l(:,2), l(:,3),'color',jetn(17,1:3));
 l(1,:) = l(1,:) + sumx;
 l(2,:) = l(2,:) + sumx;
end

corns2 = [0 0 depth; 30 0 depth];
plot3(corns2(:,1), corns2(:,2), corns2(:,3),'color',jetn(17,1:3)),hold on

lh(1,:) = corns2(1,:) + dline;
lh(2,:) = corns2(2,:) + dline;
linesh = 9;
for i=1:linesh
 plot3(lh(:,1), lh(:,2), lh(:,3),'color',jetn(17,1:3));
 lh(1,:) = lh(1,:) + dline;
 lh(2,:) = lh(2,:) + dline;
end
plot3(lh(:,1), lh(:,2), lh(:,3),'color',jetn(17,1:3));

%stations
sta = [28, -15, 0];

plot3(sta(:,1), sta(:,2), sta(:,3), 'vy','markerfacecolor','y','markersize',20)


%fault points
f1 = [0 0 depth] + dline.*5 + sumx.*3;
f1 = f1 + dline.*0.5 + sumx.*0.5;


f2 = f1 + sumx.*2;
plot3(f1(:,1), f1(:,2), f1(:,3),'pr','markerfacecolor',jetn(17,:),'markersize',20);
plot3(f2(:,1), f2(:,2), f2(:,3),'pr','markerfacecolor',jetn(4,:),'markersize',20);

text(f1(1),f1(2)+1,f1(3)-1.5, ...
     't_1', ...
     'Color', 'black', ...
     'BackgroundColor', 'white', ...
     'HorizontalAlignment', 'Left', ...
     'FontSize',20);
text(f2(1),f2(2)+1,f2(3)-1.5, ...
     't_2', ...
     'Color', 'black', ...
     'BackgroundColor', 'white', ...
     'HorizontalAlignment', 'Left', ...
     'FontSize',20);
text(sta(1)+2,sta(2)+2,sta(3), ...
     'Station', ...
     'Color', 'black', ...
     'BackgroundColor', 'white', ...
     'HorizontalAlignment', 'Left', ...
     'FontSize',20);



xlim([minx,maxx])
ylim([-maxy,0])
zlim([minz,maxz])
box on
grid on
ax = gca;
ax.BoxStyle = 'full';



xlabel('X \rightarrow East (km)')
ylabel('Y \rightarrow South (km)')
zlabel('Depth (km)')

set(gca,'Fontsize',25);

