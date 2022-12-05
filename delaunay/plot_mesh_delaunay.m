clear all
close all
clc

%s=load('fort.77');
%d=delaunay(s);
%save('-ascii','nodes.in','s');
%save('-ascii','elements.out','d');

pos1 = load('peruchile.xyz');
ver = load('triang.out');
pos = pos1(:,1:2);

%tolerance for side distance
tolerance = 1;

[nt,~] = size(ver);
[np,~] = size(pos);


figure(1)
%scatter(pos(:,1),pos(:,2),'k'),hold on
k = 1;
for i=1:nt
 p1 = [pos(ver(i,1),1) pos(ver(i,1),2)];
 p2 = [pos(ver(i,2),1) pos(ver(i,2),2)];
 p3 = [pos(ver(i,3),1) pos(ver(i,3),2)];
 p = [p1;p2;p3;p1];
 if (( norm(p1-p2) < tolerance ) && (norm(p1-p3) < tolerance ))
 plot(p(:,1),p(:,2),'k'),hold on
 triang(k,:) = ver(i,:);
 a(i) = norm(p1-p2);
 k = k+1;
 end;
end
set(gca,'FontSize',20)

save('-ascii','triangle_vert.out','triang');
