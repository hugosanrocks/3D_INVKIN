clear all
close all
clc

nodein = load('nodes.in');
nodeou = load('nodes.out');
elemin = load('elements.in');
slipin = load('slip.in');

x=0.26:0.1:34.24;
y=0.51:0.1:16.49;

k=1;
for i=1:length(y) 
  for j=1:length(x)
   nodeou(k,1:2) = [x(j), y(i)];
   k = k+1;
  end
end

sizou = k-1;

for i=1:sizou
    xw = nodeou(i,1); yw = nodeou(i,2);
    p = 0;
    j=1;
    while (( p == 0 ) && (j <= 1103))
     triang(1,:)=[nodein(elemin(j,1),1), nodein(elemin(j,1),2)];
     triang(2,:)=[nodein(elemin(j,2),1), nodein(elemin(j,2),2)];
     triang(3,:)=[nodein(elemin(j,3),1), nodein(elemin(j,3),2)];
     p = poin(xw,yw,triang(1,1),triang(1,2),triang(2,1),triang(2,2),triang(3,1),triang(3,2));
     j = j+1;
    end
    id(i) = j-1;
    val = [slipin(elemin(j-1,1)), slipin(elemin(j-1,2)),slipin(elemin(j-1,3))];
    valin(i) = interpol(xw,yw,triang(:,1),triang(:,2),val);
end

save('-ascii','slip_m.out','valin');
