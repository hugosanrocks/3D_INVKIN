clear all
close all
clc

s=load('source.out');
source=load('source_gauss_rt3p0.src');
source1=load('source_gauss_rt1p0.src');


lx=29;
ly=10;


sav1 = load('mapx.ascii');
sav2 = load('mapy.ascii');


maxa = max(max(s));
maxa = maxa*0.9;

[row,col] = size(source);
tsam = row;

%sp = s;
%for j=1:tsam
%for i=1:col
%  if (sp(j,i) < 0)
%   sp(j,i) = 0;
%  end
%end
%end

xx=1:2:58;
yy=1:2:20;

hypo = [xx(7) yy(7)];


for i=1:row
  n = find(max(source(:,i)) == source(:,i));
  source(1:n-20,i) = 0;
  source(n+20:end,i) = 0;
end

tsam = row;
j = 1;
for i=1:1:tsam
  sr=reshape(s(i,:),lx,ly)';
  src=reshape(source(i,:),lx,ly)';
  src1=reshape(source1(i,:),lx,ly)';
  colormap(jet(128))
%  subplot(2,1,1),pcolor(xx,yy,src),pause(0.4),caxis([0,maxa])
%  set(gca,'FontSize',15)%,xlabel('Along strike (km)'),
%  ylabel('Along dip (km)')
%  title('Circular rupture rt = 3 sec. (t = 5 sec)')
%  set(gca,'YDir','reverse')
%  set(gca,'XDir','reverse')
%  set(gca,'FontSize',20)%,xlabel('Along strike (km)'),

  subplot(2,1,1),pcolor(xx,yy,sr),pause(0.4),caxis([0,maxa]),hold on
  plot(hypo(1),hypo(2),'pr','linewidth',10,'markersize',20)
  for ii=1:length(y)
   for jj=1:length(x)
     if ( sr(ii,jj) > 0.5 )
       display('here')
       quiver(xx(jj),yy(ii),sr(ii,jj)*mapx(ii,jj)/1.5,sr(ii,jj)*mapy(ii,jj)/1.5)
     end
   end
  end
  set(gca,'YDir','reverse')
  set(gca,'XDir','reverse')
  set(gca,'FontSize',20)%,xlabel('Along strike (km)'),
  ylabel('Along dip (km)')
  title('Inverted slip-rate Kumamoto real data')
%  set(gca,'YDir','reverse')

  
%  c=colorbar('Position',[0.9 0.1 0.03 0.8]);
%  set(gca,'FontSize',15),ylabel(c,'Slip-rate (m/s)');
c=colorbar('Position',[0.1 0.44 0.805 0.03]);
set(c,'orientation','horizontal')
set(gca,'FontSize',20),ylabel(c,'Slip-rate (m/s)');
  message = sprintf('movie_%02d',j);
  print(message,'-dpng','-r300')

  close all
  j = j+1;
end
%  sr=reshape(s(60,:),21,9)';
%  pcolor(sr),pause(0.4),caxis([0,5])
find_angle
