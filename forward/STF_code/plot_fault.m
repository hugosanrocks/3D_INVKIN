clear all
close all

fault2d=load('fault_2d.out');
fault3d=load('fault_3d.out');
stat=load('stations_2d.out');

epi=[0 0];

sname=importdata('stations.name');

nstat=length(stat);
nsubf=length(fault2d);

stat=stat/1e3;
fault2d=fault2d./1e3;
fault3d=fault3d./1e3;

xmin=min(stat(:,1));
xmax=max(stat(:,1));
ymin=min(stat(:,2));
ymax=max(stat(:,2));

figure(1)

% Plot fault
for isubf=1:nsubf
    plot(fault2d(isubf,1),fault2d(isubf,2),'ob'); hold on
end

% Plot stations
for istat=1:nstat
    plot(stat(istat,1),stat(istat,2),'or','LineWidth',2); hold on
    text(stat(istat,1)-9,stat(istat,2)+8,sname(istat),'fontsize',10); hold on
end

% Plot epicenter
plot(epi(1),epi(2),'ok','LineWidth',2); hold on

axis image
axis([xmin-30 xmax+30 ymin-30 ymax+30])

xlabel('Longitude UTM (km)')
ylabel('Latitude UTM (km)')

filename2 = ['fault_plot_2d.eps'];    
print ('-depsc2' , filename2);


