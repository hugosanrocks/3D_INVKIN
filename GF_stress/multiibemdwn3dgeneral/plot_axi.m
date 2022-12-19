close all

dta=127.88/2047;
ta=0:dta:127.88;
dtd=127.9844/8191;
td=0:dtd:127.9844;

z=load('axi.z0001');
e=load('axi.e0001');
n=load('axi.n0001');

figure(1)
plot(ta,e),%hold on,plot(td,res.utc(:,1,1,1),'.r'),xlim([0,120])
xlabel('time (sec)'),ylabel('displacement (m)')
legend('AXITRA','DWN-IBEM')
title('East-West')
xlim([40,120])

figure(2)
plot(ta,n),%hold on,plot(td,res.utc(:,1,1,2),'.r'),xlim([0,120])
xlabel('time (sec)'),ylabel('displacement (m)')
legend('AXITRA','DWN-IBEM')
title('North-South')
xlim([40,130])

figure(3)
plot(ta,z),%hold on,plot(td,res.utc(:,1,1,3),'.r'),xlim([0,120])
xlabel('time (sec)'),ylabel('displacement (m)')
legend('AXITRA','DWN-IBEM')
title('Up-Down')
xlim([40,130])

