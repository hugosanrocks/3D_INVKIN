
tmax=127.88;
tsamd=8192;

axi=load('axi.x0001');
tsam=length(axi);
dta=tmax/tsam;
t=0:dta:tmax-dta;
t=t+0.5;

dtd=para.tmax/tsamd;
td=0:dtd:para.tmax-dtd;

plot(td,res.utc(:,1,1,1),'r'),hold on
plot(t,axi(1:tsam).*-1,'b'),xlim([36,127]),legend('DWN-IBEM','AXITRA')
