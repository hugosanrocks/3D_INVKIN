figure;hold on;
plot(DWN0.k2(indd1(end)-2:indd1(end)+3),real(tmp(indd1(end)-2:indd1(end)+3)),'.-')
plot(k20(end) ,0,'+r')
plot(k20l(end),0,'+r')

DWN0.k2 	= linspace(DWN0.k2(indd1(end)-2),DWN0.k2(indd1(end)+3),1e3);
if  pol==1
    tmp = mode_Love(para,DWN0);
else
    tmp = mode_Rayleigh_2(para,DWN0);
end
tmp(tmp==inf)=0;




plot(DWN0.k2,real(tmp),'k')