function DWN=rebuildf(DWN)
nw=length(DWN.omegac);
tmp=zeros(nw,1);
for i=1:nw;
    tmp(i)=1./abs(det(DWN.A_DWN(:,:,i)));
end;

omegac  = DWN.omegac;

Dw  = omegac(2)-omegac(1);
% 
% figure(205);hold on
% plot(DWN.omegac,tmp,'-o')

ddtmp       = abs(diff(tmp));
ddtmp       = [0;(ddtmp(1:end-1)+ddtmp(2:end))/2]/mean(ddtmp);
dtmp        = (tmp(1:end-1)+tmp(2:end))/2/mean(abs(tmp));

fac0        = max(dtmp,ddtmp);

Dwadap      = Dw*ones(nw-1,1);
% hold on;plot(Dwadap)
facm        = 80;%interpole
facp        = 40;%decime
fac         = min(max(1/facp,fac0),facm);
% navg        = 10;
% fac         = mvg_avg(fac,navg,1,2);
Dwadap      = Dwadap./fac;
% Dwadap      = mvg_avg(Dwadap,navg,0,1);

% plot(DWN.omegac(1:end-1),Dwadap,'r')
wa=omegac;
i=1;
j=1;
Dwa=omegac;%init length
Dwa(1)=Dw;
while wa(i)<omegac(end) % i<nw % 
    %trouver l indice du Dwadap anterieur
    i=i+1;
%     j1=j;
    while omegac(j)<wa(i-1)-Dwadap(j)/2 && j<nw-1
        j=j+1;
    end
    j=min(max(j-1,1),nw);
    j1=j;
    %trouver l indice du Dwadap superieur
    j2=j;
    while omegac(j2)<wa(i-1)+Dwadap(j)/2 && j2<nw-1
        j2=j2+1;
    end
    j2=min(max(j2-1,1),nw);

    dw0     = min(Dwadap(j1:j2));
    wa(i)  = wa(i-1)+Dwa(i-1)/2+dw0/2;
    Dwa(i)  = dw0;
end
% figure;hold on;plot(omegac);plot(wa,'r')
DWN.omegac =wa(1:i);
DWN.domegac=Dwa(1:i);