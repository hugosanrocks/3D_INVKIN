function DWN=rebuildk2(DWN)
nk2=length(DWN.k2);
tmp=zeros(nk2,1);
for i=1:nk2;
    tmp(i)=1./abs(det(DWN.A_DWN(:,:,i)));
end;

k2  = DWN.k2;

DK  = k2(2)-k2(1);
% 
% figure(205);hold on
% plot(DWN.k2,tmp,'-o')

ddtmp       = abs(diff(tmp));
ddtmp       = [0;(ddtmp(1:end-1)+ddtmp(2:end))/2]/mean(ddtmp);
dtmp        = (tmp(1:end-1)+tmp(2:end))/2/mean(abs(tmp));

fac0        = max(dtmp,ddtmp);

DKadap      = DK*ones(nk2-1,1);
% hold on;plot(DKadap)
facm        = 80;%interpole
facp        = 40;%decime
fac         = min(max(1/facp,fac0),facm);
% navg        = 10;
% fac         = mvg_avg(fac,navg,1,2);
DKadap      = DKadap./fac;
% DKadap      = mvg_avg(DKadap,navg,0,1);

% plot(DWN.k2(1:end-1),DKadap,'r')
k2a=k2;
i=1;
j=1;
DKa=k2;%init length
DKa(1)=DK;
while k2a(i)<k2(end) % i<nk2 % 
    %trouver l indice du DKadap anterieur
    i=i+1;
%     j1=j;
    while k2(j)<k2a(i-1)-DKadap(j)/2 && j<nk2-1
        j=j+1;
    end
    j=min(max(j-1,1),nk2);
    j1=j;
    %trouver l indice du DKadap superieur
    j2=j;
    while k2(j2)<k2a(i-1)+DKadap(j)/2 && j2<nk2-1
        j2=j2+1;
    end
    j2=min(max(j2-1,1),nk2);

    dk0     = min(DKadap(j1:j2));
    k2a(i)  = k2a(i-1)+DKa(i-1)/2+dk0/2;
    DKa(i)  = dk0;
end
% figure;hold on;plot(k2);plot(k2a,'r')
DWN.k2 =k2a(1:i);
DWN.dk2=DKa(1:i);