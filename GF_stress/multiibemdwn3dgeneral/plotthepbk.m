function plotthepbk(ki,kf,para,DWN)


DWN.k2  = linspace(ki,kf,1e4);

if  para.pol==1
    DWN     = calcul_A_DWN_SH_Ncapas_HS(para,DWN);
end
tmp=zeros(para.nf+1,1);
for i=1:1e4
    tmp(i)  = (det(DWN.A_DWN(:,:,i)));
end
figure;hold on
plot(DWN.k2,real(tmp))
plot(DWN.k2,imag(tmp),'r')
