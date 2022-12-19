function detA=zerodeAw(para,DWN,w)

DWN.omegac  = w;
fj      = DWN.omegac/(2*pi);
para = attenuation(para,fj);

if  para.pol==1
    DWN     = calcul_A_DWN_SH_Ncapas_HS(para,DWN);
else
    DWN     = calcul_A_DWN_PSV_Ncapas_HS(para,DWN);
end

nw  = length(DWN.omegac);
detA= zeros(nw,1);

for i=1:nw;
    detA(i)=real(det(DWN.A_DWN(:,:,i)));
end
