function [sign0,det1,para,DWN]=checkchgmtsign_phase_vec(wi,wf,para,DWN,n)


DWN.omegac  = linspace(wi,wf,n);
fj          = DWN.omegac/(2*pi);
para        = attenuation(para,fj);
if  para.pol==1
    det0        = calcul_A_DWN_SH_Ncapas_HS(para,DWN);
else
    det0        = calcul_A_DWN_PSV_Ncapas_HS(para,DWN);
end
det1        = zeros(n,1);
for i=1:n
    det1(i)        = (det(det0.A_DWN(:,:,i)));
end
sign0=sign(angle(det1(1)).*angle(det1(end)));
