function [sign0,sign1,det1,para,DWN]=checkchgmtsign_vec(wi,wf,para,DWN,n)


DWN.omegac  = linspace(wi,wf,n);
for ic=1:para.nsubmed
    para.reg(1).sub(ic).ksi  = DWN.omegac/para.reg(1).sub(ic).bet;
end
if  para.pol==1
    det0        = calcul_A_DWN_SH_Ncapas_HS(para,DWN);
else
    for ic=1:para.nsubmed
        para.reg(1).sub(ic).kpi  = DWN.omegac/para.reg(1).sub(ic).alpha;
    end
    det0        = calcul_A_DWN_PSV_Ncapas_HS(para,DWN);
end
det1        = zeros(n,1);
for i=1:n
    det1(i)        = (det(det0.A_DWN(:,:,i)));
end
detr       = real(det1);
deti       = imag(det1);

sign0=sign(detr(1).*detr(end));
sign1=sign(deti(1).*deti(end));