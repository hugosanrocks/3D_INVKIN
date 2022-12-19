function [sign0,sign1]=checkchgmtsign(wi,wf,para,DWN)

DWN.omegac  = wi;
n           = length(DWN.omegac);
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
deti        = zeros(n,1);
for i=1:n
    deti(i)        = (det(det0.A_DWN(:,:,i)));
end
detir       = real(deti);
detii       = imag(deti);


DWN.omegac  = wf;
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

detf        = zeros(n,1);
for i=1:n
    detf(i)        = (det(det0.A_DWN(:,:,i)));
end
detfr       = real(detf);
detfi       = imag(detf);

sign0=sign(detir.*detfr);
sign1=sign(detii.*detfi);