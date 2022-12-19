function detA=zerodeAwIR(para,DWN,w,part)
%part==1 : partie reelle
%part~=1 : partie imag

DWN.omegac  = w;
fj      = DWN.omegac/(2*pi);
para = attenuation(para,fj);

if  para.pol==1
    DWN     = calcul_A_DWN_SH_Ncapas_HS(para,DWN);
else
    DWN     = calcul_A_DWN_PSV_Ncapas_HS(para,DWN);
end

nw=length(DWN.omegac);
detA=zeros(nw,1);
if part==1
    for i=1:nw;
        detA(i)=real(det(DWN.A_DWN(:,:,i)));%-imag(det(DWN.A_DWN(:,:,i)));
    end
else
    for i=1:nw;
        detA(i)=imag(det(DWN.A_DWN(:,:,i)));%-imag(det(DWN.A_DWN(:,:,i)));
    end
end