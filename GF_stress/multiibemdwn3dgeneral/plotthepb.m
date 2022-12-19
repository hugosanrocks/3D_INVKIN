function plotthepb(wi,wf,para,DWN)

for j=1:length(wi)
    fj0         = linspace(wi(j),wf(j),1e4)/(2*pi);
    para        = attenuation(para,fj0);
    para.fj     = fj0;
    % para0       = para;
    DWN.omegac	= 2*pi*fj0;
    
    if  para.pol==1
        DWN     = calcul_A_DWN_SH_Ncapas_HS(para,DWN);
    else
        DWN     = calcul_A_DWN_PSV_Ncapas_HS(para,DWN);
    end
    tmp=zeros(1e4,1);
    for i=1:1e4
        tmp(i)  = (det(DWN.A_DWN(:,:,i)));
    end
    figure;hold on
    plot(2*pi*fj0,real(tmp))
    plot(2*pi*fj0,imag(tmp),'r')
%     plot(2*pi*fj0,angle(tmp),'c')
end