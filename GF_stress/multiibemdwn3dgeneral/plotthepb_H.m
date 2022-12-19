function plotthepb_H(wi,wf,para,DWN)

for j=1:length(wi)
    DWN.omegac	= linspace(wi(j),wf(j),1e4);
    
    if  para.pol==1
        tmp=mode_Love(para,DWN);
    else
        tmp=mode_Rayleigh_2(para,DWN);
    end
    tmp(tmp==inf)=0;
    %     figure(203);hold on;plot(DWN0.omegac,real(tmp),'k');plot(DWN0.omegac,imag(tmp),'r');...
    %         plot(DWN0.omegac,abs(tmp),'');
    %     figure(203);hold on;plot(DWN0.omegac,angle(tmp)-pi/2,'c')
    
    figure;hold on
    plot(DWN.omegac,real(tmp))
    plot(DWN.omegac,imag(tmp),'r')
    plot(DWN.omegac,angle(tmp)-pi/2,'c')
end