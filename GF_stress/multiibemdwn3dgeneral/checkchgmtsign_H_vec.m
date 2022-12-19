function [signp,sign0,sign1,det1,para,DWN]=checkchgmtsign_H_vec(wi,wf,para,DWN,n)
%Haskell

DWN.omegac  = linspace(wi,wf,n);
if  para.pol==1
    det1=mode_Love(para,DWN);
else
    det1=mode_Rayleigh_2(para,DWN);
end

det1(det1==inf)=0;
tmp1    = angle(det1)-pi/2;
signp   = sign(tmp1(1).*tmp1(end));
sign0   = sign(real(det1(1)).*real(det1(end)));
sign1   = sign(imag(det1(1)).*imag(det1(end)));
