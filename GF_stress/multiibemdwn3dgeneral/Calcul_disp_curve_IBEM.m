% function Calcul_disp_curve_IBEM(para)
%se occupa los datos ya ingresados en la interfase de Code_IBEM para
%calcular las curvas de dispersion


% wi           = linspace(0,para.fmax,100)*2*pi;
% tic
% [vg,f1,ikmax]= dispersion_curve_wfix_Haskel_4inv_adapt(para,wi);
% toc


% tic
% % for i=1:10
% [vp,f1,ikmax]=dispersion_curve_VP_wfix_Haskel_4inv(para,wi);
% % end
% toc
tic
% for i=1:10
%[vp,f1,ikmax]=dispersion_curve_VP_wfix_Haskel_4inv_test(para,wi);
[vg,vp,f1,ikmax,kx]=dispersion_curve_VP_wfix_Haskel_4inv_u2d3(para,wi);

% end
toc
tic
nmode=size(vg,1);%100;
[vp,f1,ikmax]=dispersion_curve_wfix_Wathelet(para,wi,nmode);
toc




% tic
% [~,~,k2,w0,ikmax]	= dispersion_curve_kfix_MG_fast2(para);
% toc