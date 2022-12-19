% function Calcul_disp_curve2mex2(para,wi)
%toma los parametros del IBEM y calcula las curvas de dispersion
para.sub=para.reg(1).sub;

h   = zeros(para.nsubmed,1);
vS  = zeros(para.nsubmed,1);
vP  = zeros(para.nsubmed,1);
rho = zeros(para.nsubmed,1);
nc  = para.nsubmed;
pol = para.pol;
for i=1:para.nsubmed
    h(i)    = para.sub(i).h;
    rho(i)  = para.sub(i).rho;
    vS(i)   = para.sub(i).bet;
    vP(i)   = para.sub(i).alpha;
end

wi=linspace(0.01,para.fmax*2*pi,200);

tic
[vg,f1,ikmax]=dispersion_curve_wfix_Haskel_4inv(para,wi);
%[vg,f1,ikmax]=dispersion_curve_wfix_Haskel_4inv_adapt2mex_mex(wi,nc,pol,vS,vP,rho,h);
% [vg,f1,ikmax]=dispersion_curve_wfix_Haskel_4inv_adapt2mex(wi,nc,pol,vS,vP,rho,h);
toc
figure(207);hold on;
nmode=length(ikmax);
for j=1:nmode
    plot(f1(j,1:ikmax(j)),vg(j,1:ikmax(j)),'r')
end