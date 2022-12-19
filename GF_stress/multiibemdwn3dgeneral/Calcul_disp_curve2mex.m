% function Calcul_disp_curve
%se ingresa a mano los datos de entrada y los acomoda para las funciones
%que calculan las curvas de dispersion


%parametros que se pueden cambiar
para.pol        = 2;
para.nsubmed    = 4;
para.fmax       = 3;


%parametros que no se pueden cambiar
para.sub        = struct('h',0,'rho',0,'bet',0,'alpha',0);
para.nmed       = 1;
para.geo     	= 3;
pol             = para.pol;
nc              = para.nsubmed;

h               = zeros(para.nsubmed,1);
vS              = zeros(para.nsubmed,1);
vP              = zeros(para.nsubmed,1);
rho             = zeros(para.nsubmed,1);

h(1)    = 1.1;
vS(1)   = 1;
vP(1)   = 2;
rho(1)  = 1;

h(2)    = 4;
vS(2)   = 1.5;
vP(2)   = 3;
rho(2)  = 1.5;

h(3)    = 30;
vS(3)   = 3;
vP(3)   = 6;
rho(3)  = 2;

h(4)    = 2;
vS(4)   = 5;
vP(4)   = 10;
rho(4)  = 2.2;

h(5)    = 0;
vS(5)   = 6;
vP(5)   = 12;
rho(5)  = 2.5;

% tic
% if para.pol==2
%     var_guess=zeros(4*nlay-1,1);
%     n0=4;
% else
%     var_guess=zeros(3*nlay-1,1);
%     n0=3;
% end
% 
% for im  = 1:nlay
%     var_guess(1+(im-1)*n0)  = para.reg(1).sub(im).rho;
%     var_guess(2+(im-1)*n0)  = para.reg(1).sub(im).bet;
%     if para.pol==2
%         var_guess(3+(im-1)*n0)  = para.reg(1).sub(im).alpha;
%         if im~=nlay
%             var_guess(4+(im-1)*n0)  = para.reg(1).sub(im).h;
%         end
%     else
%         if im~=nlay
%             var_guess(3+(im-1)*n0)  = para.reg(1).sub(im).h;
%         end
%     end
% end
for i=1:para.nsubmed
    para.sub(i).h      = h(i);
    para.sub(i).rho    = rho(i);
    para.sub(i).bet    = vS(i);
    para.sub(i).alpha  = vP(i);
    para.sub(i).tipoatt= vP(i);
end

%pour charger ds IBEM
for i=1:para.nsubmed
    para.sub(i).qd=1000;
    para.sub(i).tipoatts=3;
    para.sub(i).lambda=3;
    para.sub(i).mu=3;
end
para.reg(1).sub=para.sub;
% nmodem          = 2;
wi           = linspace(0,para.fmax,100)*2*pi;
% tic
% [vg,f1,ikmax]= dispersion_curve_wfix_Haskel_4inv_adapt(para,wi);
% toc

checktime_disp(para,wi)

tic
% for i=1:10
[vp,f1,ikmax]=dispersion_curve_VP_wfix_Haskel_4inv(para,wi);
% end
toc

tic
% for i=1:10
nmode=length(ikmax);

[vp,f1,ikmax]=dispersion_curve_wfix_Wathelet(para,wi,nmode);
% end
toc

tic
% for i=1:10
[vp,f1,ikmax]=dispersion_curve_VP_wfix_Haskel_4inv_test(para,wi);
% end
toc




% 
% tic
% % for i=1:10
% [vg,f1,ikmax]=dispersion_curve_wfix_Haskel_4inv(para,wi);
% % end
% toc
% figure(207);hold on;
% nmode=length(ikmax);
% for j=1:nmode
%     plot(f1(j,1:ikmax(j)),vg(j,1:ikmax(j)),'k')
% end
% 
% nc=para.nsubmed;
% tic
% % for i=1:10
% [vg,f1,ikmax]=dispersion_curve_wfix_Haskel_4inv_adapt2mex_mex(wi,nc,pol,vS,vP,rho,h);
% % end
% toc
% figure(207);hold on;
% for j=1:nmode
%     plot(f1(j,1:ikmax(j)),vg(j,1:ikmax(j)),'r')
% end
% 
% 
% tic
% % for i=1:10
% [vg,f1,ikmax]=dispersion_curve_wfix_Wathelet(para,wi,nmode);
% % end
% toc
% figure(207);hold on;
% for j=1:nmode
%     plot(f1(j,1:ikmax(j)),vg(j,1:ikmax(j)),'c')
% end