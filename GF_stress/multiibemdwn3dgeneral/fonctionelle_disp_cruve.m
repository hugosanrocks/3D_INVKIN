function err=fonctionelle_disp_cruve(f_in,vg_in,para,var2opt)

nlay    = para.nsubmed;

if para.pol==2
    n0=4;
else
    n0=3;
end

for im  = 1:nlay
    para.reg(1).sub(im).rho     = var2opt(1+(im-1)*n0);
    para.reg(1).sub(im).bet     = var2opt(2+(im-1)*n0);
    if para.pol==2
        para.reg(1).sub(im).alpha   = var2opt(3+(im-1)*n0);
        if im~=nlay
            para.reg(1).sub(im).h       = var2opt(4+(im-1)*n0);
        end
    else
        if im~=nlay
            para.reg(1).sub(im).h       = var2opt(3+(im-1)*n0);
        end
    end
end
para	= Vp2Cij(para);

% [vg,f1,k2,w0,ikmax]=dispersion_curve_kfix_MG_fast2(para);
[f_in_u,~,ifu]  = unique(f_in);
fpb=f_in_u;
fpb(end+1)=2*fpb(end)-fpb(end-1);
[vg,f1,ikmax]   = dispersion_curve_wfix_Haskel_4inv(para,2*pi*fpb);
% [vg,f1,ikmax]   = dispersion_curve_wfix_Haskel_4inv_adapt(para,2*pi*fpb);
% [vg,f1,ikmax]=dispersion_curve_wfix_Wathelet(para,2*pi*fpb,3);

nmode           = length(ikmax);


vg_it           = zeros(nmode,length(f_in));
for imode = 1:nmode
    if ikmax(imode)>1
        vgtmp               = interp1(f1(imode,1:ikmax(imode)),vg(imode,1:ikmax(imode)),f_in_u);
        vgtmp(isnan(vgtmp)) = 1e100;
    else
        vgtmp               = ones(1,length(f_in))*1e100;
    end
    vg_it(imode,:)   	= vgtmp(ifu);
end
err=0;
n0=length(vg_in);
n1=n0;
dvg = diff(vg_in);
df  = diff(f_in);
df(df>.2)=100;
dvg =abs(dvg./df);
dvg(dvg>1)=1;
dvg(dvg<.01)=.01;
dvg(end+1)=dvg(end);
for i=1:n0
    [min0,indv]=min(abs(vg_it(:,i)-vg_in(i)));
    if min0==1e100
        err=err+abs(vg_in(i))^2;
    else
%         err=err+abs(vg_it(indv,i)-vg_in(i))^2;
        err=err+dvg(i)*abs(vg_it(indv,i)-vg_in(i))^2;
    end
    
    
end
err=err/n1;