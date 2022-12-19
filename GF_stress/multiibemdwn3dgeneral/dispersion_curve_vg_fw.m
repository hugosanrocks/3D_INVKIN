function [kxw,vgw]=dispersion_curve_vg_fw(para,w)

% a partir de las curbas de dispersion, busca precisamente el k2
% correspndiendo al w dado por cada modo

ikmax   = para.mode.ikmax;
w1      = para.mode.f1*2*pi;
vg      = para.mode.vg;
kx      = para.mode.k2;

nmode   = size(w1,1);
vgw     = zeros(nmode,1);
kxw     = zeros(nmode,1);

for imode=1:nmode
    if w>=w1(imode,1)
        %approx de la solution par interpolation
        vgw(imode)  = interp1(w1(imode,1:ikmax(imode)),vg(imode,1:ikmax(imode)),real(w),'pchip','extrap');
        kxw(imode)  = interp1(w1(imode,1:ikmax(imode)),kx(imode,1:ikmax(imode)),real(w),'pchip','extrap');
    else
        vgw(imode:nmode)=[];
        kxw(imode:nmode)=[];
        break
    end
end