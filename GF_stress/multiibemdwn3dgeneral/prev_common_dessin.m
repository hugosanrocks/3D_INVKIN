nf      = para.nf;
df      = para.fmax/(para.nf/2);     %paso en frecuencia
Fq      = (0:nf/2)*df;
nfN     = nf/2+1; %Nyquist
% zerospad= para.zeropad;

% tps     = 0:(1/(df*2*(nfN+zerospad))*(2*(nfN+zerospad)/(2*(nfN+zerospad)-2))):1/df;
dt = 1/(df*para.zeropad);    %disp(['dt = ',num2str(dt)])
tps = (0:para.zeropad-1)*dt; %disp(['tmx= ',num2str(tps(end))])

para.rec.nrec   = para.rec.nrecx*para.rec.nrecy*para.rec.nrecz;
nrec    = para.rec.nrec;
para    = normalizacion(para);
sal     = para.sortie;

%valores posibles
%                 Ux Uy Uz Ut UPh USh UIh UPt USt sxx syy szz sxy sxz syz
if para.dim==1
    if para.pol==1 %SH
        sal0    = [0  1  0  1  0   1   1   0   0   0   0   0   1   0   1];
    else%PSV
        sal0    = [1  0  1  1  1   1   1   1   1   1   0   1   0   1   0];
    end
else
    ;   sal0    = [1  1  1  1  0   0   0   0   0   1   1   1   1   1   1];
end
fn0     = fieldnames(para.sortie);
for ifn0= 1:length(fn0)
    tmp = getfield(para.sortie,fn0{ifn0});
    sal0(ifn0)  = sal0(ifn0)*tmp;
end
ns0       = 0;
clear fieldV
%composante deplacement
for ifn0= 4:9
    if sal0(ifn0)
        ns0=ns0+1;
        fieldV(ns0).name=fn0{ifn0};
        fieldV(ns0).nc=sum(sal0(1:3));
        fieldV(ns0).namec(1:fieldV(ns0).nc)=fn0(logical(sal0(1:3)));
    end
end
%composante stress
if length(fn0)>9 %patch compatibilite vieux resultats
    if max(sal0(10:15))
        ns0=ns0+1;
        fieldV(ns0).name='S';
        fieldV(ns0).nc=sum(sal0(10:15));
        fieldV(ns0).namec(1:fieldV(ns0).nc)=fn0(logical([zeros(1,9),sal0(10:15)]));
    end
end
%     ns0=ns0+1;



% ns0     = double( ...
%     ((sal.Ut+sal.UPh+sal.USh+sal.UIh+sal.UPt+sal.USt)*(para.pol==2) + ...
%     (sal.Ut+sal.USh+sal.UIh)*(para.pol==1)                                  )*(para.dim==1)+...
%     (para.dim>1));
% nsxyz   = double(((sal.Ux + sal.Uz)*(para.pol==2)+(para.pol==1))*(para.dim==1)+(sal.Ux + sal.Uy + sal.Uz)*(para.dim>1));
% ns      = nsxyz*ns0;
%
%
% nssxyz  = double(((sal.sxx+sal.szz+sal.sxz)*(para.pol==2) +  (sal.sxy+sal.syz)*(para.pol==1))*(para.dim==1)+...
%     + (sal.sxx+sal.sxy+sal.sxz+sal.syy+sal.syz+sal.szz)*(para.dim>1));
% ns0     = ns0+(nssxyz>0);
% nss     = nssxyz;


% if isfield(para.rec,'xr')
%     xr      = para.rec.xr;
%     yr      = para.rec.yr;
%     zr      = para.rec.zr;
% else

if exist('cont1','var')
    para.cont1=cont1;
    tmp     = pos_rec(para);
    xr      = tmp.rec.xr;
    yr      = tmp.rec.yr;
    zr      = tmp.rec.zr;
else
    tmp     = pos_rec(para);
    xr      = tmp.rec.xr;
    yr      = tmp.rec.yr;
    zr      = tmp.rec.zr;
end