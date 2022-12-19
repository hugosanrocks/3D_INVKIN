function l=malla_contorno_fast2(cont1,para,fj)
% funcion que permite posicionar los puntos approximadamente separados de
% dr0, los puntos se situan en medio de cada segmento

fj      = abs(fj);
npplo   = para.npplo;

n       = length(cont1.vec.xc);
l       = zeros(n,1);

%inicializacion
%calculo de la velocidad minima de cada lado del contorno en el punto
%corriente
i       = 1;
m       = cont1.m;          %medio interior
m1      = cont1.mv(i);     %medio exterior al punto corriente
if m1>0
    m1      = para.subm1(m1);
end
indm    = [m;m1];
indm(indm==0)   = [];
indm    = squeeze(indm);
v       = 0*indm; %init
for k=1:length(indm)
    if para.reg(indm(k)).rho==0
        v(k) = 0;
    else
        v(k)       = para.reg(indm(k)).bet;
        if v(k)==0
            v(k)       = para.reg(indm(k)).alpha;
        end
    end
    %         if indm(k)>1
    %             v(k)       =  v(k)/30;
    %         end
end
v(v==0) = [];
v       = squeeze(v);
if isempty(v)
    l=[];
    return
end
cmin    = min(v);
lambda  = cmin/fj;
if para.dim==1
    nmin=32;%32
else
    nmin=8;
end
% nmin=1;
nbpt    = max(nmin,round(cont1.vec.r(end)/lambda*npplo));
% nbpt    = 4;%07/07/15 debug only ###
dr      = cont1.vec.r(end)/nbpt;
dr2     = dr/2;                 %primer punto
rc      = 0;
rd      = rc + dr2;             %la posicion deseada es la posicion corriente + el dr
[~,imin]= min(abs(rd-cont1.vec.r(1:end)));
l(i)    = imin;
i       = i+1;
j       = imin;
% dr0     = dr;
% test=0;
% imin0=imin;
while rd<cont1.vec.r(end)
    rc      = cont1.vec.r(imin);
    %     if cont1.vec.zc(imin)~=0 && cont1.vec.zc(imin)~=0.5 && cont1.vec.zc(imin)~=-0.2 || test==1
    %         dr=dr0/10;
    %         test=0;
    %     else
    %         dr=dr0;
    %     end
    rd      = rc + dr;          %la posicion deseada es la posicion corriente + el dr
    [~,imin]= min(abs(rd-cont1.vec.r(j:end)));
    imin    = j-1+imin;
    %     if cont1.vec.zc(imin)~=0 && cont1.vec.zc(imin)~=0.5 && cont1.vec.zc(imin)~=-0.2  && dr==dr0
    %         test=1;
    %         imin=imin0;
    %     else
    l(i)    = imin;
    i       = i+1;
    j       = imin;
%     imin0   = imin;
    %     end
end
l(i-1:end)  = [];
l           = squeeze(l);