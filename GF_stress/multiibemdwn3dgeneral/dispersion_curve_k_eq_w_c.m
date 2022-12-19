% function [k2,w0,ikmax]=dispersion_curve_k_eq_w_c(para)

%function qui calcule sur des diagonales similaires à celles des asymptotes
%peu efficace car l avantage des asymptotes est surtout que les zeros sont
%equi-espacés


tic

para.DWNomei= 0;
para        = Vp2Cij(para);
wmax        = 2*pi*para.fmax;
ncapas      = para.nsubmed;
h           = zeros(ncapas,1);
vT          = zeros(ncapas,1);

for im  = 1:ncapas-1
    h(im)   = para.reg(1).sub(im).h;
    C       = para.reg(1).sub(im).C;
    rho     = para.reg(1).sub(im).rho;
    vT(im) 	= sqrt(C(6,6)/rho);
end
im      = ncapas;
C       = para.reg(1).sub(im).C;
rho     = para.reg(1).sub(im).rho;
vT(im) 	= sqrt(C(6,6)/rho);
vmax    =max(vT);
vTtmp   =vT;
vTtmp(vT==vmax)=[];
vmax2   =max(vTtmp);
vmin    = min(vT);

j           = 1;
v           = vmax;
[k20,w00]   = dispersion_curve_k_critik_MG(para,wmax,v);
nmode       = length(k20);
npt         = 10;
k2          = zeros(npt,nmode);
w0          = zeros(npt,nmode);  
k2(1,:)     = k20;
w0(1,:)     = w00;

vscan       = linspace(vmax,vmin,npt);
dv          = vscan(1)-vscan(2);
vscan       = linspace(vmax-dv,vmin,npt-1);

for v=vscan
    if max(v==vT)==1
        v=v*(1-1e-8);
    end
    [k20,w00]   = dispersion_curve_k_critik_MG(para,wmax,v);
    nmv         = length(k20);
    j           = j+1;
    k2(j,1:nmv-1)	= k20(2:nmv);
    w0(j,1:nmv-1)	= w00(2:nmv);
end
figure(205);hold on;
i=1;
j=find(k2(:,i)==0,2,'first');
j=j-1;
plot(k2(1:j,i),w0(1:j,i),'c')
for i=1:nmode
    j=find(k2(:,i)==0,1);
    j=j-1;
    plot(k2(1:j,i),w0(1:j,i),'c')
end
