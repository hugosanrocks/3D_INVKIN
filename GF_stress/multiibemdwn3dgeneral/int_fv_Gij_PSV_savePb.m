function uydiff=int_fv_Gij_PSV(phi_fv,coordf,para,m,xr,zr,gaussian,kpi,ksi,C)

%%%%%%%%%%%%%
% unpacking %
%%%%%%%%%%%%%
ninc    = para.ninc;
nrecm   = length(xr);

%init el campo difractado
uydiff  = zeros(2,ninc,nrecm);

%posicion de las fuentes virtuales
xf      = coordf.x;
zf      = coordf.z;
dr      = coordf.dr;
phi     = coordf.phi;
nbeq    = coordf.nbeq;

%posicion de las fv que hay que tomar en cuenta
j       = 1:coordf.nbpt;
jphi    = 1:coordf.nbeq;
jjphi   = false(coordf.nbeq,1);
jjx     = coordf.indm(m).ind;
jjphi(phi(j(jjx),m)) = true(1);
jj      = jphi(jjphi);
ii      = j(jjx);
njj     = length(jj);

%calculo de la suma de cada una de las contribuciones
xij=zeros(nrecm,njj);
zij=zeros(nrecm,njj);
drj=zeros(nrecm,njj);
for i=1:nrecm
    xij(i,:)= xr(i)-xf(ii);
    zij(i,:)= zr(i)-zf(ii);
    drj(i,:)= dr(ii);
end
xij=reshape(xij,1,njj*nrecm);
zij=reshape(zij,1,njj*nrecm);
drj=reshape(drj,1,njj*nrecm);

% xij     = xr-xf(ii);
% zij     = zr-zf(ii);
% drj     = dr(ii);
rij     = sqrt(xij.^2+zij.^2);
g       = zeros(2,njj*nrecm);
g(1,:)  = xij./rij;
g(2,:)  = zij./rij;

Gij     = Gij_PSV(ksi,kpi,rij,g,C,njj*nrecm);

Gij     = reshape(Gij,2,2,njj,nrecm);
drj     = reshape(drj,1,njj,nrecm);
rij     = reshape(rij,1,njj,nrecm);

for i=1:nrecm
    %integracion gaussiana donde se requiere
    jjtmp2  = find(rij(:,i)<5*para.npplo*drj(:,i));
    Gij(:,:,jjtmp2)=Gij_PSV_r_small(coordf,xr(i),zr(i),ii(jjtmp2),ksi,kpi,gaussian,C);
    for iinc=1:ninc
        uydiff(1,iinc,i)=sum(squeeze(Gij(1,1,:,i)).*drj.'.*phi_fv(jj,iinc)+squeeze(Gij(1,2,:,i)).*drj.'.*phi_fv(jj+nbeq,iinc));
        uydiff(2,iinc,i)=sum(squeeze(Gij(2,1,:,i)).*drj.'.*phi_fv(jj,iinc)+squeeze(Gij(2,2,:,i)).*drj.'.*phi_fv(jj+nbeq,iinc));
    end
end