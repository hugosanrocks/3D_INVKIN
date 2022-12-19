function sydiff=invertT22(phi_fv,fj,coordf,para,m,xr,zr,gaussian)

%%%%%%%%%%%%%
% unpacking %
%%%%%%%%%%%%%
ninc    = para.ninc;

%posicion de las fuentes virtuales
xf      = coordf.x;
zf      = coordf.z;
dr      = coordf.dr;
phi     = coordf.phi;

%propriedad del medio
mu      = para.reg(m).Ci(6,6);
ksi     = para.reg(m).ksi;

%init el campo difractado
sydiff  = zeros(1,ninc);

%posicion de las fv que hay que tomar en cuenta
j       = 1:coordf.nbpt;
jphi    = 1:coordf.nbeq;
jjphi   = false(coordf.nbeq,1);
jjx     = coordf.indm(m).ind;
jjphi(phi(j(jjx),m)) = true(1);
jj      = jphi(jjphi);
ii      = j(jjx);

%calculo de la suma de cada una de las contribuciones
xij     = xr-xf(ii);
zij     = zr-zf(ii);
drj     = dr(ii);
rij     = sqrt(xij.^2+zij.^2);
% drdn	=(xij*vnxf+zij*vnzf)./rij;

t22      = 1i/4*ksi.*besselh(1,2,ksi.*rij);
t221     = t22.*xij./rij;
t222     = t22.*zij./rij;

%integracion gaussiana donde se requiere
jjtmp2=find((rij<para.npplo*drj).*(rij>=0));
gn22(jjtmp2)=G22_SH_r_small(coordf,xr,zr,ii(jjtmp2),ksi,gaussian,mu);

jjtmp2=find((rij<=para.npplo*drj));
t22(jj(jjtmp2)) = T22_SH_r_small(coordf,1,jr(jjtmp2),ksi,gaussian);

t22(jj) = ((m==1) - (m~=1))*t22(jj).*drj;

%cuando rij=0 %### a checar cuando hay atenuacion como ??
signo1= ((m==1) + (m~=1))*(-(coord.Xim(i,m)==1)+(coord.Xim(i,m)==2));
t22(jj(rij==0))  = t22(jj(rij==0))+signo1*0.5;
%le signe1 dépend des normales, il vaut 1 si la normale pointe vers
%l interieur ou -1 si la normale pointe vers l exterieur
%la convention de la normale est d etre toujours pointee vers le
%bas de facon a permettre une ecriture plus facile de la continuite
%des tractions normales
%en ce point, il y a discontinuite des tractions ti-ti'=phi



for iinc=1:ninc
    t22diff( 1,iinc)=sum(t22 .*drj.*phi_fv(jj,iinc).');
    t221diff(1,iinc)=sum(t221.*drj.*phi_fv(jj,iinc).');
    t222diff(1,iinc)=sum(t222.*drj.*phi_fv(jj,iinc).');
end