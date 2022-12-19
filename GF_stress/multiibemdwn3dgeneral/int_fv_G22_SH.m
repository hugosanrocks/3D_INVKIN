function [uydiff,sdiff]=int_fv_G22_SH(phi_fv,coordf,para,m,xr,zr,gaussian)
%integracion de la contribucion de todas las fuentes virtuales

%%%%%%%%%%%%%
% unpacking %
%%%%%%%%%%%%%
ninc    = para.ninc;
sal     = para.sortie;

%posicion de las fuentes virtuales
xf      = coordf.x;
zf      = coordf.z;
dr      = coordf.dr;
phi     = coordf.phi;

%init el campo difractado
uydiff  = zeros(1,ninc);
nss     = (sal.sxy + sal.syz);
sdiff   = zeros(nss,ninc);

%propriedad del medio
if m==1 && para.geo(1)==3
    mu      = para.reg(1).sub(1).Ci(6,6);
    ksi     = para.reg(1).sub(1).ksi;
else
    mmat    = para.subm1(m);
    mu      = para.reg(mmat).Ci(6,6);
    ksi     = para.reg(mmat).ksi;
end

%posicion de las fv que hay que tomar en cuenta
j       = 1:coordf.nbpt;
jphi    = 1:coordf.nbeq;

%indice (logic y natural) de los puntos de colocacion perteneciendo a m
jjx     = coordf.indm(m).ind;
ii      = j(jjx);

%indice (logic y natural) de los phi que hay que contemplar (columnas
%de la matriz)
jjphi   = false(coordf.nbeq,1);
jjphi(phi(j(jjx),m)) = true(1);
jj      = jphi(jjphi);
        
%calculo de la suma de cada una de las contribuciones
xij = xr-xf(ii);
zij = zr-zf(ii);
drj = dr(ii);
rij = sqrt(xij.^2+zij.^2);

gn22        = G22_SH(ksi,rij,mu);
jjtmp2      = find((rij<1*para.npplo*drj));%.*(rij>=0));
gn22(jjtmp2)= G22_SH_r_small(coordf,xr,zr,ii(jjtmp2),ksi,gaussian,mu);
if sal.sxy || sal.syz
    tmp    = 1i/4.*ksi.* besselh(1,2,ksi*rij);
    Sxy    = tmp.*xij./rij;
    Syz    = tmp.*zij./rij;
    if ~isempty(jjtmp2)
        coordf2.x   = coordf.x(ii(jjtmp2)).';
        coordf2.z   = coordf.z(ii(jjtmp2)).';
        coordf2.vnx = coordf.vnx(ii(jjtmp2)).';
        coordf2.vnz = coordf.vnz(ii(jjtmp2)).';
        coordf2.dr  = coordf.dr(ii(jjtmp2)).';
        [Sxy(jjtmp2),Syz(jjtmp2)]=Sxzy_SH_r_small(coordf2,xr,zr,ksi,gaussian);
        jjtmp3      = find(rij<=para.npplo*drj/2);
        signo1= (coordf.Xim(ii(jjtmp3),m)==1) - (coordf.Xim(ii(jjtmp3),m)==2);
        
        vnxpb=coordf.vnx(ii(jjtmp3)).';
        vnzpb=coordf.vnz(ii(jjtmp3)).';
        
        Sxy(jjtmp3) =0*Sxy(jjtmp3)+(signo1*0.5.*vnxpb./sqrt(vnxpb.^2+vnzpb.^2)).'/drj(jjtmp3);
        Syz(jjtmp3) =0*Syz(jjtmp3)+(signo1*0.5.*vnzpb./sqrt(vnxpb.^2+vnzpb.^2)).'/drj(jjtmp3);
    end
end

if sal.USh==1 || sal.UIh==1
    %angle entre le nouveau repere lie a l element et le repere IBEM
    vnx   	= coordf.vnx(ii);
    vnz 	= coordf.vnz(ii);
    th1     = atan2(vnz,vnx)+pi/2;
    th0     = atan2(zij,xij);
    xij1    = rij.*cos(th0-th1);
    zij1    = rij.*sin(th0-th1);
    G22r	= G22_S_L(ksi,xij1,zij1,rij,mu,drj);
    G22i    = gn22-G22r;
end


if m==1 && para.geo(1)==3 && para.nsubmed==1
    %fuente imagen
    zij             = zr+zf(ii);
    rij             = sqrt(xij.^2+zij.^2);
    gn221           = G22_SH(ksi,rij,mu);
    coordf.z(ii)	=-coordf.z(ii);
    coordf.vnz(ii)	=-coordf.vnz(ii);
    jjtmp2          = find((rij<1*para.npplo*drj).*(rij>=0));
    gn221(jjtmp2)   = G22_SH_r_small(coordf,xr,zr,ii(jjtmp2),ksi,gaussian,mu);
    gn22            = gn22+gn221;
    if sal.USh==1 || sal.UIh==1
        %angle entre le nouveau repere lie a l element et le repere IBEM
        vnx   	= coordf.vnx(ii);
        vnz 	= coordf.vnz(ii);
        th1     = atan2(vnz,vnx)-pi/2;
        th0     = atan2(zij,xij);
        xij1    = rij.*cos(th0-th1);
        zij1    = rij.*sin(th0-th1);
        %         [G22r1,G22i1]	= G22_S_L(ksi,xij1,zij1,rij,mu);
        %         G22r=G22r+G22r1;
        %         G22i=G22i+G22i1;
        G22r1   = G22_S_L(ksi,xij1,zij1,rij,mu);
        G22r    = G22r+G22r1;
        G22i    = G22i+gn221-G22r1;
    end
end
%integracion de las contribuciones
for iinc=1:ninc
    k=1;
    if sal.Ut==1
        uydiff(k,iinc)=sum(gn22.*drj.*phi_fv(jj,iinc).');
        k=k+1;
    end
    if sal.USh==1
        uydiff(k,iinc)=sum(G22r.*drj.*phi_fv(jj,iinc).');
        k=k+1;
    end
    if sal.UIh==1
        uydiff(k,iinc)=sum(G22i.*drj.*phi_fv(jj,iinc).');
    end
    k=1;
    if sal.sxy==1
        sdiff(k,iinc)=sum(Sxy.*drj.*phi_fv(jj,iinc).');
        k=k+1;
    end
    if sal.syz==1
        sdiff(k,iinc)=sum(Syz.*drj.*phi_fv(jj,iinc).');
    end
end