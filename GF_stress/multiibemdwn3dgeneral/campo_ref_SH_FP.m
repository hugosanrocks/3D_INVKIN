function [uy,tn,s]=campo_ref_SH_FP(xs,zs,coord,ks,mu,salu,salt,para)
% funcion que da los campos de desplazamientos y las tractiones
% para una fuente puntual SH en un semi-espacio
% xs,zs coordenadades del origen de la onda (fase nula)
% #r relativo a los receptores

n       = coord.nbpt;
uy     	= zeros(1,n);
s      	= 0;
tn     	= zeros(1,n);

xr  	= coord.x;
zr  	= coord.z;

x       = xr-xs;
z       = zr-zs;
r       = sqrt(x.^2+z.^2);

uy(salu)      = G22_SH(ks,r(salu),mu);

if isfield(coord,'vnx')
    vnxr  	= coord.vnx;
    vnzr  	= coord.vnz;
    drdn    =(x.*vnxr+z.*vnzr)./r;
    
    tn(salt)= 1i/4*ks.*besselh(1,2,ks*r(salt)).*drdn(salt);%+besselh(1,2,ks*r_im).*drdn_im);
    
    ipb=find(r<1*para.npplo*coord.dr);
    if ~isempty(ipb)
        gaussian   = para.gaussian;
        j=1:n;
        j1=j(ipb);
        if ~isempty(j1(salu(ipb)))
            uy(j1(salu(ipb))) = G22_SH_r_small_FP(coord,xs,zs,j1(salu(ipb)),ks,gaussian,mu);
        end
        if ~isempty(j1(salt(ipb)))
            tn(j1(salt(ipb))) = T22_SH_r_small_FP(coord,xs,zs,j1(salt(ipb)),ks,gaussian);
        end
    end
else
    tn=[];
    ipb=find(r==0);
     if ~isempty(ipb(salu(ipb)))
        uy(ipb)	= 1/(4*1i*mu);
     end
     
     nss     = (para.sortie.sxy + para.sortie.syz);
     if nss>0
         s = zeros(nss,coord.nbpt);
         if para.sortie.sxy==1
             s(1,salt)=1i/4*ks.*besselh(1,2,ks*r(salt)).*x(salt)./r(salt);
         end
         if para.sortie.syz==1
             s(nss,salt)=1i/4*ks.*besselh(1,2,ks*r(salt)).*z(salt)./r(salt);
         end
     end

end