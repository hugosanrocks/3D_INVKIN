function [u,t,S]=campo_ref_FP3D(xf,yf,zf,coord,kp,ks,C,fij,salu,salt,para)

n       = coord.nbpt;
if para.dim == 3
gaussian   = para.gaussex;
else
gaussian   = para.gaussian;
end
u       = zeros(3,n);
t       = zeros(3,n);
Gij     = zeros(3,3,n);
Tij     = zeros(3,3,n);

xr  	= coord.x;
yr  	= coord.y;
zr  	= coord.z;

x       = xr-xf;
y       = yr-yf;
z       = zr-zf;
r       = sqrt(x.^2+y.^2+z.^2);
    
g       = zeros(3,n);
g(1,:)  = x./r;
g(2,:)  = y./r;
g(3,:)  = z./r;

Gij(:,:,salu)	= Gij_3D(ks,kp,r(salu),g(:,salu),C,sum(salu));

if isfield(coord,'vnx')
    %calculo de la matriz IBEM, en los puntos de colocacion
    vn(3,:) = coord.vnz;
    vn(2,:) = coord.vny;
    vn(1,:) = coord.vnx;
    
    Tij(:,:,salt)	= Tij_3D(ks,kp,r(salt),g(:,salt),C,vn(:,salt));
    %###
    ipb=find(r<para.npplo*coord.drxz);
    %integracion gaussiana donde se requiere
    
    if ~isempty(ipb) && para.dim < 4 % AQUI ESTARÍA BIEN AGREGAR Gij_3D_r_small_FP_GEN
%         [~,ipb]=min(r);
        j=1:n;
        j1=j(ipb);
        if ~isempty(j1(salu(ipb)))
            Gij(:,:,j1(salu(ipb))) = Gij_3D_r_small_FP(coord,xf,yf,zf,j1(salu(ipb)),ks,kp,gaussian,C);
        end
        if ~isempty(j1(salt(ipb)))
            Tij(:,:,j1(salt(ipb))) = Tij_3D_r_small_FP(coord,xf,yf,zf,j1(salt(ipb)),ks,kp,gaussian,C);
        end
    end
    t(1,:)          = Tij(1,1,:)*fij(1)+Tij(1,2,:)*fij(2)+Tij(1,3,:)*fij(3);
    t(2,:)          = Tij(2,1,:)*fij(1)+Tij(2,2,:)*fij(2)+Tij(2,3,:)*fij(3);
    t(3,:)          = Tij(3,1,:)*fij(1)+Tij(3,2,:)*fij(2)+Tij(3,3,:)*fij(3);
else
    %calculo en los receptores
    t=[];
    ipb=find(r==0);
    if ~isempty(ipb(salu(ipb)))
        %cuando la fuente real coincide con el receptor, la solucion no
        %existe. Se mueve entonces la fuentes para un lado
        
        g(1,1)  = 0;
        g(2,1)  = 0;
        g(3,1)  = 0;
        lambda0 = min(abs(1/kp/5),1e-2);
        Gij(:,:,ipb)	= Gij_3D(ks,kp,lambda0,g(:,1),C,1);
    end
    if max(salt)
        [S_fx,S_fy,S_fz]= S_3D(r(salt),g(:,salt),ks,kp,fij);
        S1          = zeros(3,3,n);
        S1(:,:,:)   = S_fx*fij(1)+S_fy*fij(2)+S_fz*fij(3);
        sal      	= para.sortie;
        nss      	= sal.sxx + sal.syy + sal.szz+sal.sxy+sal.sxz+sal.syz;
        S           = zeros(nss,n);
        j   = 1;
        if sal.sxx; S(j,:,:) = S1(1,1,:,:);j=j+1; end;
        if sal.syy; S(j,:,:) = S1(2,2,:,:);j=j+1; end;
        if sal.szz; S(j,:,:) = S1(3,3,:,:);j=j+1; end;
        if sal.sxy; S(j,:,:) = S1(1,2,:,:);j=j+1; end;
        if sal.sxz; S(j,:,:) = S1(1,3,:,:);j=j+1; end;
        if sal.syz; S(j,:,:) = S1(2,3,:,:);       end;
    else
        S=0;
    end
end

u(1,:)          = Gij(1,1,:)*fij(1)+Gij(1,2,:)*fij(2)+Gij(1,3,:)*fij(3);
u(2,:)          = Gij(2,1,:)*fij(1)+Gij(2,2,:)*fij(2)+Gij(2,3,:)*fij(3);
u(3,:)          = Gij(3,1,:)*fij(1)+Gij(3,2,:)*fij(2)+Gij(3,3,:)*fij(3);