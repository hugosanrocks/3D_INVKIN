function [u,t,sout]=campo_ref_PSV_FP(xf,zf,coord,kp,ks,C,fij,salu,salt,para)
% funcion que da los campos de desplazamientos y las tractiones
% para una fuente puntual direcional (real) en un espacio completo
% xf,zf coordenadades de la fuente
% #r relativo a los receptores

n       = coord.nbpt;
gaussian   = para.gaussian;

u       = zeros(2,n);
t       = zeros(2,n);
sout       = zeros(3,n);
Gij     = zeros(2,2,n);
Tij     = zeros(2,2,n);


xr  	= coord.x;
zr  	= coord.z;

x       = xr-xf;
z       = zr-zf;
r       = sqrt(x.^2+z.^2);
npt     = size(x);

if C(6,6)~=0
    % a dentro de un medio solido
    g(1,:)  = x./r;
    g(2,:)  = z./r;
    
    Gij(:,:,salu)	= Gij_PSV(ks,kp,r(salu),g(:,salu),C,sum(salu));
    
    % tensor de esfuerzos
    vn(1,1:npt) = 1; vn(2,1:npt) = 0;
    Tij(:,:,salt) = Tij_PSV(ks,kp,r(salt),g(:,salt),C,vn(:,salt));
    if isfield(coord,'dr')
    ipb=find(r<2*para.npplo*coord.dr);
    if ~isempty(ipb)
        j=1:n;
        j1=j(ipb);
        coordaux = coord; coordaux.vnx=vn(1,1:npt); coordaux.vnz=vn(2,1:npt);
        if ~isempty(j1(salu(ipb)))
            Gij(:,:,j1(salu(ipb))) = Gij_PSV_r_small_FP(coordaux,xf,zf,j1(salu(ipb)),ks,kp,gaussian,C);
        end
        if ~isempty(j1(salt(ipb)))
            Tij(:,:,j1(salt(ipb))) = Tij_PSV_r_small_FP(coordaux,xf,zf,j1(salt(ipb)),ks,kp,gaussian,C);
        end
    end
    end
    k = 1; sal = para.sortie;
    if sal.sxx==1
    sout(k,:) = Tij(1,1,:)*fij(1)+ Tij(1,2,:)*fij(2); k=k+1;%sxx
    end
    if sal.sxz==1
    sout(k,:) = Tij(2,1,:)*fij(1)+ Tij(2,2,:)*fij(2); k=k+1; %sxz
    end
    if sal.szz==1
    vn(1,1:npt) = 0; vn(2,1:npt) = 1;
    Tij(:,:,salt) = Tij_PSV(ks,kp,r(salt),g(:,salt),C,vn(:,salt));
    if exist('ipb','var')
    if ~isempty(ipb)
        coordaux.vnx=vn(1,1:npt); coordaux.vnz=vn(2,1:npt);
        if ~isempty(j1(salt(ipb)))
            Tij(:,:,j1(salt(ipb))) = Tij_PSV_r_small_FP(coordaux,xf,zf,j1(salt(ipb)),ks,kp,gaussian,C);
        end
    end
    end
    sout(k,:) = Tij(2,1,:)*fij(1)+ Tij(2,2,:)*fij(2); %szz
    end
    
    if isfield(coord,'vnx')
        %calculo del vector de B en el IBEM, fuerzas y desplazamientos debidos
        %a las fuerzas reales en los puntos de colocacion
        vn(1,:) = coord.vnx;
        vn(2,:) = coord.vnz;
        
        Tij(:,:,salt)	= Tij_PSV(ks,kp,r(salt),g(:,salt),C,vn(:,salt));
        
        ipb=find(r<2*para.npplo*coord.dr);
        
        if ~isempty(ipb)
            %         [~,ipb]=min(r);
            j=1:n;
            j1=j(ipb);
            if ~isempty(j1(salu(ipb)))
                Gij(:,:,j1(salu(ipb))) = Gij_PSV_r_small_FP(coord,xf,zf,j1(salu(ipb)),ks,kp,gaussian,C);
            end
            if ~isempty(j1(salt(ipb)))
                Tij(:,:,j1(salt(ipb))) = Tij_PSV_r_small_FP(coord,xf,zf,j1(salt(ipb)),ks,kp,gaussian,C);
            end
        end
        t(1,:)          = Tij(1,1,:)*fij(1)+Tij(1,2,:)*fij(2);
        t(2,:)          = Tij(2,1,:)*fij(1)+Tij(2,2,:)*fij(2);
        u(1,:)          = Gij(1,1,:)*fij(1)+Gij(1,2,:)*fij(2);
        u(2,:)          = Gij(2,1,:)*fij(1)+Gij(2,2,:)*fij(2);
    else
        %calculo de la contribucion de la fuente real en los receptores
        
        t=[];
        sal     = para.sortie;
        ns      = (sal.Ux+sal.Uz)*(sal.Ut+sal.UPh+sal.USh+sal.UIh+sal.UPt+sal.USt);
        u       = zeros(ns,n);
        k       = 1;
        
        ipb=find(r==0);
        if ~isempty(ipb(salu(ipb)))
            %cuando la fuente real coincide con el receptor, la solucion no
            %existe. Se mueve entonces la fuentes para un lado
            g(1,1)      = 0;
            g(2,1)      = 0;
            lambda0     = min(abs(1/kp/5),1e-2);
            Gij(:,:,ipb)= Gij_PSV(ks,kp,lambda0,g(:,1),C,1);
        end
        
        
        if sal.Ut==1
            if sal.Ux==1
                u(k,:)          = Gij(1,1,:)*fij(1)+Gij(1,2,:)*fij(2);k=k+1;
            end
            if sal.Uz==1
                u(k,:)          = Gij(2,1,:)*fij(1)+Gij(2,2,:)*fij(2);k=k+1;
            end
        end
        
        if sal.UPh==1 || sal.USh==1 || sal.UIh==1
            % Separacion del campo de desplazamiento total en 3 contribuciones: P
            % homogeneo, SV homogeneo, y ondas inhomogeneas (k>kp y >ks).
            % Se obtiene a traves de la descomposicion de las funciones de Hankel
            % en ondas planas homogeneas y inhomogeneas
            [Gpr,Gsr]	= Gij_P_SV_R(ks,kp,x(salu),z(salu),r(salu),g(:,salu),C,zeros(sum(salu),1));
            Gi      = Gij-(Gpr+Gsr);
            if sal.UPh==1
                if sal.Ux==1
                    u(k,:)          = Gpr(1,1,:)*fij(1)+Gpr(1,2,:)*fij(2);k=k+1;
                end
                if sal.Uz==1
                    u(k,:)          = Gpr(2,1,:)*fij(1)+Gpr(2,2,:)*fij(2);k=k+1;
                end
            end
            if sal.USh==1
                if sal.Ux==1
                    u(k,:)          = Gsr(1,1,:)*fij(1)+Gsr(1,2,:)*fij(2);k=k+1;
                end
                if sal.Uz==1
                    u(k,:)          = Gsr(2,1,:)*fij(1)+Gsr(2,2,:)*fij(2);k=k+1;
                end
            end
            if sal.UIh==1
                if sal.Ux==1
                    u(k,:)          =  Gi(1,1,:)*fij(1)+ Gi(1,2,:)*fij(2);k=k+1;
                end
                if sal.Uz==1
                    u(k,:)          =  Gi(2,1,:)*fij(1)+ Gi(2,2,:)*fij(2);k=k+1;
                end
            end
        end
        if sal.UPt==1 || sal.USt==1
            % Separacion de los campos de polarizacion P y S
            Gpt	= Gij_Pt_SVt(kp,r(salu),g(:,salu),C,sum(salu));
            Gst = Gij-Gpt;
            if sal.UPt==1
                if sal.Ux==1
                    u(k,:)          = Gpt(1,1,:)*fij(1)+Gpt(1,2,:)*fij(2);k=k+1;
                end
                if sal.Uz==1
                    u(k,:)          = Gpt(2,1,:)*fij(1)+Gpt(2,2,:)*fij(2);k=k+1;
                end
            end
            if sal.USt==1
                if sal.Ux==1
                    u(k,:)          = Gst(1,1,:)*fij(1)+Gst(1,2,:)*fij(2);k=k+1;
                end
                if sal.Uz==1
                    u(k,:)          = Gst(2,1,:)*fij(1)+Gst(2,2,:)*fij(2);
                end
            end
        end
    end
else
    if isfield(coord,'vnx')
        %calculo del vector de B en el IBEM, fuerzas y desplazamientos debidos
        %a las fuerzas reales en los puntos de colocacion
        vnxi            = coord.vnx;
        vnzi            = coord.vnz;
        drdn            = (x(salu).*vnxi(salu)+z(salu).*vnzi(salu))./r(salu);
        
        Gij(1,1,salu)	= 1i/4*kp.*besselh(1,2,kp.*r(salu)).*drdn;
        
        Tij(1,1,salt)	=-C(1,1)*kp^2*G22_SH(kp,r(salt),1);
        
        
        %         ipb=find(r<2*coord.dr);
        %
        %         if ~isempty(ipb)
        %             %         [~,ipb]=min(r);
        %             j=1:n;
        %             j1=j(ipb);
        %             if ~isempty(j1(salu(ipb)))
        %                 Gij(:,:,j1(salu(ipb))) = Gij_PSV_r_small_FP(coord,xf,zf,j1(salu(ipb)),ks,kp,gauss,C);
        %             end
        %             if ~isempty(j1(salt(ipb)))
        %                 Tij(:,:,j1(salt(ipb))) = Tij_PSV_r_small_FP(coord,xf,zf,j1(salt(ipb)),ks,kp,gauss,C);
        %             end
        %         end
        
        u(1,:)	= Gij(1,1,:);%un
        t(1,:)	= Tij(1,1,:);%tn
    else
        %calculo de la contribucion de la fuente real en los receptores
        
        t=[];
        sal     = para.sortie;
        ns      = (sal.Ux+sal.Uz)*(sal.Ut+sal.UPh+sal.USh+sal.UIh+sal.UPt+sal.USt);
        u       = zeros(ns,n);
        k       = 1;

        Gij0            = 1i/4*kp.*besselh(1,2,kp.*r);
        if sal.Ut==1
            if sal.Ux==1
                u(k,:) 	= Gij0.*x./r;k=k+1;
            end
            if sal.Uz==1
                u(k,:)	= Gij0.*z./r;k=k+1;
            end
        end
        
        if sal.UPh==1 || sal.USh==1 || sal.UIh==1
            % Separacion del campo de desplazamiento total en 3 contribuciones: P
            % homogeneo, SV homogeneo, y ondas inhomogeneas (k>kp y >ks).
            % Se obtiene a traves de la descomposicion de las funciones de Hankel
            % en ondas planas homogeneas y inhomogeneas

            if sal.UPh==1
                if sal.Ux==1
                    u(k,:)          = Gij0.*x./r;k=k+1;
                end
                if sal.Uz==1
                    u(k,:)          = Gij0.*z./r;k=k+1;
                end
            end
            if sal.USh==1
                if sal.Ux==1
                    u(k,:)          = 0;k=k+1;
                end
                if sal.Uz==1
                    u(k,:)          = 0;k=k+1;
                end
            end
            if sal.UIh==1
                if sal.Ux==1
                    u(k,:)          = 0;k=k+1;
                end
                if sal.Uz==1
                    u(k,:)          = 0;k=k+1;
                end
            end
        end
        if sal.UPt==1 || sal.USt==1
            % Separacion de los campos de polarizacion P y S
            if sal.UPt==1
                if sal.Ux==1
                    u(k,:)          = Gij0.*x./r;k=k+1;
                end
                if sal.Uz==1
                    u(k,:)          = Gij0.*z./r;k=k+1;
                end
            end
            if sal.USt==1
                if sal.Ux==1
                    u(k,:)          = 0;k=k+1;
                end
                if sal.Uz==1
                    u(k,:)          = 0;
                end
            end
        end
    end
end
