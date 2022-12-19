function [Gc,Gp,Spz,Spx]=Gij_3D_DWN_polar(para,DWN,coordr,coordf,ic,mode)
%
%     Programa Gij3D-POLAR
%
%     Programa para calcular sismogramas sintéticos de Gij3D(X,0,t); Fuente
%     en el origen (0,0,0) y receptor en X=(x,y,z) en espacio completo.
%     Calculo a partir de la sol en coordenadas polares, Con integración en
%     el número de onda radial con Jo y J1,Calculo con frecuencia imaginaria
%     Se calculan sismogramas sinteticos en hasta 25 estaciones.
%
%                                                      (c)FJSS 2012
%
%mode==1, hace la integracion y considera los factores apropiados
%mode==2, no hace la suma ni agrega fgactores al fin de construir la matriz
%independientemente de la orientacion de la fuerza

%estratigrafia
wi          = DWN.omegac;

kr          = DWN.kr;
nk          = length(kr);
DK          = DWN.dkr;
kz          = zeros(2,nk);
DWN.kz      = kz;

% fuentes
xs          = coordf.xs;
ys          = coordf.ys;
zs          = coordf.zs;
nxs       	= length(xs);

% % convolucion con la geometria de la fuente
% mxspecs1=1;
% mxspecb1=1;
% mxspecs2=1;
% mxspecb2=1;
%
% if nxs==1
%     mxspec=1;
%     if isfield(coordf,'vnx')
%         lsegx=coordf.dr*coordf.vnz;
%         lsegz=coordf.dr*coordf.vnx;
%         mxspecs1=sinc((-lsegx*kr-lsegz*kzf(1,:))/2/pi);
%         mxspecb1=sinc((-lsegx*kr+lsegz*kzf(1,:))/2/pi);
%         mxspecs2=sinc((-lsegx*kr-lsegz*kzf(2,:))/2/pi);
%         mxspecb2=sinc((-lsegx*kr+lsegz*kzf(2,:))/2/pi);
%     end
% end

% receptores
xr          = coordr.x;
yr          = coordr.y;
zr          = coordr.z;
nrec        = length(xr);

if mode==1
    % inicialisacion tensor en referencia polar
    Gp          = zeros(3,3,nrec,nxs);%desplazamientos
    Spz         = zeros(3,3,nrec,nxs);%tractiones
    Spx         = zeros(3,3,nrec,nxs);%tractiones
    % Grr,Grt,Grz G11 G12 G13
    % Gtr,Gtt,Gtz G21 G22 G23
    % Gzr,Gzt,Gzz G31 G32 G33
    
    % inicialisacion tensor en referencia carthesiana
    Gc          = zeros(3,3,nrec,nxs);
    % Gxx,Gxy,Gxz G11 G12 G13
    % Gyx,Gyy,Gyz G21 G22 G23
    % Gzx,Gzy,Gzz G31 G32 G33
else
    Gp          = zeros(3,3,nrec,nxs,nk);%desplazamientos
    Spz         = zeros(3,3,nrec,nxs,nk);%tractiones
    Spx         = zeros(3,3,nrec,nxs,nk);%tractiones
    Gc          = 0;
end

for is=1:nxs
    ksi         = para.reg(1).sub(ic(is)).ksi;
    kpi         = para.reg(1).sub(ic(is)).kpi;
    kzS         = sqrt(ksi.^2-kr.^2);
    kzP         = sqrt(kpi.^2-kr.^2);
    kzS         =-kzS.*(sign(imag(kzS)).*(imag(kzS)~=0)+(imag(kzS)==0));
    kzP       	=-kzP.*(sign(imag(kzP)).*(imag(kzP)~=0)+(imag(kzP)==0));
    C66         = para.reg(1).sub(ic(is)).Ci(6,6);
    
    %la depandancia en x0 e y0 solo aparece fuera de la suma por lo que se
    %puede agrupar los terminos con el mismo z0
    for ir=1:nrec
        z0      = zr(ir)-zs(is);
        expSz 	= exp(-1i*kzS*abs(z0));
        expPz   = exp(-1i*kzP*abs(z0));

        if mode==1
            x0      = xr(ir)-xs(is);
            y0      = yr(ir)-ys(is);
            r       = sqrt(x0^2+y0^2);
            J0      = besselj(0,kr*r);
            J1      = besselj(1,kr*r);
            
            if r<1e-6
                ct  = 1/sqrt(2);
                st  = 1/sqrt(2);
            else
                ct  = x0/r;
                st  = y0/r;
            end
        end
        
        %%%%%%
        % Fz %
        %%%%%%
        if mode==1
            Grz             =-sign(z0)*sum((expPz-expSz).*kr.^2.*J1.*DK);% P&SV
            
            Gp(1,3,ir,is)   = Grz;
            %Gp(2,3)=0;
            Gc(1,3,ir,is)   = Grz*ct;%Gxz
            Gc(2,3,ir,is)   = Grz*st;%Gyz
            
            Gc(3,3,ir,is) 	=-1i*sum((kzP.*expPz+kr.^2./kzS.*expSz).*kr.*J0.*DK);% Gzz P&SV
            Gp(3,3,ir,is)   = Gc(3,3,ir,is);
            
            Spz(3,3,ir,is)  = sign(z0)*C66*sum(((kr.^2-kzS.^2).*expPz-2*kr.^2.*expSz).*kr.*J0.*DK);% P&SV
            Spz(3,1,ir,is)  = 1i*C66*sum((2*kzP.*expPz+(kr.^2-kzS.^2)./kzS.*expSz).*kr.^2.*J1.*DK);% P&SV
        else
            Gp(1,3,ir,is,:) = sign(z0)*(expPz-expSz).*kr;% P&SV Grz; (-k J1)
            %Gp(2,3)=0;
            
            Gp(3,3,ir,is,:) =-1i*(kzP.*expPz+kr.^2./kzS.*expSz).*kr;% Gzz P&SV (J0) 
            
            Spz(3,3,ir,is,:)= C66*sign(z0)*((kr.^2-kzS.^2).*expPz-2*kr.^2.*expSz).*kr;% P&SV (J0)
            Spz(3,1,ir,is,:)= C66*1i*(2*kzP.*expPz+(kr.^2-kzS.^2)./kzS.*expSz).*kr;% P&SV (k J1)
        end
        
        %%%%%%
        % Fx %
        %%%%%%
        termPSV0        = (kr.^2./kzP.*expPz+kzS.*expSz);
        termSH0         = ksi.^2./kzS.*expSz;
        
        if mode==1
            if r<1e-6
                fac2= 0.5;
                fac3= 0.5;
            else
                fac3= J1./(kr*r);
                fac2= J0-fac3;
            end
            
            
            % Ur and Ut
            Gr0             =-1i*sum((termPSV0.*fac2 + termSH0.*fac3).*kr.*DK);            % P&SV y SH 
%             Gp(1,1,ir,is)   = Gr0;
            Grx             = Gr0 * ct;%projection de la force 
            Gp(1,1,ir,is)   = Grx;
            Gry             = Gr0 * st;
            Gp(1,2,ir,is)   = Gry;
            
            Gtr             = 1i*sum((termPSV0.*fac3 + termSH0.*fac2).*kr.*DK);             % P&SV y SH 
%             Gp(2,1,ir,is)   = Gtr;
            Gtx             = Gtr * st;
            Gp(2,1,ir,is)   = Gtx;
            Gty             = Gtr *-ct;
            Gp(2,2,ir,is)   = Gty;
            
            % Ux and Uy projection des deplacements
            Gc(1,1,ir,is)   = Grx*ct-Gtx*st;%Gxx
            Gc(2,1,ir,is)   = Grx*st+Gtx*ct;%Gyx
            
            Gc(1,2,ir,is)   = Gry*ct-Gty*st;%Gxy
            Gc(2,2,ir,is)   = Gry*st+Gty*ct;%Gyy
            
            % Uz:
            Gzr             =-sign(z0)*sum((expPz-expSz).*kr.^2.*J1.*DK);  % P&SV  
            Gp(3,1,ir,is)   = Gzr;
            Gc(3,1,ir,is)   = Gzr * ct;%Gzx
            Gc(3,2,ir,is)   = Gzr * st;%Gzy
            
            Spx(3,1,ir,is)  =-sign(z0)*C66*sum(((2*kr.^2.*expPz-(kr.^2-kzS.^2).*expSz).*(J0-fac3) ...
                +ksi.^2.*expSz.*fac3).*kr.*DK)*ct;% P&SV & SH
            Spx(3,2,ir,is)  = sign(z0)*C66*sum(((2*kr.^2.*expPz-(kr.^2-kzS.^2).*expSz).*fac3 ...
                +ksi.^2.*expSz.*(J0-fac3)).*kr.*DK)*st;% P&SV & SH
            Spx(3,3,ir,is)  =-1i*C66*sum(((kr.^2-kzS.^2)./kzP.*expPz+2*kzS.*expSz).*kr.^2.*J1.*DK)*ct;% P&SV
        else
            %en este caso se ocupa Gp para el caso PSV
            % y Gp(2,2) para SH, los dos en polares
            %tambien solo se calcula las componentes necesarias a la
            %resolucion del sistema de ecuaciones
            Gp(1,1,ir,is,:)	=-1i*termPSV0;                                  % P&SV Grx ok
            Gp(3,1,ir,is,:) =-sign(z0)*(expPz-expSz).*kr.^2;                % P&SV Gzx ok
            
            Spx(3,1,ir,is,:)= C66*sign(z0)*(2*kr.^2.*expPz-(kr.^2-kzS.^2).*expSz);      % Szr P&SV ok
            Spx(3,3,ir,is,:)=-1i*C66*((kr.^2-kzS.^2)./kzP.*expPz+2*kzS.*expSz).*kr.^2;  % Szz P&SV ok
            
            Gp(2,2,ir,is,:)	=-1i*termSH0;            % SH Grr ok
            Spx(2,2,ir,is,:)= C66*sign(z0)*ksi.^2.*expSz;% SH Szr ok 
        end
    end
    if mode==1
        Gc(:,:,:,is,:)=Gc(:,:,:,is,:)/(4*pi*wi^2*para.reg(1).sub(ic(is)).rho);
    end
    Gp(:,:,:,is,:) = Gp(:,:,:,is,:) /(4*pi*wi^2*para.reg(1).sub(ic(is)).rho);
    Spx(:,:,:,is,:)= Spx(:,:,:,is,:)/(4*pi*wi^2*para.reg(1).sub(ic(is)).rho);
    Spz(:,:,:,is,:)= Spz(:,:,:,is,:)/(4*pi*wi^2*para.reg(1).sub(ic(is)).rho);
end
