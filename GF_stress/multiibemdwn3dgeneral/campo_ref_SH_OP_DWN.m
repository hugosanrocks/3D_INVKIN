function [U,t,DWN]=campo_ref_SH_OP_DWN(iinc,coord,para,DWN,kxk,kzsigno)

% funcion que da los campos de desplazamientos y las tractiones
% para una onda plana SH incidente en un multi-estarto
% la onda proviene del semi espacio
% xf,zf coordenadades del origen de la onda (fase nula)
% #r relativo a los receptores
% se incluye de una ves los receptores reales para el calculo de los
% campos incidente

n       = coord.nbpt;
t       = zeros(1,n);
U       = zeros(1,n);

indi  	= coord.indm(1,1).ind;%true-false de los puntos perteneciendo a m=1

%posicion receptores virtuales (puntos de colocacion) y reales 
nrv     = DWN.nxrv;
xr      = DWN.xr;
zr0     = DWN.zr0;
izr0    = DWN.izr0;
salu    = DWN.salu;
sals    = DWN.sals;

%posicion fuente real
coordf.xs       = para.xs(iinc);
coordf.zs       = para.zs(iinc);

if abs(kxk)<=1
    %ondas homogeneas
    [u,S]   	= calcul_US_DWN_SH_Ncapas_HS_incOP(para,xr,zr0,izr0,salu,sals,coordf,kxk,kzsigno);
else
    %ondas heterogeneas
    [u,S]       = calcul_US_DWN_SH_Ncapas_HS_incOPI(para,xr,zr0,izr0,salu,sals,coordf,kxk,kzsigno);
    vS          = para.reg(1).sub(para.nsubmed).bet;
    k21         = kxk*para.reg(1).sub(para.nsubmed).ksi;
    w           = DWN.omegac;
    vLove       = 1/para.dkdw(iinc);%w/k21;
    peso        = sqrt(2/w*vS^2/vLove/2);% el ultimo /2 es para 2 incidencias
    u           = u*peso;
    S           = S*peso;

%     w           = DWN.omegac;
%     w0          = para.mode.w0;
% %     k2          = para.mode.k2;
% %     ikmax       = para.mode.ikmax;
% 
%     k20ok       = dispersion_curve_k2_fw(para,DWN,real(w));
%     
%     for imode=1:size(w0,2)
%         if w>w0(1,imode)
%             k21         = k20ok(imode);
%             
%             %             indk2       = find(real(w)<w0(1:ikmax(imode),imode),1);
%             %             dk2         = k2(indk2,imode)-k2(indk2-1,imode);
%             %             k20         = linspace(k2(indk2-1,imode),k2(indk2,imode),1000);
%             %             k21         = interp1(w0(1:ikmax(imode),imode),k2(1:ikmax(imode),imode),real(w),'pchip','extrap');
%             %             k20         = linspace(k21-dk2/5000,k21+dk2/5000,1000);
%             %             kxk         = real(k20/para.reg(1).sub(para.nsubmed).ksi);
%             %
%             %             for i=1:1000
%             %                 [u(:,i),S]   	= calcul_US_DWN_SH_Ncapas_HS_incOPI(para,xr,zr0,izr0,salu,sals,coordf,kxk(i),kzsigno);
%             %             end
%             kxk         = real(k21/para.reg(1).sub(para.nsubmed).ksi);
%             [u,S]     = calcul_US_DWN_SH_Ncapas_HS_incOPI(para,xr,zr0,izr0,salu,sals,coordf,kxk,kzsigno);
%             
%             vS      = para.reg(1).sub(para.nsubmed).bet;
%             vLove   = k21/w;
%             peso    = 2/w*vS^2/vLove/2;% el ultimo /2 es para 2 incidencias
%             
%             u       = u*peso;
%             S       = S*peso;
% 
%             %verification de la normalisation
%             %             zr0 = linspace(0,1.5,1000);
%             %             xr  = 0*zr0;
%             %             izr0= linspace(1,1000,1000);
%             %             salu= ones(1000,1);
%             %             sals= zeros(1000,1);
%             %             [u,S]   	= calcul_US_DWN_SH_Ncapas_HS_incOPI(para,xr,zr0,izr0,salu,sals,coordf,kxk,kzsigno);
%             %             E=sum(abs(u).^2*zr0(2));
%             %             figure(3);hold on;plot(zr0,abs(u).^2*zr0(2))
%             %             toto=0;
% 
%         end
%     end
end
DWN.uy0(:,iinc) = u(nrv+1:end);
DWN.s0(:,:,iinc)= S(DWN.inds,nrv+1:end);

U(indi)         = u(1:nrv);
S               = S(:,1:nrv);
% S = [SxyKW,SzyKW](xr,zr,w)

%calculo de la traccion en los puntos de colocacion solamente
vn(1,:)	= coord.vnx(indi);
vn(2,:)	= coord.vnz(indi);

t(indi)  = S(1,:).*vn(1,:)+S(2,:).*vn(2,:);
