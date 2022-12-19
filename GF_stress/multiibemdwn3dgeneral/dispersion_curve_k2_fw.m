function [k20,dkdw]=dispersion_curve_k2_fw(para,DWN,w)

% a partir de las curbas de dispersion, busca precisamente el k2
% correspndiendo al w dado

% tic

para.DWNomei= 0;
% for ms=1:para.nsubmed
%     para.reg(1).sub(ms).tipoatts=0;
% end
% para=Vp2Cij(para);
%
% nmode   = 9;
% k2      = zeros(1,nmode);
% w0      = zeros(1,nmode);
% err1    = zeros(1,nmode);
% ikmax 	= zeros(1,nmode);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % pour k=0, il n y a que des modes de volume,c est une oscillation facile a trouver
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C       = para.reg(1).sub(1).C;
% rho     = para.reg(1).sub(1).rho;
% vT1     = sqrt(C(6,6)/rho);
% C       = para.reg(1).sub(2).C;
% rho     = para.reg(1).sub(2).rho;
% vT2     = sqrt(C(6,6)/rho);

w0      = para.mode.w0;
k2      = para.mode.k2;
ikmax   = para.mode.ikmax;
nmode   = size(w0,2);
k20     = zeros(nmode,1);
dkdw     = zeros(nmode,1);
for imode=1:nmode
    if w>w0(1,imode)
%         disp(imode);
        
        %approx de la solution par interpolation
        kguess  = interp1(w0(1:ikmax(imode),imode),k2(1:ikmax(imode),imode),real(w),'pchip','extrap');
        
        indw    = find(w0(1:ikmax(imode),imode)>w,1);
        dw      = 1/1000*(w0(indw,imode)-w0(indw-1,imode));
        kp      = interp1(w0(1:ikmax(imode),imode),k2(1:ikmax(imode),imode),real(w+dw),'pchip','extrap');
        km      = interp1(w0(1:ikmax(imode),imode),k2(1:ikmax(imode),imode),real(w-dw),'pchip','extrap');
        k20(imode)   = kguess;
%         if imode>1 && imode<nmode
%             dkp = interp1(w0(1:ikmax(imode-1),imode-1),k2(1:ikmax(imode-1),imode-1),real(w),'pchip','extrap');
%             dkp = (dkp-kguess)/1000;
%             if w>w0(1,imode+1)
%                 dkm = interp1(w0(1:ikmax(imode+1),imode+1),k2(1:ikmax(imode+1),imode+1),real(w),'pchip','extrap');
%                 dkm = (kguess-dkm)/1000;
%             else
%                 dkm=dkp;
%             end
%         elseif imode==1
%             if w>w0(1,imode+1)
%                 dkm = interp1(w0(1:ikmax(imode+1),imode+1),k2(1:ikmax(imode+1),imode+1),real(w),'pchip','extrap');
%                 dkm = (kguess-dkm)/1000;
%             else
%                 dkm=kguess/1000;
%             end
%             dkp=dkm;
%         elseif imode==nmode
%             %keep the same as imode=nmode-1
%         end
%         
%         paratmp=para;
%         paratmp.DWNomei=0;
%         DWN.omegac  = real(w);
%         paratmp1    = attenuation(paratmp,real(w)/(2*pi));
%         k20(imode)	= fzero(@(k) zerodeAk(paratmp1,DWN,k),[kguess-dkm kguess+dkp]);
%         
%         DWN.omegac  = real(w+dw);
%         paratmp1    = attenuation(paratmp,real(w+dw)/(2*pi));
%         kp      	= fzero(@(k) zerodeAk(paratmp1,DWN,k),[kp-dkm     kp+dkp]);
%         
%         DWN.omegac  = real(w-dw);
%         paratmp1    = attenuation(paratmp,real(w-dw)/(2*pi));
%         km          = fzero(@(k) zerodeAk(paratmp1,DWN,k),[km-dkm     km+dkp]);

        dkdw(imode)	= (kp-km)/(2*dw);

        
        %     %intervale de recherche
        %     if imode==1
        %         ksup= w/vT1;
        %
        %     else
        %
        %     end
        %
        %     %verifier que la partie reelle et imaginaire change de signe
        %     ki  = kguess-dk;
        %     kf  = kguess+dk;
        %     [sign0,sign1]=checkchgmtsign(wi,wf,para,DWN);
        %
        %     %pb qd l intervale de recherche est trop petit
        %     j=1;
        %     while sign0>=0 && sign1>=0 && j<6
        %         % plotthepb(wi,wf,para,DWN)
        %
        %         j   = j+1;
        %         dw0 = dw0*2;
        %         wi  = nextw0-dw0;
        %         wf  = nextw0+dw0;
        %         [sign0,sign1]=checkchgmtsign(wi,wf,para,DWN);
        %     end
        %     %pb qd l intervale de recherche est trop grand et voit 2 zeros
        %     while (sign0>=0 || sign1>=0) && j<6
        %         % plotthepb(wi,wf,para,DWN)
        %         if sign0<0
        %             dw=wf-wi;
        %             [sign0]=checkchgmtsign(wi,wi+dw/2,para,DWN);
        %             if sign0<0
        %                 nextw0=wi+dw/4;
        %             else
        %                 nextw0=wi+3*dw/4;
        %             end
        %         else
        %             dw0 = dw0*4;
        %         end
        %         j   = j+1;
        %         dw0 = dw0/2;
        %         wi  = nextw0-dw0;
        %         wf  = nextw0+dw0;
        %         [sign0,sign1]=checkchgmtsign(wi,wf,para,DWN);
        %     end
        %     if sign0<=0 && sign1<=0
        %         ik  = ik+1;
        %         w0(ik1,imode)   = fzero(@(w) zerodeAw(para,DWN,w),[wi wf]);
        %         DK              = min(1.3*DK,DK1);
        %         err             = abs(w0(ik1,imode)-nextw0);
        %         err1(ik1,imode) = err;
        %     else
        %         %si rien n a marche avant, on diminue le pas d avance
        %         DK              = DK/2;
        %     end
    end
end
% figure(207);hold on
% plot(real(w)+dkdw*0,1/dkdw,'.')
% toc

% figure(205);hold on;
% for imode=1:nmode
%     plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode)- 2*err1(1:ikmax(imode),imode),'r')
%     plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode),'k')
%     plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode)+ 2*err1(1:ikmax(imode),imode),'r')
% end
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                    cambio a grafica vp=w/k - w                          %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(206);hold on;
% for imode=1:nmode
%     vp=w0(1:ikmax(imode),imode)./k2(1:ikmax(imode),imode);
%     if imode==1
%         vp(1)=vT2;
%     end
%     plot(w0(1:ikmax(imode),imode)/2/pi,vp,'k')
% end