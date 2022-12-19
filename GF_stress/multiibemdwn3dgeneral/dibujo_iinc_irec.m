%% título
titledibujo;
titletxt=[titletxt,'at receptor ',num2str(irec),' in xr:',num2str(xr(irec))];
if para.dim >= 3
    titletxt=[titletxt,', yr:',num2str(yr(irec))];
end
titletxt    =[titletxt,', zr:',num2str(zr(irec))];

set(h101(ifig),'name',titletxt,'numbertitle','off');
figure(h101(ifig));

tmpc    =get(bouton.couleur,'string');
coul    =get(bouton.couleur,'value');
if holdtr==1
    itmpc=mod(coul+iinc-1+irec-1,length(tmpc));
    if itmpc==0;itmpc=length(tmpc);end
    couleur=tmpc{itmpc};
else
    couleur =tmpc{coul};
end

if fieldV(ifig).name(1)=='U'
    varplot=utc;
    varplotw=uw;
    ifig0=ifig;
elseif fieldV(ifig).name(1)=='S'
    varplot=stc;
    varplotw=sw;
    ifig0=1;
end

% senal de una posicion
if para.spct==0
  ha = zeros(1,fieldV(ifig).nc);
    for i=1:fieldV(ifig).nc
        if para.b_dib(1).dessinsp==1
            subplot(2,fieldV(ifig).nc,i)
            if holdtr==2;hold off;else hold on;end
            h=plot(tps,squeeze(varplot(:,irec,iinc,i+(ifig0-1)*fieldV(ifig).nc)),couleur);
            set(h,{'DisplayName'},{['t_R_',num2str(irec),'_I',num2str(iinc)]});ha(i)=gca;
        else
            subplot(1,fieldV(ifig).nc,i)
            if holdtr==2;hold off;else hold on;end
            h=plot(tps,squeeze(varplot(:,irec,iinc,i+(ifig0-1)*fieldV(ifig).nc)),couleur);
            set(h,{'DisplayName'},{['t_R_',num2str(irec),'_I',num2str(iinc)]});ha(i)=gca;
        end
        xlabel('time')
        ylabel('amplitude')
        tmpc0= fieldV(ifig).namec{i};
        title([fieldV(ifig).name,'_{',tmpc0(2:end),'}']);
    end
    linkaxes(ha,'x'); clear ha h
end

% espectro de una posicion
if para.b_dib(1).dessinsp==1
  ha = zeros(1,fieldV(ifig).nc);
    cspectre    =correction_spectre(para,nfN,df);%pour supprimer la convolution temporelle de la source
    for i=1:fieldV(ifig).nc
        if para.spct==0
            uisubp = subplot(2,fieldV(ifig).nc,fieldV(ifig).nc+i);
        else
            uisubp = subplot(1,fieldV(ifig).nc,i);
        end
        if holdtr==2;hold off;else hold on;end
%         plot(Fq,imag(varplotw(1:nfN,irec,iinc,i+(ifig0-1)*fieldV(ifig).nc)),...
%                [couleur,'--'],'LineWidth',2)
%         hold on;
        espe = varplotw(1:nfN,irec,iinc,i+(ifig0-1)*fieldV(ifig).nc).*cspectre.';
        hr=plot(Fq,real(espe),[couleur '-']); hold on
        hi=plot(Fq,imag(espe),[couleur '--']);
        set(hr,{'DisplayName'},{['f_R_',num2str(irec),'_I',num2str(iinc)]});
        set(hi,{'DisplayName'},{['f_R_',num2str(irec),'_I',num2str(iinc)]});ha(i)=gca;
%         legend('pulse response','source response');
        xlabel('frequency')
        ylabel('amplitude')
    end
    linkaxes(ha,'x'); clear ha h espe
end

%H/V 2D
% para.ninc   = 2;
% para.gam    = [90   180];
% para.fuente = 2;
% para.xs     = [0     0];
% para.zs     = [1     1]*1e-6;
% para.pol    = 2;
% figure(108);hold on
% if para.dim==1 && para.pol==2
%     plot(Fq,sqrt(abs(imag(uw(:,1,1,1))./imag(uw(:,1,2,2)))))
%     plot(Fq,sqrt(abs(imag(uw(:,2,3,1))./imag(uw(:,2,4,2)))))
%     plot(Fq,sqrt(abs(imag(uw(:,3,5,1))./imag(uw(:,3,6,2)))))
% else
%     plot(Fq,sqrt(abs(2*imag(uw(:,1,1,1))./imag(uw(:,1,2,3)))))
% end