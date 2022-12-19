nf      = para.nf;
df      = para.fmax/(para.nf/2);     %paso en frecuencia
Fq      = (0:nf/2)*df;
nfN     = nf/2+1; %Nyquist
zerospad= para.zeropad;
tps     = 0:(1/(df*2*(nfN+zerospad))*...
          (2*(nfN+zerospad)/(2*(nfN+zerospad)-2))):1/df;
tps     = para.pulso.b+tps;

[n1,n2,n3,n4]=size(utc);
nt=n1;
listinc={'P','S','R'};
nrecx=para.rec.nrecx;
nrecy=para.rec.nrecy;
nrecz=para.rec.nrecz;

xr      = para.rec.xri+para.rec.dxr*(0:(nrecx-1));
zr      = para.rec.zri+para.rec.dzr*(0:(nrecz-1));
Ut      = zeros(para.rec.nrecx,para.rec.nrecz);
% Up     = zeros(para.rec.nrecx,para.rec.nrecz);
% Us     = zeros(para.rec.nrecx,para.rec.nrecz);
% Ui     = zeros(para.rec.nrecx,para.rec.nrecz);

m3max=max(max(max(max(utc))));
clear F
fig=figure;
set(fig,'DoubleBuffer','on');
name=[para.nomrep,para.nomrep(1),'filmtmp','_000'];
while exist([name,'.avi'],'file')==2
    compt=str2double(name(length(name)-2:length(name)));
    name=[name(1:length(name)-3),num2str(compt+1,'%.3i')];
end

idt = 10;
i0  = 1;
iend= 800;
%mov = avifile(name,'compression','none');
mov = VideoWriter(name);%,'LosslessCompression','true');
open(mov);

for iinc=1
%     name=['C:\Users\mperton\Desktop\TMP_calculos\Cuna_inc',num2str(iinc),'_'];
%  name='C:\Users\mperton\Desktop\TMP_calculos\diff_field';
%  name='C:\Users\Mathieu\Desktop\TMP_calcul\current';
    hold off
    if para.pol==1
        u  = sign(utc(i0:idt:iend,:,iinc,1)).*log(abs(utc(i0:idt:iend,:,iinc,1)).*1e3+1);
%         us = sign(utc(i0:idt:iend,:,iinc,2)).*log(abs(utc(i0:idt:iend,:,iinc,2)).*1e6+1);
%         ui = sign(utc(i0:idt:iend,:,iinc,3)).*log(abs(utc(i0:idt:iend,:,iinc,3)).*1e6+1);
    else
%        u =sign(utc(1:idt:end,:,iinc,1)).*log(abs(utc(1:idt:end,:,iinc,2)).*1e7+1);

         u =log((sqrt(utc(1:idt:end,:,iinc,1).^2+utc(1:idt:end,:,iinc,2).^2)).*1e9+1);
%         up=log((sqrt(utc(1:idt:end,:,iinc,3).^2+utc(1:idt:end,:,iinc,4).^2)).*1e5+1);
%         us=log((sqrt(utc(1:idt:end,:,iinc,5).^2+utc(1:idt:end,:,iinc,6).^2)).*1e5+1);
%         ui=log((sqrt(utc(1:idt:end,:,iinc,7).^2+utc(1:idt:end,:,iinc,8).^2)).*1e5+1);
    end
    
    m3max=max(max(max(max(abs(u)))));
    
    iit=0;
    [ntf,nes]=size(u);
    for it= 1:ntf%)%1:20:ntf%(1:3:ntf)%[(4600:idt:nt) ]
        iit=iit+1;
        for iz=1:nrecz
            for iy=1:nrecy
                for ix=1:nrecx
                    ies         = ix+(iy-1)*nrecx+(iz-1)*nrecx*nrecy;
                    if max(abs(u(:,ies)))==0
                        Ut(ix,iz)   = nan;
                    else
                        Ut(ix,iz)   = u(it,ies) /m3max;%min((u(it,ies) /m3max),1);
                    end
                    
                    %                     Up(ix,iz)  = 1.9*(up(it,ies)/m3max)-.9;
%                     Us(ix,iz)  = us(it,ies)/m3max;
%                     Ui(ix,iz)  = ui(it,ies)/m3max;
                end
            end
        end
        Ut(1,1)   = 1;
        Ut(2,1)   =-1;
        %         Up(1,1)  = 1;
        %         Up(2,1)  =-1;
%         Us(1,1)  = 1;
%         Us(2,1)  =-1;
%         Ui(1,1)  = 1;
%         Ui(2,1)  =-1;
        Ut(isnan(Ut))=0;
        %         a1=subplot(2,2,1);
          tmph = surf(xr,zr,Ut.');hold on;plot_contour;title('||U||')
%         a1 = subplot(3,1,1); tmph = surf(xr,zr,Ut.');hold on;plot_contour;title('u_2')
%         a2 = subplot(3,1,2); tmph= surf(xr,zr,Us.');hold on;plot_contour;title('u_2 S')
%         a3 = subplot(3,1,3); tmph= surf(xr,zr,Ui.');hold on;plot_contour;title('u_2 evanescentes')
        %         a2=subplot(2,2,2);tmphp= surf(zr,xr,Up);hold on;plot_contour;title('||u|| P')
        %         a3=subplot(2,2,3);tmphs= surf(zr,xr,Us);hold on;plot_contour;title('||u|| S')
        %         a4=subplot(2,2,4);tmphi= surf(zr,xr,Ui);hold on;plot_contour;title('||u|| evanescente')
        
        
%         h=annotation('textbox',[.35 .0 .3 .1]);
%         set(h,'FitHeightToText','on',...
%             'string',['onda ',listinc{para.tipo_onda(iinc)},' incidente con un angulo de ',num2str(para.gam(iinc)),'º'],...
%             'HorizontalAlignment','center','FontWeight','bold','LineStyle','none')
%         
        

        %         F(iit) = getframe;
        %         F = getframe;
        %         mov = addframe(mov,F);
        tmpstr=num2str(iit);
        while length(tmpstr)<floor(log10(ntf/idt)+1)
            tmpstr=['0',tmpstr];
        end
        frame = getframe;
        writeVideo(mov,frame);
%         saveas(gcf,[name,tmpstr],'jpg')
        pause(.1)
    end
end
% mov = close(mov);
close(mov);


% movie2avi(F,name,'compression','Indeo5')
% movie(F)