figure(h201(ifig));
%%%%%%%%%%%%%%%%%%%%%%%
% dibujo de la fuente %
%%%%%%%%%%%%%%%%%%%%%%%
if para.fuente==1
    para.b_dib(ifig).xzs=plot(para.b_dib(ifig).hg,para.xs(iinc),para.zs(iinc),'r.');
    xarrow  =-.25*cos(para.gamr(iinc));
    yarrow  =-.25*sin(para.gamr(iinc));
    x       = [xarrow -xarrow];
    y       = [yarrow -yarrow];
    para.b_dib(ifig).xzb= plot(para.b_dib(ifig).hg,x+para.xs(iinc),y+para.zs(iinc),'r');
    
    rarrow  = 0.4;
    xarrow  = rarrow*sin(para.gamr(iinc));
    yarrow  =-rarrow*cos(para.gamr(iinc));
    para.b_dib(ifig).harrow  = quiver(para.b_dib(ifig).hg,para.xs(iinc),para.zs(iinc),xarrow,yarrow,'r');
    set(para.b_dib(ifig).harrow,'MaxHeadSize',800);
    set(para.b_dib(ifig).harrow,'ShowArrowHead','off'); %bug matlab de mise a jour
    set(para.b_dib(ifig).harrow,'ShowArrowHead','on');
    %     end
else %if para.fuente==2
    %fuerza vectorial
    para.b_dib(ifig).xzs=plot(para.b_dib(ifig).hg,para.xs(iinc),para.zs(iinc),'r.');
    
    if para.dim==1
        if para.pol==1
            para.b_dib(ifig).xzs=plot(para.b_dib(ifig).hg,para.xs(iinc),para.zs(iinc),'ro');
        else
            rarrow  = 0.4;
            xarrow  = rarrow*sin(para.gamr(iinc));
            yarrow  =-rarrow*cos(para.gamr(iinc));
            para.b_dib(ifig).harrow  = quiver(para.b_dib(ifig).hg,para.xs(iinc),para.zs(iinc),xarrow,yarrow,'r');
            set(para.b_dib(ifig).harrow,'MaxHeadSize',800);
            set(para.b_dib(ifig).harrow,'ShowArrowHead','off'); %bug matlab de mise a jour
            set(para.b_dib(ifig).harrow,'ShowArrowHead','on');
        end
    else
        rarrow  = 0.4;
        xarrow  = rarrow*sin(para.gamr(iinc));
        yarrow  =-rarrow*cos(para.gamr(iinc));
        para.b_dib(ifig).harrow  = quiver(para.b_dib(ifig).hg,para.xs(iinc),para.zs(iinc),xarrow,yarrow,'r');
        set(para.b_dib(ifig).harrow,'MaxHeadSize',800,'ShowArrowHead','off');%bug matlab de mise a jour
        set(para.b_dib(ifig).harrow,'ShowArrowHead','on');
     end
end

%%%%%%%%%%
% titulo %
%%%%%%%%%%
titledibujo;
set(h201(ifig),'name',titletxt);

%%%%%%%%%%%%%%%%%%%%%
% senales en tiempo %
%%%%%%%%%%%%%%%%%%%%%
%parametros de normalizacion estandares
if para.rec.nrecz>1 && para.rec.nrecx==1
    dxr=para.rec.dzr;
    xr0=xr;
    xr=zr;
else
    dxr=para.rec.dxr;
end

it0=find(tps>=min(0,tps(1)),1,'first');
it1=find(tps>=min(10,tps(end)),1,'first');
% it1=find(tps>=min(80,tps(end)),1,'first');

%superposicion o no de las senales
tmpc    =get(bouton.couleur,'string');
coul    =get(bouton.couleur,'value');
if holdtr==2
    for j=1:fieldV(ifig).nc
        hold(para.b_dib(ifig).ht(j),'off');
    end
    couleur =tmpc{coul};
else
    hold on
    itmpc=mod(coul+iinc-1,length(tmpc));
    if itmpc==0;itmpc=length(tmpc);end
    couleur=tmpc{itmpc};
end

if fieldV(ifig).name(1)=='U'
    varplot=utc;
    ifig0=ifig;
elseif fieldV(ifig).name(1)=='S'
    varplot=stc;
    ifig0=1;
end
nrec = size(varplot,2);
wiggle = get(para.b_dib(ifig).wiggle,'value'); %disp(['wiggle=' num2str(wiggle)])
if para.dim < 3
    if para.rec.nrecy==1 && para.rec.nrecz==1
        %caso tipico de los receptores en superficie
        for j=1:fieldV(ifig).nc
            if para.b_dib(1).normalise==1
                mmax=max(max(abs(varplot(it0:it1,:,iinc,j+(ifig0-1)*fieldV(ifig).nc))));
                if mmax==0;mmax=1;end
            else
                mmax=1;
            end
            if (wiggle == 2);axes(para.b_dib(ifig).ht(j));cla;end
            for i=1:nrec
                if para.b_dib(1).normalise==1
                    if (wiggle == 2)
                    plot(para.b_dib(ifig).ht(j),tps,squeeze(varplot(:,i,iinc,j+(ifig0-1)*fieldV(ifig).nc))/mmax*dxr+xr(i),couleur);
                    else
                        plotwiggle(para.b_dib(ifig).ht(j),tps,squeeze(varplot(:,i,iinc,j+(ifig0-1)*fieldV(ifig).nc))/mmax*dxr+xr(i));
                    end
                else
                    if (wiggle == 2)
                    plot(para.b_dib(ifig).ht(j),tps,squeeze(varplot(:,i,iinc,j+(ifig0-1)*fieldV(ifig).nc)),couleur);
                    else
                        plotwiggle(para.b_dib(ifig).ht(j),tps,squeeze(varplot(:,i,iinc,j+(ifig0-1)*fieldV(ifig).nc)));
                    end
                end
                hold(para.b_dib(ifig).ht(j),'on');
            end
        end
    elseif para.rec.nrecy==1
        for j=1:size(varplot,4)%para.pol
            if para.b_dib(1).normalise==1
                mmax=max(max(abs(varplot(it0:it1,:,iinc,j+(ifig0-1)*fieldV(ifig).nc))));
                if mmax==0;mmax=1;end
            else
                mmax=1;
            end
            for k=1:para.rec.nrecz
                itmpc=mod(coul+iinc-1+k-1,length(tmpc));
                if itmpc==0;itmpc=length(tmpc);end
                couleur=tmpc{itmpc};
                for i=1:para.rec.nrecx
                    if (wiggle == 2)
      plot(para.b_dib(ifig).ht(j),tps,squeeze(varplot(:,i+(k-1)*para.rec.nrecx,iinc,j+(ifig0-1)*fieldV(ifig).nc))/mmax*dxr+xr(i),couleur)
                    else
plotwiggle(para.b_dib(ifig).ht(j),tps,squeeze(varplot(:,i+(k-1)*para.rec.nrecx,iinc,j+(ifig0-1)*fieldV(ifig).nc))/mmax*dxr+xr(i));
                    end
                    hold(para.b_dib(ifig).ht(j),'on');
                end
            end
        end
    end
else
    %3D
    if para.rec.nrecx > 2; mindis = abs(xr(2) - xr(1)); 
    elseif para.rec.nrecy >2; mindis = abs(yr(2) - yr(1)); 
    elseif para.rec.nrecz >2; mindis = abs(zr(2) - zr(1)); 
    else mindis = 1;
    end
    if para.b_dib(1).normalise==1
      rav=zeros(para.rec.nrecx*para.rec.nrecy*para.rec.nrecz,1);
      for iz=1:para.rec.nrecz
        for iy=1:para.rec.nrecy
          for ix=1:para.rec.nrecx
            ies = ix+(iy-1)*para.rec.nrecx+(iz-1)*para.rec.nrecx*para.rec.nrecy; 
  rav(ies) = ((xr(ix)-para.xs(iinc)).^2 + ...
              (yr(iy)-para.ys(iinc)).^2 + ...
              (zr(iz)-para.zs(iinc)).^2).^0.5;
          end
        end
      end
        varplot(it0:it1,rav<mindis*0.2,iinc,:) = NaN;
        mmax=max(max(max(abs(varplot(it0:it1,rav>=mindis*0.2,iinc,:)))));
        if mmax==0;mmax=1;end
    else
        mmax=1;
    end
    itmpc=mod(coul+iinc-1,length(tmpc));
    if itmpc==0;itmpc=length(tmpc);end
    couleur=tmpc{itmpc};
    
    for j=1:fieldV(ifig).nc%cada componente
        for iz=1:para.rec.nrecz
            for iy=1:para.rec.nrecy
                for ix=1:para.rec.nrecx
                  
%                   r = ((xr(ix)-para.xs(iinc)).^2 + ...
%                        (yr(iy)-para.ys(iinc)).^2 + ...
%                        (zr(iz)-para.zs(iinc)).^2).^0.5;
%                   if (r(1) < mindis)
%                     continue
%                   end
                    if (wiggle == 2)
      plot(para.b_dib(ifig).ht(j),tps,squeeze(varplot(:,ix+(iy-1)*para.rec.nrecx+(iz-1)*para.rec.nrecx*para.rec.nrecy,iinc,j))/mmax*dxr+xr(ix+(iy-1)*para.rec.nrecx+(iz-1)*para.rec.nrecx*para.rec.nrecy),couleur)
                    else
plotwiggle(para.b_dib(ifig).ht(j),tps,squeeze(varplot(:,ix+(iy-1)*para.rec.nrecx+(iz-1)*para.rec.nrecx*para.rec.nrecy,iinc,j))...
        /mmax*dxr+xr(ix+(iy-1)*para.rec.nrecx+(iz-1)*para.rec.nrecx*para.rec.nrecy));
                    end
                    hold(para.b_dib(ifig).ht(j),'on');
                end
            end
        end
    end
end
if para.b_dib(1).normalise==1
    ylimg=get(para.b_dib(ifig).hg,'xlim');
    ylimt=get(para.b_dib(ifig).ht(1),'ylim');
    ylimok=[min(ylimg(1),ylimt(1)) max(ylimg(2),ylimt(2))];
    set(para.b_dib(ifig).hg,'xlim',ylimok);
    for j=1:fieldV(ifig).nc
        set(para.b_dib(ifig).ht(j),'ylim',ylimok);
    end
end

if para.rec.nrecz>1 && para.rec.nrecx==1
    xr=xr0;
    for j=1:fieldV(ifig).nc%cada componente
        set(para.b_dib(ifig).ht(j),'ydir','reverse');
    end
end

% ligar los ejes
linkaxes(para.b_dib(ifig).ht(1:fieldV(ifig).nc),'xy');