function filmoscopio3D(para,utc,~,film,iinc0)
%filmoscopio Make a movie, receiv~ers on a grid and on the boundary.
%   The movie es made with data either from a grid of receivers or
%   from receivers on the boundary of a 2D surface.
%   Four styles of plot are available

filmeRange   = film.filmeRange;
filmStyle    = film.filmStyle;
filmeMecElem = film.filmeMecElem;
cmd_takingNames;
mecelemlist = film.strmecelemlist;
%if (para.rec.resatboundary); warpBoundaries = true; else warpBoundaries = false; end
if (filmStyle == 1); filmStyle = 3; warning ('Using filmSyle = grid with shadow');end

nf      = para.nf;
df      = para.fmax/(nf/2);     %paso en frecuencia
nfN     = nf/2+1; %Nyquist
zerospad= para.zeropad;
dt = (1/(df*2*(nfN+zerospad))*...
  (2*(nfN+zerospad)/(2*(nfN+zerospad)-2)));
tps     = 0:dt:1/df;
tps     = para.pulso.b+tps;

% tomamos las características de la vista previa en main:
[az,el]=view; ds = daspect; xl = xlim; yl = ylim; zl = zlim;

fig=figure('Position', [0, 0, 800, 800],'Color',[1 1 1]);
set(fig,'DoubleBuffer','on'); set(gcf,'Renderer','zbuffer')
[parentdir,~] = fileparts(pwd);
name=[parentdir,parentdir(1),'out',parentdir(1),'filmtmp','_000'];

nrecx=para.rec.nrecx; % inicializar variables
nrecy=para.rec.nrecy;
nrecz=para.rec.nrecz;
xr      = para.rec.xri+para.rec.dxr*(0:(nrecx-1));
yr      = para.rec.yri+para.rec.dyr*(0:(nrecy-1));
zr      = para.rec.zri+para.rec.dzr*(0:(nrecz-1));
mxWarp = abs(xr(end) - xr(1))/15;
if (mxWarp == 0)
  mxWarp = abs(yr(end) - yr(1))/15;
end
if (filmStyle ~= 4) % not a quiver plot
  [MX,MY,MZ] = meshgrid(xr,yr,zr);
end

while exist([name,'.avi'],'file')==2 %nombre del video
  compt=str2double(name(length(name)-2:length(name)));
  name=[name(1:length(name)-3),num2str(compt+1,'%.3i')];
end

mov = VideoWriter(name);%,'LosslessCompression','true');
mov.FrameRate = film.fps;
open(mov);
nam = char(mecelemlist(filmeMecElem));

for iinc=iinc0 %Para la incidenica actual
  cmd_iflist_loadthelast; %Cargar la subdivición de los elementos más reciente si hizo falta
  dibujo_conf_geo(para,gca);
  xlim('manual');ylim('manual');zlim('manual')
  hold on
  view(az,el);
  daspect(ds);
  if para.rec.nrecx > 2; mindis = abs(xr(2) - xr(1)); 
    elseif para.rec.nrecy >2; mindis = abs(yr(2) - yr(1)); 
    else mindis = abs(zr(2) - zr(1)); 
  end
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
  utc(filmeRange,rav<mindis*0.2,iinc,:) = NaN;
      
  m3max=max(max(max(max(utc(filmeRange,rav>=mindis*0.2,iinc,:)))))/mxWarp;
  u =utc(filmeRange,:,iinc,1)/m3max;
  v =utc(filmeRange,:,iinc,2)/m3max;
  w =utc(filmeRange,:,iinc,3)/m3max;
  
  minw = min(min(min(min(w))));
  r = (u.^2+v.^2+w.^2).^0.5;
%   m3maxS=max(max(max(max(abs(stc(filmeRange,:,iinc,1))))));
%   r = abs(stc(filmeRange,:,iinc,1)); % Para los esfuerzos
%   r = r / mean(mean(r,'omitnan'));
%    [ntf,nes]=size(r);
%   for it= 1:ntf
%     for ies = 1:nes
%   r(it,ies) = comprimir(r(it,ies),0.4);
%     end
%   end
  
  maxr = max(max(r));
  mx = 0.3*mean(maxr); %<-- Entre menor la constante, mayor la compresión
  txtiempo = squeeze(tps(filmeRange));
  
  iit=0;
  [ntf,nes]=size(u);
  for it= 1:ntf
    iit=iit+1;
      if ((filmStyle == 2) || (filmStyle == 3)); MXi=MX;MYi=MY;MZi=MZ;C=MX.*0;end
      if (filmStyle == 4); x=zeros(1,nes);y=x;z=x;mag=x;n=zeros(3,nes);end
      for iz=1:nrecz
        for iy=1:nrecy
          for ix=1:nrecx
            ies = ix+(iy-1)*nrecx+(iz-1)*nrecx*nrecy;
            if ((filmStyle == 2) || (filmStyle == 3))
              % mesh
              if (max(abs(sqrt(u(it,ies)^2+v(it,ies)^2+w(it,ies)^2)))==0)
                MZi(iy,ix,iz)= nan; %No traza el grid
              else
                MXi(iy,ix,iz)=MXi(iy,ix,iz)+u(it,ies);
                MYi(iy,ix,iz)=MYi(iy,ix,iz)+v(it,ies);
                MZi(iy,ix,iz)=MZi(iy,ix,iz)+w(it,ies);
                C(iy,ix,iz)  = -r(it,ies);
                %                 MZi(ix,iz)=MZi(ix,iz)-r(it,ies);
              end
            elseif (filmStyle == 4)
              % quiver
              x(1,ies) = xr(ix);
              y(1,ies) = yr(iy);
              z(1,ies) = zr(iz);
              n(1,ies) = u(it,ies);
              n(2,ies) = v(it,ies);
              n(3,ies) = w(it,ies);
              mag(ies) = comprimir(magnitud(n(1:3,ies)),mx);
            end
          end
        end
      end
      if (filmStyle == 2)
        tmph = mesh(squeeze(MXi),squeeze(MYi),squeeze(MZi),squeeze(C));
%         caxis([-40 1]); %todo blanco
        set(tmph,'FaceAlpha',0)
      elseif (filmStyle == 3)
        tmph = surf(squeeze(MXi),squeeze(MYi),squeeze(MZi),squeeze(C),'FaceColor','interp',...
          'FaceLighting','phong','AmbientStrength',.9,'DiffuseStrength',.8,...
          'SpecularStrength',.9,'SpecularExponent',25,...
          'BackFaceLighting','unlit');%,...
%           'EdgeColor','none','LineStyle','none');
%         caxis([-0.1*maxr 0.005*maxr]);
        alpha(tmph,0.6);
        set(tmph,'FaceAlpha',0.4);
      elseif (filmStyle == 4)
        iv = mag>0.6*mean(mag,'omitnan');
        tmph = quiver3(x(iv),y(iv),z(iv),n(1,iv),n(2,iv),n(3,iv),1.0*max(mag),'k');
      end
    colormap(gray);
    cm = colormap.^4;
    colormap(cm);
%     if it==1
%     axis tight
%     xl = xlim; yl = ylim; zl = zlim;
%     else
    xlim(xl); ylim(yl); zlim(zl); %garantizar que se vean todos los receptores
%     end
    zlim([minw,max(get(gca,'zlim'))]);
    if (filmStyle == 1)
      t_title = '|| U ||';
    else
      if strcmp(nam(1:1),'U')
        t_title = '|| U ||';
      else
        t_title = ['|| ' nam ' ||'];
      end
%       shading faceted
    end
    title(t_title)
    % Marca de tiempo
    tmpstr=[t_title '  ' num2str(txtiempo(iit)),' s'];
    %             h=annotation('textbox',[.35 .9 .3 .3]);
    %             set(h,'FitHeightToText','on',...
    %             'string',tmpstr,...
    %             'HorizontalAlignment','center','FontWeight','bold','LineStyle','none')
    
    h = uicontrol('Style','text','String',tmpstr,...
      'Units','normalized','Fontsize',12,...
      'Position',[0.35 0.15 0.3 0.035],'BackgroundColor',[1 1 1]);
    frame = getframe(gcf);
    writeVideo(mov,frame);
    pause(.1)
    delete(h)
    delete(tmph)
    %     hold off
  end
end
close(mov);
close(fig);
implay([name,'.avi'])
end

function [s] = comprimir(s,mx)
p = 15; %8
s =  log((1. + exp(p)*abs(s))) / (log(exp(p)+1.));
%disp(s)
s = s / mx;
end

function [r] = magnitud(s)
r = sqrt(sum((s(1:3)).^2));
end


