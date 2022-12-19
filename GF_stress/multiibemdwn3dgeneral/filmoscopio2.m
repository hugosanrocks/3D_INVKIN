function filmoscopio2(para,RESULT,iinc0)
utc = RESULT.utc;
stc = RESULT.stc;
film = para.film;
%filmoscopio Make a movie, receivers on a grid and on the boundary.
%   The movie es made with data either from a grid of receivers,
%   from receivers on the boundary (The boundary appears deformed)
%   or both.
if (para.dim >= 3) 
  filmoscopio3D(para,utc,stc,film,iinc0);return
end
filmeRange = film.filmeRange;
filmStyle    = film.filmStyle;
filmeMecElem = film.filmeMecElem;
cmd_takingNames;
mecelemlist = film.strmecelemlist;
malla = false; warpBoundaries = false;
if (para.rec.resatboundary); warpBoundaries = true;end
if (para.recpos == 2); malla = true;end
if (filmStyle == 4); malla = false;end
if ((~malla && ~warpBoundaries) || ~exist('utc','var'))
  msgbox('No results to plot in a snapshot'); return;end
if para.pol==1; filmStyle = 1; warpBoundaries = false;end % SH

nf      = para.nf;
df      = para.fmax/(nf/2);     %paso en frecuencia
Fq      = (0:nf/2)*df;
nfN     = nf/2+1; %Nyquist
zerospad= para.zeropad;
tps     = 0:(1/(df*2*(nfN+zerospad))*...
  (2*(nfN+zerospad)/(2*(nfN+zerospad)-2))):1/df;
tps     = para.pulso.b+tps;

[n1,n2,n3,n4]=size(utc);
nt=n1;
listinc={'P','S','R'};
fig=figure('Position', [0, 0, 800, 800],'Color',[1 1 1]);
set(fig,'DoubleBuffer','on');
[parentdir,~] = fileparts(pwd);
name=[parentdir,parentdir(1),'out',parentdir(1),'filmtmp','_000'];
% cd ..; cd out
% thisdir = pwd;
% cd ..; cd multi-dwn-ibem.matlab/
% name=[thisdir,thisdir(1),'filmtmp','_000'];
cont1 = para.cont1;

nresAtBoundary = para.rec.nresAtBoundary;
nrecx=para.rec.nrecx-nresAtBoundary;

nrecy=para.rec.nrecy;
nrecz=para.rec.nrecz;

xr      = para.rec.xri+para.rec.dxr*(0:(nrecx-1));
zr      = para.rec.zri+para.rec.dzr*(0:(nrecz-1));
Ut      = zeros(nrecx,nrecz);
mxWarp = abs(xr(end) - xr(1))/15;
if (mxWarp == 0)
  mxWarp = abs(max(max(cont1(:).vec.xc)) - min(min(cont1(:).vec.xc)))/15;
end
if ((filmStyle == 2) || (filmStyle == 3)) % estilo tipo mesh y mesh con sombreado
  [MX,MY] = meshgrid(para.rec.xri:para.rec.dxr:para.rec.xri+para.rec.dxr*(nrecx-1),...
    para.rec.zri:para.rec.dzr:para.rec.zri+para.rec.dzr*(nrecz-1));
  MZ = zeros(nrecx,nrecz);
end
while exist([name,'.avi'],'file')==2
  compt=str2double(name(length(name)-2:length(name)));
  name=[name(1:length(name)-3),num2str(compt+1,'%.3i')];
end

mov = VideoWriter(name);%,'LosslessCompression','true');
open(mov);
nam = char(mecelemlist(filmeMecElem));
if strcmp(nam(1:1),'s')
  filmeMecElem = filmeMecElem - size(utc,4);
  filmeMecElem = min(filmeMecElem,3); %filmeMecElem=1;
end

for iinc=iinc0
  hold off
  % u,v : señal en el rango de tiempo y todos los receptores
  if para.pol==1 % SH
    u =sign(utc(filmeRange,:,iinc,1)).*log(abs(utc(filmeRange,:,iinc,1)).*1e3+1);
    m3max=max(max(max(max(abs(u)))));
    u = u /m3max;
  else % P-SV
    if (filmStyle == 1)
      u =log((sqrt(utc(filmeRange,:,iinc,1).^2+utc(filmeRange,:,iinc,2).^2)).*1e9+1);
      m3max=max(max(max(max(abs(u)))));
      u = u /m3max;
    elseif (filmStyle == 4 && ~strcmp(nam(1:1),'U')) % mec elem & esfuerzo
        m3max=max(max(max(max(stc(filmeRange,:,iinc,:)))))/mxWarp;
        u = stc(filmeRange,:,iinc,filmeMecElem); % uok = u;
        u = u/m3max;
        r = u;
        maxr = max(max(r));
    else %filmStyle = 2 o 3 (malla)
      m3max=max(max(max(max(utc(filmeRange,:,iinc,:)))))/mxWarp;
      u =utc(filmeRange,:,iinc,1)/m3max;
      v =utc(filmeRange,:,iinc,2)/m3max;
      r = (u.^2+v.^2).^0.5;
      maxr = max(max(r));
    end
  end
  txtiempo = squeeze(tps(filmeRange));
  
  iit=0;
  [ntf,nes]=size(u);
  hbou = zeros(1,ntf);
  for it= 1:ntf
    iit=iit+1;
    if (malla)
      if ((filmStyle == 2) || (filmStyle == 3)); MXi=MX;MYi=MY;MZi=MZ;end
      for iz=1:nrecz
        for iy=1:nrecy
          for ix=1:nrecx
            ies = ix+(iy-1)*nrecx+(iz-1)*nrecx*nrecy;
            
            if (filmStyle == 1)
              if max(abs(u(:,ies)))==0
                Ut(ix,iz)   = nan;
              else
                Ut(ix,iz)   = u(it,ies);
              end
            elseif ((filmStyle == 2) || (filmStyle == 3))
              if (max(abs(sqrt(u(it,ies)^2+v(it,ies)^2)))==0)
                MZi(ix,iz)= nan; %No traza el grid
              else
                MXi(iz,ix)=MXi(iz,ix)+u(it,ies);
                MYi(iz,ix)=MYi(iz,ix)+v(it,ies);
                MZi(ix,iz)=MZi(ix,iz)-r(it,ies);
              end
            end
          end
        end
        
        if (filmStyle == 1)
          Ut(isnan(Ut))=0;
          tmph = surf(xr,zr,Ut.');
          caxis([-1 1]);
        elseif (filmStyle == 2)
          tmph = surf(MXi,MYi,MZi');
          % caxis([-1 1]); % todo gris
          caxis([-40 1]); %todo blanco
        elseif (filmStyle == 3)
          tmph = surf(MXi,MYi,MZi');
          caxis([-1*maxr 0.05*maxr]); % fondo blanco, color en los frentes de onda
        end
      end
    end
    colormap(gray);
    
    hold on;
    plot_contour;
    
    if (filmStyle == 1)
      t_title = '|| U ||';
    else
      if strcmp(nam(1:1),'U')
        t_title = '|| U ||';
      else
        t_title = ['|| ' nam ' ||'];
      end
      shading faceted
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
    frame = getframe;
    writeVideo(mov,frame);
    pause(.1)
    delete(h)
  end
end
close(mov);
close(fig);
implay([name,'.avi'])
end

