function dibujo_conf_geo(para,axe_conf_geo)
try
  bouton = evalin('base','bouton');
catch 
  bouton = [];
end
if isfield(bouton,'rafraichiEveryTime')
   refrescar = get(bouton.rafraichiEveryTime,'value');
else
   refrescar = 2; 
end
if refrescar == 1 % no refrescar
  return
end
if refrescar == 3 % refrescar sólo una vez
    if isfield(bouton,'rafraichiEveryTime')
      set(bouton.rafraichiEveryTime,'value',1);
    end
    refrescar = 2;
end
if refrescar == 2
axes(axe_conf_geo); cla; hold on
disp('refrescando imagen'); %drawnow limitrate; 
end
%-----------------------%
% parametros geometrico %
%-----------------------%
% reescritura de parametros de geometria simetrica
para.cont(1,1).za	= 0;
for m=1:para.nmed
  para.cont(m,2).a	= para.cont(m,1).a ;
  para.cont(m,2).xa	= para.cont(m,1).xa;
  para.cont(m,2).za	= para.cont(m,1).za;
  para.cont(m,1).th   = para.cont(m,1).th *pi/180;
  para.cont(m,2).th	= para.cont(m,1).th;
end

%% receptores
% if isfield(para,'cont1')
%     xr      = para.rec.xr;
%     zr      = para.rec.zr;
% else
para    = pos_rec(para);
xr      = para.rec.xr;
yr      = para.rec.yr;
zr      = para.rec.zr;
xranrec = [min(xr) max(xr)];
yranrec = [min(yr) max(yr)];
% zranrec = [min(zr) max(zr)];
% end
if isfield(bouton,'verReceptores')
    verRecepotores = get(bouton.verReceptores,'value');
else
    verRecepotores = true;
end
if verRecepotores
if para.dim == 4
  colFondo = 'k'; arrcol = ['k','b','r','g','y','m','c','w']; 
  if length(xr) > 10
    marcador = '.';
    markerTamano = 5;
  else
    marcador = '^';
    markerTamano = 15;
  end
%   if nargin == 2; marcador = 'none'; end
  
  for im = 1:para.nmed
    ip = find(para.rec.medio==im); 
    if isempty(ip); continue; end
    if im == 1
      plot3(axe_conf_geo,xr(ip),yr(ip),zr(ip),...
        ['.' colFondo],'MarkerSize',10,'Marker',marcador,'MarkerSize',markerTamano)
    else
      plot3(axe_conf_geo,xr(ip),yr(ip),zr(ip),...
      ['.' arrcol(para.cont(para.rec.medio(ip),1).piece{1}.ColorIndex+1)],...
      'MarkerSize',10,'Marker',marcador,'MarkerSize',markerTamano)
    end
  end
  view(3);
elseif para.dim == 3
  plot3(axe_conf_geo,xr,yr,zr,'.b','MarkerSize',10);view(3)
else
  plot(axe_conf_geo,xr,zr,'.b'); view(2)
end
end

%% fronteras de la irregularidad
% if nargin == 4; EA = 0.2; elseif nargin == 2; EA = 0.; else EA = 0.5;end
EA = 0.2;
for m=2:para.nmed
  if para.dim == 4 % 3Dgeneral
    view(3)
    if isfield(bouton,'verGeometria')
      verGeometria = get(bouton.verGeometria,'value'); 
    else
      verGeometria = true;
    end
    if verGeometria
    for p = 1:para.cont(m,1).NumPieces % para cada pieza del contorno
      if (size(para.cont(m,1).piece{p}.fileName,2)>1) % se cargó algo válido
        if size(para.cont(m,1).piece{p}.geoFileData,2)>0 % y hay datos
          % de archivo STL
          thisGF = para.cont(m,1).piece{p}.geoFileData;
%           if nargin > 2;
          if isfield(bouton,'verNormales')
            if ~get(bouton.verNormales,'value'); 
              thisGF = rmfield(thisGF,'N'); 
            end
          else
              thisGF = rmfield(thisGF,'N'); 
          end
          plotSTL(axe_conf_geo,thisGF,...
                               para.cont(m,1).piece{p}.ColorIndex,EA);
          clear thisGF
        end
      end
    end
    end
  elseif para.dim == 3 % 3Daxisim
    % axisimétrico
    if (para.geo(m) ~= 5) % todo menos STL
      %monte
      [xm,zm] = visu_courbe(para.geo(m),para.cont(m,1),1000);
      plot3(axe_conf_geo,xm,xm.*0,zm,'k');
      
      %valle
      [xv,zv] = visu_courbe(para.geo(m),para.cont(m,2),1000);
      plot3(axe_conf_geo,xv,xv.*0,zv,'k');
    end
  else % 2D ou 2.5D
    if (para.geo(m) ~= 5) % todo menos STL
      %monte
      [xm,zm] = visu_courbe(para.geo(m),para.cont(m,1),1000);
      plot(axe_conf_geo,xm,zm,'k');
      
      %valle
      [xv,zv] = visu_courbe(para.geo(m),para.cont(m,2),1000);
      plot(axe_conf_geo,xv,zv,'k');
    end
  end
end

%% fuente
if isfield(bouton,'inc')
iinc=get(bouton.inc,'value');
else
  if isfield(bouton,'iinc')
    iinc = bouton.iinc;
  else
    iinc = 1;
  end
end

if isfield(bouton,'gam')
  gam = str2double(get(bouton.gam,'string')); 
  gam=gam*pi/180;
  phi = str2double(get(bouton.phi,'string')); 
  phi=phi*pi/180;
else
  gam = para.gam(iinc)*pi/180;
  phi = para.phi(iinc)*pi/180;
end

%     hcent   = plot(axe_conf_geo,para.xs(iinc),para.zs(iinc),'r.');
if para.dim >= 3
  plot3(axe_conf_geo,para.xs,para.ys,para.zs,'r.','MarkerSize',20);
else
  plot(axe_conf_geo,para.xs,para.zs,'r.','MarkerSize',20);
end

if para.fuente==1
  %     for iinc=1:min(para.ninc,5)
  xarrow  =-.25*cos(para.gam(iinc)*pi/180);
  yarrow  =-.25*sin(para.gam(iinc)*pi/180);
  x       = [xarrow -xarrow];
  y       = [yarrow -yarrow];
  if para.dim >= 3
    plot3(axe_conf_geo,x+para.xs(iinc),x.*0,y+para.zs(iinc),'r');
  else
    plot(axe_conf_geo,x+para.xs(iinc),y+para.zs(iinc),'r');
  end
  
  rarrow  = 0.4;
  xarrow  = rarrow*sin(para.gam(iinc)*pi/180);
  yarrow  =-rarrow*cos(para.gam(iinc)*pi/180);
  if para.dim >= 3
    harrow  = quiver3(axe_conf_geo,para.xs(iinc),para.ys(iinc),para.zs(iinc),xarrow,0,yarrow,'r');
  else
    harrow  = quiver(axe_conf_geo,para.xs(iinc),para.zs(iinc),xarrow,yarrow,'r');
  end
  set(harrow,'MaxHeadSize',800);
  set(harrow,'ShowArrowHead','off'); %bug matlab de mise a jour
  set(harrow,'ShowArrowHead','on');
  %     end
elseif  para.fuente==2
  if para.dim==1
    if para.pol==1
      plot(axe_conf_geo,para.xs,para.zs,'ro');
    elseif para.pol==2
      rarrow  = 0.4;
      for iinc=1:min(para.ninc,5)
        xarrow  = rarrow*sin(para.gam(iinc)*pi/180);
        yarrow  =-rarrow*cos(para.gam(iinc)*pi/180);
        if para.dim >= 3
          harrow  = quiver3(axe_conf_geo,para.xs(iinc),para.ys(iinc),para.zs(iinc),xarrow,0,yarrow,'r');
        else
          harrow  = quiver(axe_conf_geo,para.xs(iinc),para.zs(iinc),xarrow,yarrow,'r');
        end
        set(harrow,'MaxHeadSize',800);
        set(harrow,'ShowArrowHead','off'); %bug matlab de mise a jour
        set(harrow,'ShowArrowHead','on');
      end
    end
  else
%     rarrow  = 0.4;
    extent = get(gca,'xlim');
    rarrow = 0.2*(max(extent)-min(extent));
%     for iinc=1:min(para.ninc,5)
      if para.dim >= 3
 fij     = [(sin(gam).*cos(phi)).' ...
            (sin(gam).*sin(phi)).' ...
            -cos(gam).'];
  xarrow  = rarrow*fij(1);
  yarrow  = rarrow*fij(2);
  zarrow  = rarrow*fij(3);
        
        
        harrow  = quiver3(axe_conf_geo,para.xs(iinc),para.ys(iinc),para.zs(iinc),xarrow,yarrow,zarrow,'r');
      else
      xarrow  = rarrow*sin(para.gam(iinc)*pi/180);
      yarrow  =-rarrow*cos(para.gam(iinc)*pi/180);
        harrow  = quiver(axe_conf_geo,para.xs(iinc),para.zs(iinc),xarrow,yarrow,'r');
      end
      set(harrow,'MaxHeadSize',800);
      set(harrow,'ShowArrowHead','off'); %bug matlab de mise a jour
      set(harrow,'ShowArrowHead','on');
%     end
  end
end
if (para.dim >= 3)
  set(axe_conf_geo,'zdir','reverse','dataaspectratio',[1 1 1]);
  axes(axe_conf_geo)
  xlabel('x')
  ylabel('y')
  zlabel('z')
else
  set(axe_conf_geo,'ydir','reverse','dataaspectratio',[1 1 1]);
  axes(axe_conf_geo)
  xlabel('x')
  ylabel('z')
  zlabel('')
end

if para.dim ==4
  set(axe_conf_geo,'Projection','perspective','Box','off')
  grid on
  light('Position',[2*min(get(axe_conf_geo,'xlim')) ...
    2*max(get(axe_conf_geo,'ylim')) ...
    2*max(get(axe_conf_geo,'zlim'))],'Style','infinite'); %ambient
  axis tight
  
  [limx] = xlim; dlimX=limx(2)-limx(1); 
  [limy] = ylim; dlimY=limy(2)-limy(1); 
  [limz] = zlim; dlimZ=limz(2)-limz(1);
  dlimH = min(dlimX,dlimY);
  if dlimH < 0.3 *dlimZ; dlimH = 0.3 *dlimZ; end
  xlim([min(limx(1),mean(limx)-dlimH*0.5) max(limx(2),mean(limx)+dlimH*0.5)])
  ylim([min(limy(1),mean(limy)-dlimH*0.5) max(limy(2),mean(limy)+dlimH*0.5)])
  
  
%   if ~(nargin > 2); zlim([ -1 max(get(axe_conf_geo,'zlim'))]); end
  set(axe_conf_geo, 'XColor', 'r');set(axe_conf_geo, 'YColor', 'g');set(axe_conf_geo, 'ZColor', 'b')
else
  set(axe_conf_geo,'Projection','orthographic','Box','on');
  grid off
  xlabel('X');ylabel('Z');zlabel('');
  set(axe_conf_geo,'dataaspectratio',[1 1 1])
  axis image
  set(axe_conf_geo, 'XColor', 'k');set(axe_conf_geo, 'YColor', 'k');set(axe_conf_geo, 'ZColor', 'k')
end

%% background
 Xl = get(axe_conf_geo,'xlim'); Xl(1) = min(xranrec(1),Xl(1)); Xl(2) = max(xranrec(2),Xl(2));
 Yl = get(axe_conf_geo,'ylim'); Yl(1) = min(yranrec(1),Yl(1)); Yl(2) = max(yranrec(2),Yl(2));
%  Zl = get(axe_conf_geo,'zlim'); Zl(1) = min(zranrec(1),Zl(1)); Zl(2) = max(zranrec(2),Zl(2));
if para.geo(1)==3 %Medio estratificado DWN
  % construir vector xp
  xmin=0; xmax=0;
  for i=2:para.nmed
    xmin=min(xmin,para.cont(i,1).xa);
    xmax=max(xmax,para.cont(i,1).xa+2*para.cont(i,1).a);
  end
  xp = linspace(xmin-2,xmax+2,100);
  
  % superficie libre
  zp = 0*xp;
  if para.dim == 3 % 3Daxisim
    plot3(axe_conf_geo,xp,xp.*0,zp,'g');
  elseif para.dim == 4 % 3Dgen
    if ~isempty(xp)
        plot3(axe_conf_geo,Xl,[Yl(1) Yl(1)],zp(1:2),'g');
        plot3(axe_conf_geo,Xl,[Yl(2) Yl(2)],zp(1:2),'g');
        plot3(axe_conf_geo,[Xl(1) Xl(1)],Yl,zp(1:2),'g');
        plot3(axe_conf_geo,[Xl(2) Xl(2)],Yl,zp(1:2),'g');
    end
  else %2D
    plot(axe_conf_geo,xp,zp,'g');
  end
  
  % estratos
  for i=1:para.nsubmed-1
    xp = linspace(xmin-2,xmax+2,100);
    zp = zp+para.reg(1).sub(i).h;
    if para.dim == 3 % 3Daxisim
      plot3(axe_conf_geo,xp,xp.*0,zp,'g');
    elseif para.dim == 4 % 3Dgen
      if ~isempty(xp)
        plot3(axe_conf_geo,Xl,[Yl(1) Yl(1)],zp(1:2),'g');
        plot3(axe_conf_geo,Xl,[Yl(2) Yl(2)],zp(1:2),'g');
        plot3(axe_conf_geo,[Xl(1) Xl(1)],Yl,zp(1:2),'g');
        plot3(axe_conf_geo,[Xl(2) Xl(2)],Yl,zp(1:2),'g');
      end
    else % 2D
      plot(axe_conf_geo,xp,zp,'g');
    end
  end
else %semiespacio o espacio completo
  [xp,zp]=visu_courbe_m1(para.geo(1),para.cont(1,1),1000);
  if para.dim >= 3
    plot3(axe_conf_geo,xp,xp.*0,zp,'g');
    if ~isempty(xp)
%     Xl = get(axe_conf_geo,'xlim');
%     Yl = get(axe_conf_geo,'ylim');
    plot3(axe_conf_geo,Xl,[Yl(1) Yl(1)],zp(1:2),'g');
    plot3(axe_conf_geo,Xl,[Yl(2) Yl(2)],zp(1:2),'g');
    plot3(axe_conf_geo,[Xl(1) Xl(1)],Yl,zp(1:2),'g');
    plot3(axe_conf_geo,[Xl(2) Xl(2)],Yl,zp(1:2),'g');
    end
  else
    plot(axe_conf_geo,xp,zp,'g');
  end
end

hold off
drawnow update
% end
%% Graficar la estratificación
% if nargin == 4  %Medio estratificado DWN
if isfield(bouton,'axe_estrDWN')
  axe_estrDWN = bouton.axe_estrDWN;
else
  figure(32947202)
  set(gcf,'Name','Estratificacion','numberTitle','off');
  axe_estrDWN = gca;
end
if para.geo(1)==3 && ishandle(axe_estrDWN)
  axes(axe_estrDWN); cla
%   hfigProf = figure(321);
%   set(hfigProf,...
%     'Units','normalized','name','Estratificación',...%'color',[1 1 1], ...
%     'numberTitle','off','DockControls','off');
%     'position',[0.005 0.01 .51 .9], ...
     ParaRegSub = para.reg(1).sub(:);
     if para.dim == 1
       ShouldUseAlf = false;
     else
       ShouldUseAlf = true;
     end
  plotProfile(ParaRegSub,ShouldUseAlf);
  
  clear ParaRegSub ShouldUseAlf
  axes(axe_conf_geo)
end

end
