%% Parametros del IBEM v0.2
%%%%%%%%%%%%%%%%%%%%%%

clear;
%%
%calculo normal =0, calculo con equipartition =1
para.meth_PS=0;
para.zeropad=5*512;

% UI defaults
% para.smallscreen = true;
set(0,'defaultUicontrolFontSize',12)
set(0,'defaultUicontrolFontName','Times')
c_1 = [0.9 0.9 0.9];   
c_2 = [0.85 0.85 0.9]; %[0.3 0.8 1];
c_3 = [0.8 0.8 0.8];   %[1 .1 .1];
c_4 = [0.6 0.6 0.6];
%par uso con ventanas
para.siDesktop=true;
%sino para uso con matlab -nodesktop, hay que escribir: main_nodesktop

% format short g
format short e

if exist('chemin','file')==2
    fid = fopen('chemin');
    t = textscan(fid,'%s');
    t = t{1};
    fclose(fid);clear fid;
    dossierok=0;
    for i=1:length(t)
        para.nomrep=t{i};
        if exist(para.nomrep,'dir')==7
            dossierok=1;
            break
        end
    end
    if dossierok==0
        para.nomrep=chercherep;
        fid = fopen('chemin','a');
        fprintf(fid,'\n%s',para.nomrep);
        fclose(fid);
    end
    clear dossierok
else
  para.nomrep=pwd;
  nomrep=para.nomrep;
end

h0fig1=figure(100); clf
uicg = uipanel('parent',gcf,'Title','Configuracion geometrica','HighlightColor',[1 1 1],...
    'BackgroundColor',c_1,'Position',[.005 .005 .99 .49]);
para.bar = uicontrol('parent',uicg,'Style','text',...
    'BackgroundColor',[.9 .9 .9],'Units','normalized','Fontsize',14,...
    'position',[0 0 1 .06],'string','');
cmd_fig_1=['if(ishandle(para.bar));',...
'set(para.bar,''BackgroundColor'',[.9 .9 .9],''string'','''');',...
'end;rafraichi;'];
% ';%['pause(.1);if(ishandle(para.bar));',...
% 'set(para.bar,''BackgroundColor'',[.9 .9 .9],''string'','''');',...
% 'end;rafraichi;'];

set(h0fig1,...
    'Units','normalized','name','parametros IBEM','color',[0.2 .1 .9], ...
    'position',[0.005 0.01 .51 .9],'WindowButtonDownFcn',cmd_fig_1, ...
    'numberTitle','off','DockControls','off','MenuBar','none',...
    'CloseRequestFcn','cmd_cerrarCallback; if closeMe; closereq; end');

h0rtDad  = uitabgroup('Parent',h0fig1,...
    'Position',[.005 .5 0.99 0.495],'TabLocation','top');
h0rt(1)= uitab('parent',h0rtDad,'title','Modelo');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% propriedades de los materiales isotropos %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uipm = uipanel('parent',h0rt(1),'Title','Propriedades de los medios',...
    'BackgroundColor',c_1,'Position',[.005 .6 .99 .4]);

% dimension del problema
para.dim    = 1;
bouton.dim=uicontrol('parent',uipm,'Style','popupmenu',...
    'BackgroundColor',c_2,'Units','normalized',...
    'position',[0.01 0.75 .15 .2],'string',{'2D','2.5D','3D axisimetric (z)','3D irregular geometry'},...
    'Callback','cmd_fixboutonGeo;para.dim=get(bouton.dim,''value'');cmd_dim;');

% numero de medios
para.nmed   = 2; %para.reg=zeros(para.nmed,1);
uicontrol('parent',uipm,'Style','text','BackgroundColor',c_3,...
    'Units','normalized','position',[0.01 0.5 .18 .15],...
    'string','Numero de');
uicontrol('parent',uipm,'Style','text','BackgroundColor',c_3,...
    'Units','normalized','position',[0.01 0.25 .1 .20],...
    'string','medios:');
uicontrol('parent',uipm,'Style','text','BackgroundColor',c_3,...
    'Units','normalized','position',[0.20 0.5 .1 .15],...
    'string','selección');
bouton.nmed=uicontrol('parent',uipm,'Style','edit',...
    'BackgroundColor',c_2,'Units','normalized',...
    'position',[0.11 0.25 .08 .2],'string',para.nmed,...
    'Callback','cmd_nmed;');

% se trata el medio #
strmed      = 1:para.nmed;
bouton.med=uicontrol('parent',uipm,'Style','popupmenu',...
    'BackgroundColor',c_2,'Units','normalized',...
    'position',[0.20 0.25 .1 .2],'string',strmed,...
    'Callback','cmd_med;');

% se trata el sub medio #
strsubmed    = 1:2;
info.submed  = uicontrol('parent',uipm,'Style','text',...
    'BackgroundColor',c_3,'Units','normalized',...
    'position',[0.01 0.01 .1 .20],'string','estratos:');
bouton.submed= uicontrol('parent',uipm,'Style','popupmenu',...
    'BackgroundColor',c_2,'Units','normalized',...
    'position',[0.20 0.0 .1 .2],'string',strsubmed,...
    'Callback','cmd_submed;');
% numero de sub-medios
para.nsubmed   = 2;
bouton.nsubmed = uicontrol('parent',uipm,'Style','edit',...
    'BackgroundColor',c_2,'Units','normalized',...
    'position',[0.11 0.0 .08 .2],...
    'string',para.nsubmed,'Callback','cmd_nsubmed;');
  
% rho
uicontrol('parent',uipm,'Style','text','BackgroundColor',c_3,...
    'Units','normalized','position',[0.31 0.55 .1 .4],...
    'string','r','fontname','symbol');
for i=1:para.nmed
    para.reg(i).rho  = 1;
end
bouton.rho=uicontrol('parent',uipm,'Style','edit',...
    'BackgroundColor',c_2,'Units','normalized',...
    'position',[0.31 0.15 .1 .4],'string',para.reg(1).rho,...
    'Callback','para.reg(get(bouton.med,''value'')).rho=str2num(get(bouton.rho,''string''));cmd_lambda_mu;');

% velocidad P
info.alpha=uicontrol('parent',uipm,'Style','text',...
    'BackgroundColor',c_3,'Units','normalized',...
    'position',[0.42 0.75 .1 .2],'string','a','fontname','symbol');
for i=1:para.nmed
    para.reg(i).alpha= 4;
end
bouton.alpha=uicontrol('parent',uipm,'Style','edit','BackgroundColor',c_2,...
    'Units','normalized','position',[0.42 0.55 .1 .2],...
    'string',para.reg(1).alpha,'Callback','para.reg(get(bouton.med,''value'')).alpha=str2num(get(bouton.alpha,''string''));cmd_lambda_mu;');

% velocidad S
uicontrol('parent',uipm,'Style','text','BackgroundColor',c_3,...
    'Units','normalized','position',[0.53 0.75 .1 .2],...
    'string','b','fontname','symbol');
for i=1:para.nmed
    para.reg(i).bet= 1;
end
bouton.bet=uicontrol('parent',uipm,'Style','edit',...
    'BackgroundColor',c_2,'Units','normalized','position',[0.53 0.55 .1 .2],...
    'string',para.reg(1).bet,'Callback','para.reg(get(bouton.med,''value'')).bet=str2num(get(bouton.bet,''string''));cmd_lambda_mu;');

% lambda
for i=1:para.nmed
    para.reg(i).lambda	= para.reg(i).rho*(para.reg(i).alpha^2-2*para.reg(i).bet^2);
end
info.lambda     = uicontrol('parent',uipm,'Style','text',...
    'BackgroundColor',c_3,'Units','normalized','position',[0.42 0.25 .1 .2],...
    'string','l','fontname','symbol');
bouton.lambda   = uicontrol('parent',uipm,'Style','text',...
    'BackgroundColor',c_2,'Units','normalized','position',[0.42 0.05 .1 .2]);

% mu
for i=1:para.nmed
    para.reg(i).mu      = para.reg(i).rho*para.reg(i).bet^2;
end
                 uicontrol('parent',uipm,'Style','text',...
    'BackgroundColor',c_3,'Units','normalized','position',[0.53 0.25 .1 .2],...
    'string','m','fontname','symbol');
bouton.mu       = uicontrol('parent',uipm,'Style','text',...
    'BackgroundColor',c_2,'Units','normalized','position',[0.53 0.05 .1 .2]);

% attenuation model
for i=1:para.nmed
para.reg(i).tipoatts=1;
end
lstatt={'Q','KV','no at'};
              uicontrol('parent',uipm,'Style','text',...
    'BackgroundColor',c_3,'Units','normalized','position',[0.64 0.05 .1 .9],...
    'string','modelo de atenuacion');
bouton.lstatts=uicontrol('parent',uipm,'Style','popupmenu',...
    'BackgroundColor',c_2,'Units','normalized','position',[0.64 0.05 .1 .4],...
    'string',lstatt,'Callback','para.reg(get(bouton.med,''value'')).tipoatts=get(bouton.lstatts,''value'');cmd_att;');

% Q model
info.att=uicontrol('parent',uipm,'Style','text','string','Q factor',...
    'BackgroundColor',c_3,'Units','normalized','position',[0.75 0.55 .1 .4]);
for i=1:para.nmed
    para.reg(i).qd= 1000;
end
bouton.Q=uicontrol('parent',uipm,'Style','edit','BackgroundColor',c_2,...
    'Units','normalized','position',[0.75 0.15 .1 .4],...
    'string',para.reg(1).qd ,'Callback','para.reg(get(bouton.med,''value'')).qd =str2num(get(bouton.Q  ,''string''));');

% espesor sub estrato
para.reg(1).sub(1).h=1;
info.subh    = uicontrol('parent',uipm,'Style','text',...
    'BackgroundColor',c_3,'Units','normalized',...
    'position',[0.86 0.55 .1 .4],'string','h');
bouton.subh  = uicontrol('parent',uipm,'Style','edit',...
    'BackgroundColor',c_2,'Units','normalized',...
    'position',[0.86 0.05 .1 .4],'string',para.reg(1).sub(1).h,...
    'Callback','para.reg(1).sub(get(bouton.submed,''value'')).h=str2num(get(bouton.subh,''string''));');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% geometria de los medios 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
uigm = uipanel('parent',h0rt(1),'Title','Geometria del medio #',...
    'HighlightColor' , [1 1 1],...
    'BackgroundColor',c_1,'Position',[.005 .0 .99 .6]);

dx  =1/12;
dxi =9/10*dx;
x0  =0.005;
info.parpos =uicontrol('parent',uigm,'Style','text',...
    'BackgroundColor',[.7 .7 .7]  ,'Units','normalized',...
    'position',[x0+0.95*dx 0.05 3*dx+1.1*dxi 0.95],'string',...
    'parametros de posicion del contorno');
info.rug   =uicontrol('parent',uigm,'Style','text',...
    'BackgroundColor',[.7 .7 .7]  ,'Units','normalized',...
    'position',[x0+8.95*dx 0.05 2*dx+1.1*dxi .95],'string','rugosidad');
  
bouton.geoFilePanel = uipanel('parent',uigm,'Title','',...
    'BackgroundColor',c_1,'Position',[x0+0.95*dx 0.25 0.91 0.75]);
bouton.gfPreview = axes('parent',bouton.geoFilePanel,'Units','normalized',...
    'OuterPosition',[0.4 0 0.6 1],'Layer','top',...
    'HitTest','off','zdir','reverse','dataaspectratio',[1 1 1]);
uicontrol('parent',bouton.geoFilePanel,'Style','text','BackgroundColor',c_3,'Units','normalized',...
    'position',[0.0 0.58 0.4 0.35],'string','Numero de Piezas          pieza #');
bouton.gfNumPieces = uicontrol('parent',bouton.geoFilePanel,'Style','edit','BackgroundColor',c_2,'Units','normalized',...
    'position',[0.0 0.60 0.2 0.15],'string',1,...
    'Callback','cmd_dimension_gfNumPieces;');
info.ThisPiece = 1;
para.cont(1,1).NumPieces = 1;
para.cont(1,1).piece=cell(1); 
para.cont(1,1).piece{1}.fileName = '';
para.cont(1,1).piece{1}.kind = 0;
para.cont(2,1).piece=cell(1); 
para.cont(2,1).NumPieces = 1;
para.cont(2,1).piece{1}.fileName = '';
para.cont(2,1).piece{1}.kind=1;
para.cont(2,1).piece{1}.continuosTo=1;
para.cont(2,1).piece{1}.ColorIndex=1;
bouton.gfThisPiece = uicontrol('parent',bouton.geoFilePanel,'Style','popupmenu','BackgroundColor',c_3,'Units','normalized',...
    'position',[0.2 0.50 0.2 0.25],'string',{'1'},...
    'Callback','cmd_geoThisPiece;');
bouton.gfThisPieceKind = uicontrol('parent',bouton.geoFilePanel,'Style','popupmenu','BackgroundColor',c_3,'Units','normalized',...
    'position',[0.0 0.25 0.4 0.25],'string',{'Free surface for the medium #:','Continuity with the medium #:','Just to close the polyhedron'},...  
    'Callback','cmd_geoKind;');
bouton.gfThisPieceContinuosTo = uicontrol('parent',bouton.geoFilePanel,'Style','edit','BackgroundColor',c_2,'Units','normalized',...
    'position',[0.0 0.1 0.2 0.15],'string','1','visible','on',...
    'Callback',['val = str2double(get(bouton.gfThisPieceContinuosTo,''string''));',...
                ' if(val>para.nmed); warning('' The value should be <= nmed ''); else ',...
                'para.cont(get(bouton.med,''value''),1).piece{info.ThisPiece}.continuosTo = val; clear val; end']);
bouton.gfThisPieceColor = uicontrol('parent',bouton.geoFilePanel,'Style','popupmenu','BackgroundColor',c_3,'Units','normalized',...
    'position',[0.2 0.0 0.2 0.25],'string',{'b','r','g','y','m','c','none'},...
    'Callback','para.cont(get(bouton.med,''value''),1).piece{info.ThisPiece}.ColorIndex = get(bouton.gfThisPieceColor,''value'');');
bouton.geoFileSelect =uicontrol('parent',uigm,'Style','pushbutton','BackgroundColor',[1 1 1],...
    'Units','normalized','position',[x0+4.3*dxi 0.05 8.7*dxi .2],...
    'string','Load file','Callback','cmd_geoFile');
bouton.gfPrevPopOut = uicontrol('parent',bouton.geoFilePanel,'Style','pushbutton','Units','normalized', ...
    'position',[0.9 0.9 .1 .1],'string','pop out', ...
    'Callback',['h=allchild(bouton.gfPreview);',...
    'cmd_copyFig;']);
  
para.geo(1) = 2;
para.geo(2) = 1;
para.chggeo=0;

% escoge la geometria
bouton.strgeo1 ={'Espacio completo','Semi Espacio','N_layer/HS'};
bouton.strgeo2 ={'2 contornos horizontales','placa ilimitada','semi-placa L','semi-placa R'};
bouton.strgeo3 ={'Contornos 2D (X<0:Z_revolve)',...
          'placa ilimitada axisimétrico',...
          'semi-placa L axisimétrico',...
          'semi-placa R axisimétrico'};
bouton.strgeo4 ={'From STL file(s)'};
uicontrol('parent',uigm,'Style','text'     ,'BackgroundColor',c_3  ,'Units','normalized','position',[x0 0.05 dxi .9],'string','');
info.geo    =uicontrol('parent',uigm,'Style','text'     ,'BackgroundColor',c_3  ,'Units','normalized','position',[x0 0.05 dxi .7],'string','El medio 1 se define como');
uicontrol('parent',uigm,'Style','text'     ,'BackgroundColor',c_3  ,'Units','normalized','position',[x0 0.05 4.2*dxi .2],'string',' ');
bouton.geo  =uicontrol('parent',uigm,'Style','popupmenu','BackgroundColor',c_2,'Units','normalized','position',[x0 0.1 4.2*dxi .1],...
    'string',bouton.strgeo1,'value',2,'Callback','para.geo(get(bouton.med,''value''))=get(bouton.geo,''value'');para.chggeo=1;cmd_med;');

% posicion del punto izquierdo del contorno
para.cont(1,1).xa=-3;
para.cont(1,1).a = 6;
para.cont(1,1).th= 0;

para.cont(2,1).a = 1;
para.cont(2,1).xa=-1;
para.cont(2,1).za= 0;
para.cont(2,1).th= 0;

% para.cont(1,1).geo_if3D_axi0_STL1 = [];
% para.cont(2,1).geo_if3D_axi0_STL1 = [];
info.xa  =uicontrol('parent',uigm,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[x0+dx 0.55 dxi .3],'string','xa');
bouton.xa=uicontrol('parent',uigm,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[x0+dx 0.3  dxi .25],...
    'string',para.cont(2,1).xa,'Callback','para.cont(get(bouton.med,''value''),1).xa=str2num(get(bouton.xa,''string''));');
info.za  =uicontrol('parent',uigm,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[x0+2*dx 0.55 dxi .3],'string','za');
bouton.za=uicontrol('parent',uigm,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[x0+2*dx 0.3  dxi .25],...
    'string',para.cont(2,1).za,'Callback','para.cont(get(bouton.med,''value''),1).za=str2num(get(bouton.za,''string''));');

% longitud real de la mitad del ancho horizontal
info.a   =uicontrol('parent',uigm,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[x0+3*dx 0.55 dxi .3],'string','1/2 ancho /x (a)');
bouton.a =uicontrol('parent',uigm,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[x0+3*dx 0.3 dxi .25],...
    'string',para.cont(2,1).a,'Callback','para.cont(get(bouton.med,''value''),1).a=str2num(get(bouton.a,''string''));');

% longitud real de la mitad del ancho horizontal
info.th  =uicontrol('parent',uigm,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[x0+4*dx 0.55 dxi .3],'string','angulo/z [-90,90º]');
bouton.th=uicontrol('parent',uigm,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[x0+4*dx 0.3 dxi .25],...
    'string',para.cont(2,1).th,'Callback','para.cont(get(bouton.med,''value''),1).th=str2num(get(bouton.th,''string''));');

% parte de contorno
info.partes = uicontrol('parent',uigm,'Style','text','BackgroundColor',[.7 .7 .7]  ,'Units','normalized','position',[x0+5*dx 0.05 3*dx+0.95*dxi 0.95],'string','partes del contorno');

% escoge el contorno
strcont     = {'arriba','abajo'};
info.cont   = uicontrol('parent',uigm,'Style','text'     ,'BackgroundColor',c_3  ,'Units','normalized','position',[x0+5*dx 0.65 2.2*dxi .2],'string',' ');
bouton.cont = uicontrol('parent',uigm,'Style','popupmenu','BackgroundColor',c_2,'Units','normalized','position',[x0+5*dx 0.55 2.2*dxi .25],...
    'string',strcont,'Callback','cmd_med;');

% geometria de los contornos
for i=1:2
    para.cont(1,i).geom	= 1;
    para.cont(1,i).ba   = 0.25;
    para.cont(1,i).h    = 0.25;
end
strcontgeo= ...
    {'Smooth Gaussiana', ...
    'Parabola', ...
    'Triangulo', ...
    'Coseno', ...
    'Elipce', ...
    'Elipce asimetrica', ...
    'Trapecio', ...
    'Base en los bordes'};
for i=1:2
    para.cont(2,i).geom  = 7;
end
info.contgeo    =uicontrol('parent',uigm,'Style','text'     ,'BackgroundColor',c_3  ,'Units','normalized','position',[x0+5*dx 0.30 2.2*dxi .3],'string','Tipo de contorno');
bouton.contgeo  =uicontrol('parent',uigm,'Style','popupmenu','BackgroundColor',c_2,'Units','normalized','position',[x0+5*dx 0.20 2.2*dxi .25],...
    'string',strcontgeo,'Callback','para.cont(get(bouton.med,''value''),get(bouton.cont,''value'')).geom=get(bouton.contgeo,''value'');cmd_med;');

% base
for i=1:2
    para.cont(2,i).ba  = 0.6;
end
info.base   =uicontrol('parent',uigm,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[x0+7*dx 0.55 dxi .3],'string','1/2 ancho de base');
bouton.base =uicontrol('parent',uigm,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[x0+7*dx 0.3 dxi .25],...
    'string',para.cont(2,1).ba,'Callback','para.cont(get(bouton.med,''value''),get(bouton.cont,''value'')).ba=str2num(get(bouton.base,''string''));');

% altura (hauteur)
para.cont(2,1).h  =-0.2; %monte
para.cont(2,2).h  = 0.5; %valle
info.haut   =uicontrol('parent',uigm,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[x0+8*dx 0.55 dxi .28],'string','altura');
bouton.haut =uicontrol('parent',uigm,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[x0+8*dx 0.3 dxi .25],...
    'string',para.cont(2,1).h,'Callback','para.cont(get(bouton.med,''value''),get(bouton.cont,''value'')).h=str2num(get(bouton.haut,''string''));');

% geometria de la rugosidad

strrug= {'Plano','Sinus','Triangulo'};
for j=1:2
    for i=1:2
        para.cont(j,i).ruggeo= 1;
        para.cont(j,i).rba   = 0.25;
        para.cont(j,i).rh    = 0.25;
    end
end
info.ruggeo    =uicontrol('parent',uigm,'Style','text'     ,'BackgroundColor',c_3  ,'Units','normalized','position',[x0+9.1*dx 0.72 3*dxi .16],'string','');
bouton.ruggeo  =uicontrol('parent',uigm,'Style','popupmenu','BackgroundColor',c_2,'Units','normalized','position',[x0+9.1*dx 0.72 3*dxi .15],...
    'string',strrug,'Callback','cmd_rug;cmd_med;');

% rugosidad / base
info.rugbase   =uicontrol('parent',uigm,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[x0+9.5*dx 0.55 dxi .1],'string','ancho');
bouton.rugbase =uicontrol('parent',uigm,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[x0+9.5*dx 0.3 dxi .25],...
    'string',para.cont(1,1).rba,'Callback','para.cont(get(bouton.med,''value''),get(bouton.cont,''value'')).rba=str2num(get(bouton.rugbase,''string''));cmd_med;');

% rugosidad / hauteur
info.rughaut   =uicontrol('parent',uigm,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[x0+10.5*dx 0.55 dxi .1],'string','altura');
bouton.rughaut =uicontrol('parent',uigm,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[x0+10.5*dx 0.3 dxi .25],...
    'string',para.cont(2,1).rh,'Callback','para.cont(get(bouton.med,''value''),get(bouton.cont,''value'')).rh=str2num(get(bouton.rughaut,''string''));');

%parametros DWN
info.DWN       = uicontrol('parent',uigm,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.1 0.5 .25 .4],'string','parametros DWN');

% Ahora kmax se hace para cada frecuencia 1.5 veces el polo de Rayleigh ó
% 1/3 de correspondiente a la fmax (lo que de más).
para.DWNkmax   = 20;%ke
% info.DWNkmax   = uicontrol('parent',uigm,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.1 0.7 .1 .1],'string','kx max');
% bouton.DWNkmax = uicontrol('parent',uigm,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.1 0.55 .1 .15],...
%     'string',para.DWNkmax,'Callback','para.DWNkmax=str2num(get(bouton.DWNkmax,''string''));cmd_DWN;');

para.DWNxl = 1000;
info.DWNxl   = uicontrol('parent',uigm,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.1 0.7 .1 .1],'string','xl');
bouton.DWNxl = uicontrol('parent',uigm,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.1 0.55 .1 .15],...
    'string',para.DWNxl,'Callback','para.DWNxl=str2num(get(bouton.DWNxl,''string''));');

para.DWNnbptkx  = 12000;
info.DWNnbptkx 	= uicontrol('parent',uigm,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.25 0.7 .1 .1],'string','nbpt kx');
bouton.DWNnbptkx= uicontrol('parent',uigm,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.25 0.55 .1 .15],...
    'string',para.DWNnbptkx,'Callback','para.DWNnbptkx=str2num(get(bouton.DWNnbptkx,''string''));');

% DX              = pi/para.DWNkmax;
% xmax            = DX*para.DWNnbptkx/2;
% DK              = para.DWNkmax/(para.DWNnbptkx*pi);
% strinfoDWN      = {['DX=',num2str(DX,2),'  Xmax=',num2str(xmax,2)];['DK=',num2str(DK,2)]};

% strinfoDWN = '';
% info.infoDWN 	= uicontrol('parent',uigm,'Style','text','BackgroundColor',c_4  ,'Units','normalized',...
%     'position',[0.1 0.26 .25 .27],'string',strinfoDWN);
info.estrDWN  = uipanel('parent',uigm,'title','Estratificación',...
    'HighlightColor',[1 1 1],'BackgroundColor',c_1,...
    'Position',[.351 0 .649 0.99]);
bouton.axe_estrDWN=axes('parent',info.estrDWN,'Units','normalized',...
    'OuterPosition',[0 0 1 1],'Layer','top',...
    'HitTest','off');
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% caracteristica fuente 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
h0rt(2)= uitab('parent',h0rtDad,'title','Fuente y Receptores');
uif     = uipanel('parent',h0rt(2),'Title','Caracteristica de la fuente',...
    'HighlightColor',[1 1 1],...
    'Position',[.005 .5 .99 .5]);

% tipo de fuente
uicontrol('parent',uif,'Style','text','BackgroundColor',c_3,...
    'Units','normalized','position',[0.01 0.85 .15 .15],'string','Tipo de fuentes:');
para.fuente=1;
strfuente ={'Ondas Planas','Fuentes Puntuales'};
bouton.fuente=uicontrol('parent',uif,'Style','popupmenu',...
    'BackgroundColor',c_2,'Units','normalized',...
    'position',[1.5146e-01 7.5704e-01 2.5000e-01 2.4000e-01],...
    'string',strfuente,'Callback',...
    'para.fuente=get(bouton.fuente,''value'');para.ninc=1;cmd_fuente;');
% Numero de fuentes
info.ninc     =uicontrol('parent',uif,'Style','text','BackgroundColor',c_3,...
  'Units','normalized','position',[0.4 0.85 .10 .15],'string','# fuentes:');
para.ninc     =1;
bouton.ninc   =uicontrol('parent',uif,'Style','edit','BackgroundColor',c_2,...
  'Units','normalized','position',[0.5 0.85 .07 .15],...
    'string',para.ninc,'Callback','para.ninc=str2num(get(bouton.ninc,''string''));cmd_ninc;');
% selecion de la onda # i
info.inc    =uicontrol('parent',uif,'Style','text','BackgroundColor',c_3,...
  'Units','normalized','position',[0.01 0.5 .07 .24],'string','Fuente #:');
para.strinc ='1';
bouton.inc  =uicontrol('parent',uif,'Style','popupmenu','BackgroundColor',c_2,...
  'Units','normalized','position',[0.00 0.25 .10 .23],...
    'string',para.strinc,'Callback','cmd_inc;');
% tipo de onda plana
para.tipo_onda=1;
info.tipo_onda=uicontrol('parent',uif,'Style','text','BackgroundColor',c_3,...
  'Units','normalized','position',[0.09 0.5 .07 .24],'string','Tipo de onda');
strtipo_onda ={'P','SV','R'};
bouton.tipo_onda=uicontrol('parent',uif,'Style','popupmenu','BackgroundColor',...
  c_2,'Units','normalized','position',[0.08 0.25 .15 .23],...
    'string',strtipo_onda,'Callback','para.tipo_onda(get(bouton.inc,''value''))=get(bouton.tipo_onda,''value'');');
% polarisation de ondas en anlalisis 2D  (TODO: remover y usar tipo_onda)
para.pol    = 1;
strpol ={'SH','P/SV'};
bouton.pol=uicontrol('parent',uif,'Style','popupmenu',...
    'BackgroundColor',c_2,'Units','normalized',...
    'position',[0.08 0.25 .15 .23],'string',strpol,...
    'Callback','para.pol=get(bouton.pol,''value'');cmd_pol;');

% Posicion de la fuente o origen de la fase de las ondas
info.xzs       =uicontrol('parent',uif,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.26  0.375 .15 .12],'string','Origen fases');
               uicontrol('parent',uif,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.26  0.25  .15 .12],'string','');
               uicontrol('parent',uif,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.26  0.25  .05 .12],'string','Xs');
info.ys        =uicontrol('parent',uif,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.31  0.25  .05 .12],'string','Ys','visible','off');
               uicontrol('parent',uif,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.36  0.25  .05 .12],'string','Zs');
para.xs(1)= 0;
para.ys(1)= 0;
para.zs(1)= 0;
bouton.xs      =uicontrol('parent',uif,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.26  0.01  .05 .23],...
    'string',para.xs(1),'Callback','para.xs(get(bouton.inc,''value''))=str2num(get(bouton.xs,''string'')); ');
bouton.ys      =uicontrol('parent',uif,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.31 0.01  .05 .23],'visible','off', ...
    'string',para.ys(1),'Callback','para.ys(get(bouton.inc,''value''))=str2num(get(bouton.ys,''string'')); ');
bouton.zs      =uicontrol('parent',uif,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.36 0.01  .05 .23],...
    'string',para.zs(1),'Callback','para.zs(get(bouton.inc,''value''))=str2num(get(bouton.zs,''string'')); ');

info.orient =uicontrol('parent',uif,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.43 0.375 .15 .12],'string','Orientacion');
% Angulo de las ondas planas incindentes
info.gam    =uicontrol('parent',uif,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.43 0.25 .07 .12],'string','q','fontname','symbol');
para.gam(1) =0;
bouton.gam  =uicontrol('parent',uif,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.43 0.01 .07 .23],...
    'string',para.gam(1),'Callback','para.gam(get(bouton.inc,''value''))=str2num(get(bouton.gam,''string''));');
% Angulo de las ondas planas incindentes
info.phi    =uicontrol('parent',uif,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.51 0.25 .07 .12],'string','f','fontname','symbol');%,'visible','off');
para.phi(1) =0;
bouton.phi  =uicontrol('parent',uif,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.51 0.01 .07 .23],...
    'string',para.phi(1),'Callback','para.phi(get(bouton.inc,''value''))=str2num(get(bouton.phi,''string''));','visible','off');


% para escoger si se pasa al tiempo las senales (resultados mas volumicos)
para.spct=0;
bouton.spct = uicontrol('parent',uif,'Style','checkbox',...
    'BackgroundColor',c_2,'Units','normalized','position',[0.705 0.01 .25 .2],...
    'string','Solo espectros','Callback','para.spct=get(bouton.spct,''value'');cmd_spec;');

% tipo profil temporal de la fuente
info.ptf    = uicontrol('parent',uif,'Style','text'     ,...
    'BackgroundColor',c_3  ,'Units','normalized',...
    'position',[0.705 0.65 .25 .3],'string','Amplitude function');
strpulsotps ={'File (~,t0,~)',...
              'Ricker (tp,t0,~) [duracion total]',...
              'Ricker (tp,ts,t0) [Characteristic T]',...
              'Gaussiana (tp,ts,f_break)',...
              'Gaussiana (rt,t0,~)',...
              'butterworth (n,Wn)',...
              'Dirac (~,t0,~)'};
para.pulso.tipo   = 3;
bouton.pulsotps =uicontrol('parent',uif,'Style','popupmenu',...
    'BackgroundColor',c_2,'Units','normalized','position',[0.705 0.6 .25 .2],...
    'string',strpulsotps,'Callback','para.pulso.tipo=get(bouton.pulsotps,''value'');');

% duracion del pulso
uicontrol('parent',uif,'Style','text'     ,...
    'BackgroundColor',c_3  ,'Units','normalized',...
    'position',[0.705 0.4 .07 .2],'string','(a)');
para.pulso.a   = 1;
bouton.Ricker_tp  =uicontrol('parent',uif,'Style','edit',...
    'BackgroundColor',c_2,'Units','normalized',...
    'position',[0.705 0.2 .07 .2],...
    'string',para.pulso.a,'Callback',...
    'para.pulso.a=str2num(get(bouton.Ricker_tp,''string''));');

% duracion del retraso
uicontrol('parent',uif,'Style','text',...
    'BackgroundColor',c_3  ,'Units','normalized',...
    'position',[0.785 0.4 .07 .2],'string','(b)');
para.pulso.b = 2;
bouton.delais   =uicontrol('parent',uif,'Style','edit',...
    'BackgroundColor',c_2,'Units','normalized',...
    'position',[0.785 0.2 .07 .2],...
    'string',para.pulso.b,'Callback',...
    'para.pulso.b=str2num(get(bouton.delais,''string''));');
  
% frecuencia inicial
uicontrol('parent',uif,'Style','text',...
    'BackgroundColor',c_3  ,'Units','normalized',...
    'position',[0.865 0.4 .07 .2],'string','(c)');
para.pulso.c = 0;
bouton.GaussSta  =uicontrol('parent',uif,'Style','edit',...
    'BackgroundColor',c_2,'Units','normalized',...
    'position',[0.865 0.2 .07 .2],...
    'string',para.pulso.c,'Callback',...
    'para.pulso.c=str2num(get(bouton.GaussSta,''string''));');

bouton.pulsoShow = uicontrol('parent',uif,'Style',...
    'pushbutton','Units','normalized', ...
    'position',[0.94 0.4 .06 .2],'string','plot', ...
    'Callback','showPulso(para);');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% posicion de los receptores 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uirt0(1)     = uipanel('parent',h0rt(2),'Title','Receptores: Posiciones','HighlightColor' ,...
    [1 1 1],'BackgroundColor',c_1,'Position',[.005 .0 .99 .5]);

para.rec.resatboundary=0;
para.rec.resatboundaryDecimate=20;
bouton.resatboundary = uicontrol('parent',uirt0(1),'Style','checkbox',...
    'BackgroundColor',c_2,'Units','normalized','position',[0.01 0.7 .5 .20],...
    'string','Receptores en puntos de la frontera    Decimate fact:',...
    'Callback',['para.rec.resatboundary=get(bouton.resatboundary,''value'');',...
    ' ']);
bouton.resatboundaryDecimate  =uicontrol('parent',uirt0(1),'Style','edit',...
    'BackgroundColor',c_2,'Units','normalized','position',[0.51 0.7 .1 .20],...
    'string',para.rec.resatboundaryDecimate,'Callback',[...
    'para.rec.resatboundaryDecimate=str2num(get(bouton.resatboundaryDecimate,''string''));',...
    ' ']);
uicontrol('parent',uirt0(1),'Style','text','BackgroundColor',c_2  ,...
    'Units','normalized','position',[0.61 0.71 .07 .16],'string',' Scale:');
para.rec.resatboundaryScale = 1;
bouton.resatboundaryScale =uicontrol('parent',uirt0(1),'Style','edit',...
    'BackgroundColor',c_2,'Units','normalized','position',[0.68 0.7 .05 .2],...
    'string',para.rec.resatboundaryScale,'Callback',[...
    'para.rec.resatboundaryScale=str2num(get(bouton.resatboundaryScale,''string''));',...
    ' ']);
uicontrol('parent',uirt0(1),'Style','text','BackgroundColor',c_2  ,...
    'Units','normalized','position',[0.73 0.71 .1 .16],'string',' for Media:');
para.rec.resatboundaryScaleMediaRange = '1:end';
bouton.resatboundaryScaleMediaRange =uicontrol('parent',uirt0(1),'Style','edit',...
    'BackgroundColor',c_2,'Units','normalized','position',[0.83 0.7 .15 .20],...
    'string',para.rec.resatboundaryScaleMediaRange,'Callback',[...
    'para.rec.resatboundaryScaleMediaRange=get(bouton.resatboundaryScaleMediaRange,''string'');',...
    ' ']);

% calculo de las posiciones de los receptores
para.recpos =1;
para.chgrec      =0;
strrecpos ={'En la superficie','Malla constante','Posicion libre'};
bouton.recpos=uicontrol('parent',uirt0(1),'Style','popupmenu','BackgroundColor',c_2,'Units','normalized','position',[0.01 0.45 .22 .22],...
    'string',strrecpos,'Callback','para.chgrec=1;para.recpos=get(bouton.recpos,''value'');cmd_recpos;');

info.recx  =uicontrol('parent',uirt0(1),'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.01 0.3  .22 .15],'string','X');
info.recy  =uicontrol('parent',uirt0(1),'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.01 0.15 .22 .15],'string','Y');
info.recz  =uicontrol('parent',uirt0(1),'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.01 0.0  .22 .15],'string','Z');

info.rec0  =uicontrol('parent',uirt0(1),'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.25 0.45 .22 .15],'string','Posicion primer receptor');
info.recdr =uicontrol('parent',uirt0(1),'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.5  0.45 .22 .15],'string','dr entre receptores');
info.recnr  =uicontrol('parent',uirt0(1),'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.75 0.45 .22 .15],'string','numero de estaciones');

% posicion del primer receptor
para.rec.xri=-1.5;       %posicion del primer receptor
bouton.xri  =uicontrol('parent',uirt0(1),'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.25 0.3  .22 .15],...
    'string',para.rec.xri,'Callback','para.rec.xri=str2num(get(bouton.xri,''string'')); ');

para.rec.yri= 0;       %posicion del primer receptor
bouton.yri  =uicontrol('parent',uirt0(1),'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.25 0.15 .22 .15],...
    'string',para.rec.yri,'Callback','para.rec.yri=str2num(get(bouton.yri,''string'')); ');

para.rec.zri= 0;       %posicion del primer receptor
bouton.zri  =uicontrol('parent',uirt0(1),'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.25 0.0  .22 .15],...
    'string',para.rec.zri,'Callback','para.rec.zri=str2num(get(bouton.zri,''string'')); ');

% paso entre los receptores
para.rec.dxr	= 0.3;
bouton.dxr  =uicontrol('parent',uirt0(1),'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.5 0.3  .22 .15],...
    'string',para.rec.dxr,'Callback','para.rec.dxr=str2num(get(bouton.dxr,''string'')); ');

para.rec.dyr	= 0.2;
bouton.dyr  =uicontrol('parent',uirt0(1),'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.5 0.15 .22 .15],...
    'string',para.rec.dyr,'Callback','para.rec.dyr=str2num(get(bouton.dyr,''string'')); ');

para.rec.dzr	= 0.2;
bouton.dzr  =uicontrol('parent',uirt0(1),'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.5 0.0  .22 .15],...
    'string',para.rec.dzr,'Callback','para.rec.dzr=str2num(get(bouton.dzr,''string'')); ');

% numero de estaciones
para.rec.nrecx= 11;
para.chgnrecx = 0;
bouton.nrecx  =uicontrol('parent',uirt0(1),'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.75 0.3  .22 .15],...
    'string',para.rec.nrecx,'Callback','para.rec.nrecx=str2num(get(bouton.nrecx,''string''));para.chgnrecx=1;cmd_recpos; ');

para.rec.nrecy= 1;
bouton.nrecy  =uicontrol('parent',uirt0(1),'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.75 0.15 .22 .15],...
    'string',para.rec.nrecy,'Callback','para.rec.nrecy=str2num(get(bouton.nrecy,''string'')); ');

para.rec.nrecz= 1;
bouton.nrecz  =uicontrol('parent',uirt0(1),'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.75 0.0  .22 .15],...
    'string',para.rec.nrecz,'Callback','para.rec.nrecz=str2num(get(bouton.nrecz,''string'')); ');

% se trata el receptor #
irec    = 1;
strrec	= 1:para.rec.nrecx;
info.irec =uicontrol('parent',uirt0(1),'Style','text','BackgroundColor',c_3,'Units','normalized','position',[0.75 0.15 .22 .15],'string','rec #');
bouton.irec=uicontrol('parent',uirt0(1),'Style','popupmenu','BackgroundColor',c_2,'Units','normalized','position',[0.75 0.0 .22 .15],...
    'string',strrec,'Callback','cmd_recpos;');

%%%%%%%%%%%%%%%%%%%%%%%
%% variables de salida 
%%%%%%%%%%%%%%%%%%%%%%%
  
%deplacements
para.sortie.Ux	= 1;
para.sortie.Uy  = 1;
para.sortie.Uz  = 1;

para.sortie.Ut  = 1;
para.sortie.UPh = 0;
para.sortie.USh = 0;
para.sortie.UIh = 0;
para.sortie.UPt = 0;
para.sortie.USt = 0;

% %vitesses
% para.sortie.Vx  =0;
% para.sortie.Vy  =0;
% para.sortie.Vz  =0;

%contrainte
para.sortie.sxx =1;
para.sortie.syy =1;
para.sortie.szz =1;
para.sortie.sxy =1;
para.sortie.sxz =1;
para.sortie.syz =1;

% %pression
% para.sortie.P	=0;
% %deformation
% para.sortie.exx =0;
% para.sortie.eyy =0;
% para.sortie.ezz =0;
% % para.sortie.exy =0;
% para.sortie.exz =0;
% % para.sortie.eyz =0;

fn0     = fieldnames(para.sortie);
eps     = 1/100;
n       = length(fn0);


h0rt(3)= uitab('parent',h0rtDad,'title','Parametros de calculo');
uipc(1)     = uipanel('parent',h0rt(3),'Title','Variables de Salida',...
    'HighlightColor' , [1 1 1],...
    'BackgroundColor',c_1,'Position',[.005 .5 .99 .5]);
%deplacements 
%unite cf echelle espace
boutonsal.Ux = uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[0/8 2/4 1/8-eps 1/4-eps],'value',1);
boutonsal.Uy = uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[1/8 2/4 1/8-eps 1/4-eps],'value',1);
boutonsal.Uz = uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[2/8 2/4 1/8-eps 1/4-eps],'value',1);
% boutonsal.G  = uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[3/8 2/4 1/8-eps 1/4-eps],'value',1);

boutonsal.Ut = uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[0/8 3/4 1/8-eps 1/4-eps],'value',1,'TooltipString','campo total','Visible','off');
boutonsal.UPh= uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[1/8 3/4 1/8-eps 1/4-eps],'value',0,'TooltipString','campo P homogeneo; ! calculacion lenta');
boutonsal.USh= uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[2/8 3/4 1/8-eps 1/4-eps],'value',0,'TooltipString','campo S homogeneo; ! calculacion lenta');
boutonsal.UIh= uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[3/8 3/4 1/8-eps 1/4-eps],'value',0,'TooltipString','campo P&S inhomogeneo; ! calculacion lenta');
boutonsal.UPt= uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[5/8 3/4 1/8-eps 1/4-eps],'value',0,'TooltipString','campo P total');
boutonsal.USt= uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[6/8 3/4 1/8-eps 1/4-eps],'value',0,'TooltipString','campo S total');


% %vitesses
% %unite cf echelle espace / echelle tps
% boutonsal.Vx=uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[3/8 2/4 1/8-eps 1/4-eps]);
% boutonsal.Vy=uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[4/8 2/4 1/8-eps 1/4-eps]);
% boutonsal.Vz=uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[5/8 2/4 1/8-eps 1/4-eps]);

%contrainte
%unite en GPa directement
boutonsal.sxx=uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[0/8 1/4 1/8-eps 1/4-eps]);
boutonsal.syy=uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[1/8 1/4 1/8-eps 1/4-eps]);
boutonsal.szz=uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[2/8 1/4 1/8-eps 1/4-eps]);
boutonsal.sxy=uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[3/8 1/4 1/8-eps 1/4-eps]);
boutonsal.sxz=uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[4/8 1/4 1/8-eps 1/4-eps]);
boutonsal.syz=uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[5/8 1/4 1/8-eps 1/4-eps]);
% boutonsal.S  =uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[6/8 1/4 1/8-eps 1/4-eps]);
%pression
% boutonsal.P  =uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[7/8 1/4 1/8-eps 1/4-eps]);
%deformation
% boutonsal.exx=uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[0/8 0/4 1/8-eps 1/4-eps]);
% boutonsal.eyy=uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[1/8 0/4 1/8-eps 1/4-eps]);
% boutonsal.ezz=uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[2/8 0/4 1/8-eps 1/4-eps]);
% % boutonsal.exy=uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[3/8 0/4 1/8-eps 1/4-eps]);
% boutonsal.exz=uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[4/8 0/4 1/8-eps 1/4-eps]);
% % boutonsal.eyz=uicontrol('parent',uipc(1),'Style','checkbox','Units','normalized','position',[5/8 0/4 1/8-eps 1/4-eps]);

for ifn0=1:n
    tmp=fn0(ifn0);
    nameb0 =['boutonsal.',tmp{1}];
    nameb2=['para.sortie.',tmp{1}];
    strcllbck=[nameb2,'=get(',nameb0,',''value'');'];
    set(eval(nameb0),'BackgroundColor',c_2,'string',tmp{1},'callback',strcllbck);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parametros de calculo 
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%h0rt(3)= uitab('parent',h0rtDad,'title','Parametros de calculo');
uipc     = uipanel('parent',h0rt(3),'Title','Parametros de calculo',...
    'HighlightColor' , [1 1 1],...
    'BackgroundColor',c_1,'Position',[.005 .0 .99 .5]);

% numero de puntos por longitud de onda a emplear en las discretizaciones.
             uicontrol('parent',uipc,'Style','text'     ,'BackgroundColor',c_3  ,'Units','normalized','position',[0.01 0.55 .15 .4],'string','# de puntos / longitud de onda');
para.npplo	 =6;       %posicion del primer receptor
bouton.npplo =uicontrol('parent',uipc,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.01 0.05 .15 .4],...
    'string',para.npplo,'Callback','para.npplo=str2num(get(bouton.npplo,''string''));');

% frecuencia maxima
             uicontrol('parent',uipc,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.21 0.55 .15 .4],'string','Fq max');
para.fmax	 = 3;
bouton.fmax  =uicontrol('parent',uipc,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.21 0.05 .15 .4],...
    'string',para.fmax,'Callback','para.fmax=str2num(get(bouton.fmax,''string'')); cmd_tmax');

% numero de frecuencias en el espectro a calcular
             uicontrol('parent',uipc,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.41 0.55 .15 .4],'string','N Fq');
para.nf      = 200;
bouton.nf    =uicontrol('parent',uipc,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.41 0.05 .15 .4],...
    'string',para.nf,'Callback','para.nf=str2num(get(bouton.nf,''string'')); cmd_tmax;');

  % ahora DWNomei se estima a partir del la ventana de tiempo de interés
  % como en Bouchon 2003. Se pregunta la ventana de graficación
para.DWNomei	= 0.001;
% info.DWNomei 	= uicontrol('parent',uipc,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.61 0.55 .15 .4],'string','wi');
% bouton.DWNomei  = uicontrol('parent',uipc,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.61 0.05 .15 .4],...
%     'string',para.DWNomei,'Callback','para.DWNomei=str2num(get(bouton.DWNomei,''string''));');

para.tmax = (para.zeropad-1)/(para.fmax/(para.nf/2)*para.zeropad);
para.tmaxinteres = para.tmax;
info.tmax = uicontrol('parent',uipc,'Style','text','BackgroundColor',c_3  ,'Units','normalized','position',[0.61 0.55 .15 .4],'string',...
  {['tmax = ' num2str(para.tmax)];'';'T interes:'});

bouton.tmax  = uicontrol('parent',uipc,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.61 0.05 .15 .4],...
    'string',para.tmaxinteres,'Callback','para.tmaxinteres=str2num(get(bouton.tmax,''string''));');


% para.matfor='fortran';
% lstmatfor={'matlab','fortran'};
% bouton.lstmatfor=uicontrol('parent',uitint0(2),'Style','listbox','BackgroundColor',c_2,'Units','normalized','position',[0.01 0.1 .35 .98],...
%     'string',lstmatfor,'value',2,'Callback','para.matfor=lstmatfor{get(bouton.lstmatfor,''value'')};');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dibujo de la configuracion geometrica 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


bouton.axe_conf_geo=axes('parent',uicg,'Units','normalized',...
    'OuterPosition',[0 0 1 1],'ButtonDownFcn','rafraichi');
% para.han.uicg = uicg;
% para.han.axe_conf_geo = axe_conf_geo;

uicontrol('parent',uicg,'Style','pushbutton','Units','normalized', ...
    'position',[0.9 0.9 .1 .1],'string','pop out', ...
    'Callback',['h=allchild(gca);',...
    'cmd_copyFig;']);
  
bouton.rotDespConfGeom = uicontrol('parent',uicg,'Style','popupmenu','Units','normalized', ...
    'position',[0  0.9 .1 .1],'string',{'---' 'rotar' 'desplazar'}, ...
    'Callback','cmd_rotDespConfGeom');

bouton.rafraichiEveryTime = uicontrol('parent',uicg,'Style','popupmenu','Units','normalized',...
    'position',[0.4 0.9 0.2 0.1],'string',{'Do not refresh','Refresh everytime','Refresh now'}, ...
    'value',2,...
    'Callback',['tmp=get(bouton.rafraichiEveryTime,''value'');',...
                'if(tmp==2) rafraichi; end']);

bouton.verNormales = uicontrol('parent',uicg,'Style','checkbox',...
    'BackgroundColor',c_2,'Units','normalized','position',[0 0.8 0.05 0.1],...
    'string','n','value',0);

bouton.verReceptores = uicontrol('parent',uicg,'Style','checkbox',...
    'BackgroundColor',c_2,'Units','normalized','position',[0 0.7 0.05 0.1],...
    'string','p','value',1);

bouton.verGeometria = uicontrol('parent',uicg,'Style','checkbox',...
    'BackgroundColor',c_2,'Units','normalized','position',[0 0.6 0.05 0.1],...
    'string','g','value',1);

drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%
%% Entorno de trabajo 
%%%%%%%%%%%%%%%%%%%%%%%%
uiwk= uitab('parent',h0rtDad,'title','Entorno de trabajo');

%para escoger el repertorio de trabajo
uicontrol('parent',uiwk,'Style','pushbutton','BackgroundColor',c_2,'Units','normalized', ...
    'position',[0.01 0.75 .2 .2],'string','ruta de trabajo', ...
    'Callback',['para.nomrep=chercherep;nomrep=para.nomrep;set(bouton.breptex,''string'',para.nomrep);', ...
    'fid = fopen(''chemin'',''w'');fprintf(fid,''%s'',para.nomrep);fclose(fid);']);

bouton.breptex	=uicontrol('parent',uiwk,'Style','edit','BackgroundColor',[1 1 1],...
    'Units','normalized','position',[0.25 0.75 .74 .2],'string',para.nomrep);
clear chemin;
if exist('chemin','file')==2
    fid = fopen('chemin');
    t   = textscan(fid,'%s'); t = t{1};
    fclose(fid);clear fid;
    for i=1:length(t)
        para.nomrep=t{i};
        if exist(para.nomrep,'dir')==7
            nomrep=para.nomrep;
            set(bouton.breptex,'string',t{i});
            break
        end
    end
    
end; clear t

%pour enregistrer une configuration de parametre
uicontrol('parent',uiwk,'Style','pushbutton','BackgroundColor',c_2,'Units','normalized','position',[0.01 .52  .4 .1],'string','grabar parametros','Callback','savepara;');
uicontrol('parent',uiwk,'Style','pushbutton','BackgroundColor',c_2,'Units','normalized','position',[0.01 .42  .4 .1],'string','cargar parametros default','Callback','pdef=1;chercheconfig;');
uicontrol('parent',uiwk,'Style','pushbutton','BackgroundColor',c_2,'Units','normalized','position',[0.01 .32  .4 .1],'string','cargar parametros desde archivo','Callback','pdef=0;chercheconfig;');
%uicontrol('parent',uiwk,'Style','pushbutton','BackgroundColor',c_2,'Units','normalized','position',[0.01 .22  .4 .1],'string','¿seguro de los parametros ?','Callback','cmd_fig_1;');

%%%%%%%%%%%%%%%%
%% Resultados 
%%%%%%%%%%%%%%%%
bouton.uiwr= uitab('parent',h0rtDad,'title','Resultados');
%pour afficher le parametre de compression
para.compression=20;
uicontrol('parent',bouton.uiwr,'Style','text','BackgroundColor',c_3,'Units','normalized','position',[0.01 0.9 .15 .09],'string','compresion');
bouton.comp=uicontrol('parent',bouton.uiwr,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.01 0.8 .15 .1], ...
    'string',num2str(para.compression),'callback','para.compression=str2num(get(bouton.comp,''string''));');
%pour afficher la densite de point ds le spectre
para.densitept=5;
uicontrol('parent',bouton.uiwr,'Style','text','BackgroundColor',c_3,'Units','normalized','position',[0.2 0.9 .15 .09],'string','densidad pt');
bouton.densitept=uicontrol('parent',bouton.uiwr,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.2 0.8 .15 .1], ...
    'string',num2str(para.densitept),'callback','para.densitept=str2num(get(bouton.densitept,''string''));');

%pour afficher le spectre
para.b_dib(1).dessinsp=1;
bouton.dessinsp=uicontrol('parent',bouton.uiwr,'Style','checkbox','BackgroundColor',c_2,'Units','normalized','position',[0.01 0.45 .35 .12], ...
    'string','dibujo del espectro','value',para.b_dib(1).dessinsp,'callback','para.b_dib(1).dessinsp=get(bouton.dessinsp,''value'');');
%pour afficher en normalise
para.b_dib(1).normalise=1;
bouton.normalise=uicontrol('parent',bouton.uiwr,'Style','checkbox','BackgroundColor',c_2,'Units','normalized','position',[0.01 0.60 .35 .12], ...
    'string','dibujo normalizado','value',para.b_dib(1).normalise,'callback','para.b_dib.normalise=get(bouton.normalise,''value'');');

%pour tracer un autre resultat
bouton.binitname =uicontrol('parent',bouton.uiwr,'Style','pushbutton','BackgroundColor',c_2,'Units','normalized','position',[0.12 0.25 .28 .12], ...
    'string','Dibujar otros resultados','Callback','bdibcheck;[RESULT,para.b_dib]=callbackdessin(para,para.nomcarpeta,bouton,para.b_dib);');
%pour exporter en txt                                          utc,uw,stc,sw,name,cont1
bouton.bexpok	=uicontrol('parent',bouton.uiwr,'Style','pushbutton','BackgroundColor',c_2,'Units','normalized','position',[0.01 0.25 .1 .12], ...
    'string',{'-> ascii'},'Callback','if (exist(''para'') && exist(''RESULT'')); sortietext(para,RESULT.utc);else disp(''rien a exporter'');end;');

%en differentes couleurs
strcoul={'b','k','r','m','c','g'};
uicontrol('parent',bouton.uiwr,'Style','text','BackgroundColor',c_3,'Units','normalized',                     'position',[0.5 0.79 .49 .20]);
uicontrol('parent',bouton.uiwr,'Style','text','BackgroundColor',c_3,'Units','normalized',                     'position',[0.6 0.79 .1 .17],'string','color');
bouton.couleur=uicontrol('parent',bouton.uiwr,'Style','popupmenu','BackgroundColor',c_2,'Units','normalized','position',[0.6 0.73 .1 .15],'string',strcoul);
%pour retracer le resultat en cours
bouton.btracok	=uicontrol('parent',bouton.uiwr,'Style','pushbutton','BackgroundColor',c_2,'Units','normalized','position',[0.71 0.82 .265 .12], ...
    'string','Trazas de Resultados','Callback','if(exist(''RESULT'',''var''));para.redraw=1;bdibcheck;para = dibujo(para,bouton,RESULT);end');
%pour tracer le filme
para.film.filmeRange = eval('1:10:800');
uicontrol('parent',bouton.uiwr,'Style','text','BackgroundColor',c_3,'Units','normalized','position',[0.5 0.34 .49 .43]);
bouton.bfilmeRangeBt =uicontrol('parent',bouton.uiwr,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.50 0.66 .20 0.08],...
    'string','1:10:800','Callback','cmd_filmeRange;');
info.filmeRangeTime = uicontrol('parent',bouton.uiwr,'Style','text','BackgroundColor',c_2,'Units','normalized','position',[0.50 0.56 .20 0.08],'string','');
strfilmstyle = {'style = color','style = grid','style = grid+shadow','Quiver Plot'};
para.film.filmStyle = 1;
para.film.filmeMecElem = 1;
para.film.fps = 30;
para.film.BoundaryWarpRange = '1:end';
bouton.bfilmeStyleDrop =uicontrol('parent',bouton.uiwr,'Style','popupmenu','BackgroundColor',c_2,'Units','normalized','position',[0.5 0.44 .2 .1],...
    'string',strfilmstyle,'Callback','para.film.filmStyle = get(bouton.bfilmeStyleDrop,''value''); cmd_sel_film_mecElem;');
bouton.bfilmeRunBt =uicontrol('parent',bouton.uiwr,'Style','pushbutton','BackgroundColor',c_2,'Units','normalized','position',[0.71 0.56 .265 .19], ...
    'string','Snapshots','Callback','if(exist(''RESULT'',''var''));filmoscopio2(para,RESULT,bouton.inc.Value);end;');
uicontrol('parent',bouton.uiwr,'Style','text','BackgroundColor',c_3,'Units','normalized','position',[0.70 0.46 .08 0.08],...
    'string','FPS:');
bouton.bfilmeFPS =uicontrol('parent',bouton.uiwr,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.76 0.46 .08 0.08],...
    'string','30','Callback','para.film.fps = eval(get(bouton.bfilmeFPS,''string''));');

% para hacer el cálculo
uicontrol('parent',bouton.uiwr,'Style','pushbutton','BackgroundColor',[0 1 0],'Units','normalized','position',[0.01 0.01 .48 .19],'string','calculo del espectro e inversion', ...
    'Callback','para.espyinv=1;bdibcheck;[RESULT,para]=calculo(para);para.name=RESULT.name;'); %utc,uw,stc,sw,name,cont1

% para solamente invertir
uicontrol('parent',bouton.uiwr,'Style','pushbutton','BackgroundColor',c_3,'Units','normalized','position',[0.51 0.01 .48 .19],'string','solo inversion del ultimo espectro', ...
    'Callback','para.espyinv=0;bdibcheck;[RESULT,para]=calculo(para);');


% para gardar en memoria la direccion de la carpeta del programma
para.nomcarpeta=pwd;

cmd_fig_1;

%% ini 
clear c_1 c_2 c_3 c_4 h0fig1 i dx dxi x0
clear  DK DX alpha ans 
clear bet bexpok bfilmeFPS bfilmeRangeBt 
clear bfilmeRunBt bfilmeStyleDrop binitname 
clear btracok cmd_fig_1 
clear eps fn0 h0rt h0rtDad i icont ifn0 irec j lambda 
clear lstatt med mu n nameb0 nameb2 nomrep 
clear onoff onoff1 onoff2 onoff3 onoff4 onoff5 rho strcllbck 
clear strcont strcontgeo strcoul strfilmstyle strfuente 
clear strinfoDWN strmed 
clear strpol strpulsotps strrec strrecpos strrug strsubmed 
clear strtipo_onda tatt tmp uicg uif uigm uipc uipm uirt0 uiwk xmax
clc
rafraichi;
if isempty(gcp('nocreate'))
  parpool
end
set(bouton.rafraichiEveryTime,'value',1)