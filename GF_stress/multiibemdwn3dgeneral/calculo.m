function [RESULTADO,para]=calculo(para)
% programa para calcular la respuesta en frecuencia de una geometria
% compleja con el metodo IBEM
tstart  = tic;
clearfield;
cparll  = ~isempty(gcp('nocreate'));
para.GraficarCadaDiscretizacion = true;
if isempty(gcp('nocreate'))
    % parpool(12) %Max permitido en Tonatiuh
    while cparll==0
%        pause(5)
        cparll  = parpool(4);
    end
end

nametmp        = namefile(para);
para.nametmp   = nametmp;

if ~isfield(para,'siDesktop'); para.siDesktop = false; end

if para.siDesktop
  h = waitbarSwitch(0,'Preparacion integracion de Gauss','CreateCancelBtn','closereq;');
else
  h = 0;
  disp('Preparacion integracion de Gauss');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% abscisas y pesos de integracion gaussiana %
%  quadrature de Gauss-Legendre              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if para.geo(1) ~= 3 % El fondo NO es estratificado con DWN
  if para.dim==1 % SH
    ngau        = 21;
    gaussian       = Gauss_Legendre(ngau);
    para.gaussian  = gaussian;
    gaussex=gaussian;
  elseif para.dim == 4 % 3D gen
    % variable  para.cubature  por renglones: [coordSimplex1 coordSimplex2 W]
    
    % para rij == 0
    % (los siguientes esquemas de cubatura no incluyen el centroide)
%     cmd_setCubatureTriangle6p;
     para.gaussian = [];
     cmd_setCubatureTriangle7p;
%     cmd_setCubatureTriangle12p;
%     para.cubatureex = para.cubature;
%     para.gaussian.ngauex = para.gaussian.ngau;
    
    % para rij ~= 0
%     cmd_setCubatureTriangle4p;
%      cmd_setCubatureTriangle7p;
%      cmd_setCubatureTriangle16p;
%      cmd_setCubatureTriangle24p;
%     gaussex = para.gaussian;
%     para.gaussex = gaussex;
    
    gaussian       = Gauss_Legendre(11);
    gaussex = para.gaussian;
%     para.gaussian  = gaussian;
  else % P-SV y 3D axi
    ngau        = 21;
    gaussex  	= Gauss_Legendre(ngau);
    para.gaussex= gaussex;
    
    ngau        = 11;
    gaussian       = Gauss_Legendre(ngau);
    para.gaussian  = gaussian;
  end
  
if para.geo(1)==3 % Estratificado con DWN
  if para.dim == 4
%     cmd_setCubatureTriangle4p;
%     gaussex = para.gaussian;
%     para.gaussex = gaussex;
  else
    ngau        = 3;
    gaussDWN    = Gauss_Legendre(ngau);
    para.gDWN   = gaussDWN;
  end
end
%%%%%%%%%%%%%%%%%%
%% normalizacion %
%%%%%%%%%%%%%%%%%%
% se normaliza unas variables 
para    = normalizacion(para); 
sal     = para.sortie;
nf      = para.nf;
% para.df = 1/para.tmaxinteres;
para.df = para.fmax/(nf/2);     %paso en frecuencia
df      = para.df;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculo de los espectros e inversion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if para.espyinv==1
  %% discretización inicial de la geometría 2D
  if para.nmed~=1 || para.nmed==1 && para.geo(1)==2
    if para.dim < 4 % 2D,2.5D,3Daxisim
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % dicretizacion de los contornos originales 2D %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if para.siDesktop
        waitbarSwitch(0,h,'Calculo del contorno fino');
      end
      
      %discretiza los contornos iniciales
      para = malla_fina_init_2(para,h);
      
      %calcula los verdaderos contornos, buscando los puntos de intersecion,
      %los contornos que se quedan y los que se quitan
      para = recalculate_contornos(para,h);
      
      %se calcula las coordenadas a mas alta frecuencia para visualizar la
      %continuidad de los phi y probar su interpolacion, cf linea de resample_phi_fv
      % [coord0,listc]  = malla_geom(para);
      % para.listc      = listc;
    elseif para.dim == 4 % 3Dgeneral
      % conectividad de los contornos: para.cont1
      para = getConnectivityStruct(para);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define si usar fuente imagen %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if para.siDesktop
      waitbarSwitch(0.75,h,'Calculo caracteristicas de la fuente');
    end
    if para.geo(1)==2 %semi-espacio
      %a priori se puede tomar en cuenta una fuente imagen
      para.fuenteimagen=1;
      if para.cont(1).ruggeo~=1
        %la superficie del semi-espacio tiene regosidad
        para.fuenteimagen=0;
      else
        % si hay estratos infinitos, no hay fuente imagen
        for m=2:para.nmed
          if para.geo(m)>=2 %presencia de estratos
            para.fuenteimagen=0;
          end
        end
      end
      
      if para.fuente==2 && ( (para.dim==1 && para.pol==2) || para.dim>=3 )
        % para ondas P SV, se considera no mas OP incidentes imagen, no las FP
        para.fuenteimagen=0;
      end
    else
      para.fuenteimagen=0;
    end
  else
    %caso de para.nmed==1 && para.geo(1)=={1,3}
    %en este caso no hay contornos, se va a calcular el campo via DWN o
    %campo incidente en FS no mas
    if isfield(para,'cont1')
      para.cont1=[];
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% identificacion de la fuente %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  para.xzs  :  el medio al que pertenece la fuente
  if para.dim < 4 %2D,2.5D,3Daxisim
    if para.fuente==1 %OP
      para.xzs=ones(1,para.ninc)*inclusiontest(0,1e6,para,1);
    elseif para.fuente==2 %FP
      para.xzs=ones(1,para.ninc);%init
      for iinc=1:para.ninc
        para.xzs(iinc)=inclusiontest(para.xs(iinc),para.zs(iinc),para,1);
      end
    end
  else %3Dgen
    if para.fuente==1 %OP
      % la onda plana surge del medio de fondo siempre:
      para.xzs=ones(1,para.ninc);%*inclusiontest3G(0,0,1e6,para);
    elseif para.fuente==2 %FP
      para.xzs=inclusiontest3G(...
        para.xs(1:para.ninc),...
        para.ys(1:para.ninc),...
        para.zs(1:para.ninc),para);
      disp('Medio de en el que está cada fuente:')
      disp(para.xzs.')
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%
  %% modos de dispersion %
  %%%%%%%%%%%%%%%%%%%%%%%%
  if para.geo(1)==3 && para.fuente==1 && para.nsubmed>1 && ...
      (para.dim==1 && ((max(para.tipo_onda==2)==1 && para.pol==1) || ...
      (max(para.tipo_onda==3)==1 && para.pol==2))) || ...
      (para.dim>=3 && (max(para.tipo_onda==4)==1 || max(para.tipo_onda==5)==1))
    if para.siDesktop
      waitbarSwitch(0,h,'Calculo de las curvas de dispersion');
    end
    
    %         [~,~,k2,w0,ikmax]	= dispersion_curve_kfix_MG_fast2(para);
    j               = 0:nf/2;
    wi              = 2*pi*(j*df + 0.01*df*(j==0));
    %         tic;[vg,f1,ikmax2]   = dispersion_curve_wfix_Haskel_4inv_adapt(para,wi);toc
    [vg,~,f1,ikmax,kx]=dispersion_curve_VP_wfix_Haskel_4inv_u2d2(para,wi);
    
    para.mode.w0    = f1*2*pi;
    para.mode.k2    = kx;
    para.mode.ikmax = ikmax;
    para.mode.vg    = vg;
    %         para.mode.f1    = f1;
    %         para.mode.ikmax2= ikmax2;
    para.ninc       = size(w0,2);
  end
  
  if  para.geo(1)==3 && para.fuente==2 && para.meth_PS==1 && para.nsubmed>1
    %calcul des courbes de dispersion en vue du calcul des IMGIJ dans
    %un multicouche par la methode des correlations
    if para.siDesktop
      waitbarSwitch(0,h,'Calculo de las curvas de dispersion');
    end
    
    j               = 0:nf/2;
    wi              = 2*pi*(j*df + 0.01*df*(j==0));
    %         tic
    %         [vg,~,f1,ikmax,kx]=dispersion_curve_VP_wfix_Haskel_4inv_u2d3(para,wi);
    %         toc
%     tic
    [vg,~,f1,ikmax,kx]=dispersion_curve_VP_wfix_Haskel_4inv_d2u3(para,wi);
%     toc
    %se suprima el ultimo modo si tiene no mas un punto
    indpb           = find(ikmax==1);
    kx(indpb,:)     = [];
    ikmax(indpb,:)  = [];
    vg(indpb,:)     = [];
    f1(indpb,:)     = [];
    %se guarda los resultados en para.mode
    para.mode.k2    = kx;
    para.mode.ikmax = ikmax;
    para.mode.vg    = vg;
    para.mode.f1    = f1;
    
    %         %numero initial de ondas planas para la convergencia
    %         para.nincBWOP   = 200;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% calculo de las posiciones normalizadas de los receptores %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if para.siDesktop
    waitbarSwitch(0.25,h,'Identificacion de las propiedades de los receptores');
  else
    disp('Identificacion de las propiedades de los receptores');
  end
  para = pos_rec(para);
%   if para.geo(1)~=3 || (para.geo(1)==3 && para.nmed>1)
%     if para.siDesktop
%       if ~para.smallscreen
%         figure(1);plot(para.rec.xr,para.rec.zr,'.b');hold on
%       end
%     end
%   end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% init variable de salida %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if para.dim==1 %2D
    if para.pol==1 %SH
      ns  = sal.Ut+sal.USh+sal.UIh+sal.USt;
      uw  = zeros(nf/2+1,para.rec.nrec,para.ninc,ns);
      nss = sal.sxy+sal.syz;
      inds= [(sal.sxy==1)*1 (sal.syz==1)*2];
      sw  = zeros(nf/2+1,para.rec.nrec,para.ninc,nss);
    elseif para.pol==2 %PSV
      ns  = (sal.Ux + sal.Uz)*(sal.Ut+sal.UPh+sal.USh+sal.UIh+sal.UPt+sal.USt);
      uw  = zeros(nf/2+1,para.rec.nrec,para.ninc,ns);
      nss = sal.sxx+sal.szz+sal.sxz;
      inds= [(sal.sxx==1)*1 (sal.szz==1)*2 (sal.sxz==1)*3];
      sw  = zeros(nf/2+1,para.rec.nrec,para.ninc,nss);
    end
  else %2.5D,3Daxisim,3D
    ns = (sal.Ux + sal.Uy + sal.Uz)*(sal.Ut);
    uw = zeros(nf/2+1,para.rec.nrec,para.ninc,ns);
    nss= sal.sxx + sal.syy + sal.szz+sal.sxy+sal.sxz+sal.syz;
    sw = zeros(nf/2+1,para.rec.nrec,para.ninc,nss);
    %indice de los componentes del tensor de esfuerzos a guardar
    inds            = [(sal.sxx==1)*1 (sal.syy==1)*5 (sal.szz==1)*9 (sal.sxy==1)*2 (sal.sxz==1)*3 (sal.syz==1)*6];
  end
  inds(inds==0)   = [];
  DWN.inds        = inds;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% inicialisacion campos incidente DWN %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  para.nkmaxKW = 10000;
  if para.geo(1)==3 && para.nmed~=1 %solo cuando se ocupa el DWN en el IBEM
    % indepedientemente de la frecuencia, se alloca el espacio para el
    % campo de desplazamiento incidente en los receptores
    para.rec.m(1).nr = length(para.rec.m(1).ind); %cantidad de receptores
    if para.dim==1
      if para.pol==1 && para.nsubmed>1
        DWN.uy0     = zeros(para.rec.m(1).nr,para.ninc);
        DWN.s0      = zeros(nss,para.rec.m(1).nr,para.ninc);
      elseif para.pol==1 && para.nsubmed==1
        DWN.no      = 0;
      else%PSV
        DWN.uxz0  	= zeros(2,para.rec.m(1).nr,para.ninc);
        DWN.sxz0    = zeros(nss,para.rec.m(1).nr,para.ninc);
      end
    else
      DWN.U0          = zeros(3,para.rec.m(1).nr,para.ninc);
      DWN.S0          = zeros(nss,para.rec.m(1).nr,para.ninc);
    end
  elseif para.geo(1)==3 && para.nmed==1
    %solamente para dar una idea del espectro (kx,w) o (kr,w)
    %se initializa una matriz del espectro solo cuando hay pocos puntos
    %en k
    if para.DWNnbptkx/2+1<=para.nkmaxKW && para.dim==1
      UKW0=zeros(nf/2+1,para.DWNnbptkx/2+1);
    else
      %no se pintara el espectro (demasiado memoria)
      UKW0=zeros(nf/2+1,1);
    end
  else
    DWN=[];
  end
  
  % parte imaginaria de la frecuencia
  para.tmax = (para.zeropad-1)/(para.fmax/(para.nf/2)*para.zeropad);
  para.tmaxinteres = max(para.tmaxinteres,para.tmax);
  para.DWNomei = 1.0* pi/para.tmaxinteres; % puende ser 2.0 en lugar de 1.0
  DWNomei = para.DWNomei;
  
  % Paso en k
  %   máxima distancia horizontal entre fuente y receptores
%   xl = para.DWNxl;
%   pil = 2*pi/xl;
  
  
  if (~isfield(DWN,'U0') && para.geo(1)~=3)
    DWNomei = 0;
  else
    if para.siDesktop
      waitbarSwitch(0,h,['omega_i = ' num2str(DWNomei)]);
      waitbarSwitch(0,h,['tmax de interes = ' num2str(para.tmaxinteres)]);
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  dicretizacion de las superficies con revolucion (axisimetria)  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if     para.dim==3 && (para.nmed>1 || (para.nmed==1 && para.geo(1)==2))
    % caso considerado a partir del modelo 2D,
    % se ocupa no mas la discretisacion hacia el centro del objeto
    % y se hara una revolucion alrededor del eje z
    para = axi_contornos(para);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                                                                                          %
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
  %%                                 ciclo en frecuencias                                 % %
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
  %                                                                                          %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  geo     =para.geo(1);
  nmed    =para.nmed;
  pol     =para.pol;
  dim     =para.dim;
  
  para.usingparfor = true;
  if para.siDesktop
    if para.usingparfor
    waitbarSwitch(0,h,'Calculo de los espectros en paralelo ...');
    else
    waitbarSwitch(0,h,'Calculo de los espectros');
    end
  end
%   for j=0:nf/2
 parfor (j=0:nf/2,12)
    coord = j;    phi_fv = j;    DWNtmp = j;
    fj          = j*df + 0.01*df*(j==0) - 1i*DWNomei/2/pi; 

    %</ inicio de paratmp
    paratmp     = attenuation(para,fj);
    paratmp.j   = j;
    paratmp.fj  = fj;
    
    if ~paratmp.usingparfor
      str = sprintf('[ %d / %d]',round(j),round(nf/2));
      disp(str);
    end
     % kmax a cada frecuencia
    if geo == 3
    minbeta = 100000000000000000000000;
    for im = 1:paratmp.nsubmed
    minbeta = min(minbeta,paratmp.reg(1).sub(im).bet);
    end
%     paratmp.DWNkmax = 0.9*(2*pi*max(0.3*paratmp.fmax,real(fj)))/minbeta * 1.5;
    paratmp.DWNkmax = 0.9*(2*pi*max(0.15*paratmp.fmax,real(fj)))/minbeta * 1.5;
    end
    
    if nmed==1 && geo==1 
      %campo incidente solo en espacio completo
      if paratmp.fuente==2 && paratmp.meth_PS==1
        %calculo a partir de las correlaciones y teorema de
        %equiparticion
        [uw(j+1,:,:,:),sw(j+1,:,:,:)]	= Pure_correlation_equipartition(paratmp,fj);
      else
        [auw,asw] = campo_inc(paratmp);
        if max(max(max(auw)))~=0;  uw(j+1,:,:,:)=auw; end
        if max(max(max(asw)))~=0;  sw(j+1,:,:,:)=asw; end
%         clear auw asw
      end
    elseif nmed==1 && geo==3
      % campo incidente sólo en medio estratificado
      if paratmp.fuente==2 && paratmp.meth_PS==1
        %calculo a partir de las correlaciones y teorema de
        %equiparticion
        [uw(j+1,:,:,:),sw(j+1,:,:,:)]	= Pure_correlation_equipartition(paratmp,fj);
      else
        %Pure DWN
        [uw(j+1,:,:,:),UKW0(j+1,:),sw(j+1,:,:,:)]	= Pure_DWN(paratmp,fj);
      end
    else
      %IBEM
      %-------------------------------------%
      % calculo de las densidades de fuerza %
      %-------------------------------------%
      
      % calculo de las densidades de fuerzas "phi" para cada fuente virtual
      % la coordenades de cada fuerza esta dada en coord y estas cambian por
      % cada frecuencia
      if dim==1 %2D
        [phi_fv,coord,DWNtmp]  = sol_dfv(paratmp,fj,DWN);
      elseif dim>=3 %3D
        [phi_fv,coord,DWNtmp,paratmp]  = sol3D_dfv(paratmp,fj,DWN,j);
      end
      %         [phi_fv1,coord1]  = resample_phi_fv(coord,phi_fv,para,coord0);
      %         plot_phi_fv
      
      %--------------------------%
      % calculo en cada receptor %
      %--------------------------%
      %suma de las contribuciones de cada fuente virtual en cada receptor
      if dim==1 && pol==1
        %el cuarto ":" es por la paralelizacion pero en realidad el tamano es uw(j+1,:,:)
        [uw(j+1,:,:,:),sw(j+1,:,:,:)] 	= inversion_SH_k(paratmp,coord,phi_fv,gaussian,DWNtmp); %uy
      elseif dim==1 && pol==2
        %hay una ecuacion para phix y otra para phiz
        [uw(j+1,:,:,:),sw(j+1,:,:,:)]	= inversion_PSV_k(paratmp,coord,phi_fv,gaussian,DWNtmp); %u(1)=ux, u(2)=uz
      else
        %hay una ecuacion para cada uno de phix, phiy y phiz
        if para.dim == 3
        [uw(j+1,:,:,:),sw(j+1,:,:,:)] 	= inversion_3D_k(paratmp,coord,phi_fv,gaussex,DWNtmp); %u(1)=ux,u(2)=uy,u(3)=uz
        else
        [uw(j+1,:,:,:),sw(j+1,:,:,:)] 	= inversion_3D_k(paratmp,coord,phi_fv,gaussex,DWNtmp); %u(1)=ux,u(2)=uy,u(3)=uz
        end
      end
    end
    
    % /> fin de paratmp
    
    %         if toc(tstart)>1*60
%     Save_tmp_res(para,uw(j+1,:,:,:),j)
    %         end
  end % end for
  uw(isnan(uw))=0;
  sw(isnan(sw))=0;
  
  %dibujo espectro k-w del DWN
  if para.nmed==1 && para.geo(1)==3 && para.dim==1 % &&  para.DWNnbptkx/2+1<=para.nkmaxKW 
    %     if para.DWNnbptkx<=1000 && para.dim==1
    j           = 0:nf/2;
    wj          = 2*pi*(j*df + 0.01*df*(j==0));
    nk          = para.DWNnbptkx;
    nk0         = nk;
    DK          = 2*para.DWNkmax/(nk0);
    %xmax doit etre > vmax*tmax=alpha/dfe
    k2          = (0:(fix(nk0/2)))*DK;
    k2(1)       = k2(2)/1000;
    if para.siDesktop
      figure(837);surf(k2,wj,log(1+1e8*abs(UKW0)));shading flat;view([0 0 1]);
    end
    %         w0=40;
    %         indk2s=find(k2>w0/2,1,'first');
    %         hold on;plot3(k2(1:60:indk2s),w0,20,'ko')
    %         n00=length(k2(1:60:indk2s));
    %         dth=pi/2/(n00-1);
    %         kxmax=w0/2; kth=kxmax*sin((0:(n00-1))*dth);
    %         hold on;plot3(kth,w0+1,20,'k.')
    %     end
  end
else %if para.espyinv==0
  %% si hay que post-convolucionar unos espectros
  para0           = para;
  load(para.name,'para','uw','sw');
  uw(isnan(uw))   = 0; %#ok<NODEF>
  if exist('sw','var')
    sw(isnan(sw))   = 0; %#ok<NODEF>
  end
  para.pulso.tipo 	= para0.pulso.tipo ;
  para.pulso.a	= para0.pulso.a;
  para.pulso.b  = para0.pulso.b;
  para.pulso.c  = para0.pulso.c;
  para.spct     = para0.spct;
  para.siDesktop   	= para0.siDesktop;
end
%% variables de salida
name        = nametmp;
para.name   = name;

if isfield(para,'cont1')
  cont1 = para.cont1;
end
if exist('cont1','var')
  para.cont1  = cont1;
else
  cont1=[];
end
%% Guardar variable en frecuencia
RESULTADO.uw=uw;
RESULTADO.sw=sw;
save(name,'para','uw','sw','cont1');
  if para.siDesktop
    waitbarSwitch(0,h,'Resultados uw sw para (en frecuencia) guardados en archivo');
    waitbarSwitch(0,h,name);
  else
    disp('Resultados uw sw para (en frecuencia) guardados en archivo');
    dips(name);
  end
%%%%%%%%%%%%%%%%
%% inversion w %
%%%%%%%%%%%%%%%%
if para.spct==0
%   save([name,'tmp'],'para','uw','sw');%por si acaso
  
  if para.siDesktop
    waitbarSwitch(0,h,'Inversion de los espectros');
  else
    disp('Inversion de los espectros');
  end
  if para.zeropad < para.nf; para.zeropad = para.nf; end
  [utc,stc]   = inversion_w(uw,sw,para);
  
  if para.siDesktop
    waitbarSwitch(0,h,'Guardando los resultados ...');
  else
    disp('Guardando los resultados ...');
  end
else
  if para.siDesktop
    waitbarSwitch(0,h,'Guardando los resultados ...');
  else
    disp('Guardando los resultados ...');
  end
  utc         = 0;
  stc         = 0;
end
txttime     = conv_tps(tstart);
save([name(1:end-4) 't' name(end-4:end)],'para','utc','stc','cont1','txttime');
delete([name,'tmp*']);

%%%%%%%%%%%%%%%%%%%%%%
%% dibujo resultados %
%%%%%%%%%%%%%%%%%%%%%%
% if para.siDesktop
%   waitbarSwitch(0,h,'Dibujando los resultados ...');
% else
%   disp('no display of results available');
% end
% para.redraw     = 0;
% if para.siDesktop
%   para.b_dib(1).name   = name;
%   % b_dib           = dibujo(para,bdessin,utc,uw,stc,sw,cont1,b_dib);
%   close(h)
% end

%%%%%%%%%%%%%%%%%%%%%%
%% tiempo de computo %
%%%%%%%%%%%%%%%%%%%%%%
% if para.smallscreen
  if(ishandle(para.bar));set(para.bar,'BackgroundColor',[0 1 0],...
      'string',['calculation time: ',txttime]);end
% else
%   if para.siDesktop
%     msgbox(['calcul time: ',txttime]);
%   end
% end
disp(['calcul time: ',txttime]);

%% parpol
% cparll  = parpool('size');
% if cparll~=0
%     parpool close
% end

%% salida
RESULTADO.utc=utc;     clear utc
RESULTADO.uw=uw;       clear uw
RESULTADO.stc=stc;     clear stc
RESULTADO.sw=sw;       clear sw
RESULTADO.name=name;   clear name
RESULTADO.cont1=cont1; clear cont1
end
