para.rafraichi = 1;

if(ishandle(bouton.breptex))
  tmp=get(bouton.breptex,'string');
else
  tmp='';
end
if isfield(para,'nomrep')
if strcmp(para.nomrep,tmp)==0
    para.nomrep=tmp;
end
end

if (ishandle(bouton.dim)); set(bouton.dim      ,'value'    ,para.dim); 
else warning('Handles to UI are missing'); return; end
cmd_dim;

set(bouton.nmed     ,'string'   ,para.nmed);
set(bouton.pol      ,'value'    ,para.pol);
cmd_pol;

strmed =1:para.nmed;
set(bouton.med,'string',strmed);
if get(bouton.med,'value')>para.nmed
    set(bouton.med,'value',1);
end
cmd_med;
if med==1 && para.geo(1) == 3
    cmd_submed
    cmd_att_sub;
    cmd_lambda_mu_sub;
else
    cmd_att;
    cmd_lambda_mu;
end

set(bouton.geo      ,'value'    ,para.geo(get(bouton.med,'value'))');

if get(bouton.cont,'value')>2
   set(bouton.cont,'value',1);
end
cmd_cont;

set(bouton.fuente   ,'value'    ,para.fuente);
cmd_fuente;

set(bouton.ninc,'string',para.ninc);
if para.ninc~=length(para.strinc)
    para.strinc =[];
    for i=1:para.ninc
        tmpstr=num2str(i);
        while length(tmpstr)<floor(log10(para.ninc)+1)
            tmpstr=['0',tmpstr];
        end
        para.strinc =[para.strinc;tmpstr];
    end
end
if get(bouton.inc,'value')> para.ninc
    set(bouton.inc,'string',para.strinc,'value',1);
else
    set(bouton.inc,'string',para.strinc);
end
cmd_inc;

set(bouton.pulsotps ,'value'    ,para.pulso.tipo);
set(bouton.Ricker_tp,'string'   ,para.pulso.a);
set(bouton.delais   ,'string'   ,para.pulso.b);
set(bouton.GaussSta ,'string'   ,para.pulso.c);

%%%%%%%%%%%%%%%%%%%%%
% update receptores %
%%%%%%%%%%%%%%%%%%%%%
strrec	= 1:para.rec.nrecx;
if get(bouton.irec,'value')>para.rec.nrecx
    set(bouton.irec,'value',1);
end
set(bouton.irec,'string',strrec);
cmd_recpos;

%%%%%%%%%%%%%%%%%
% update sortie %
%%%%%%%%%%%%%%%%%
fn0	= fieldnames(para.sortie);
patchsortie;
n	= length(fn0);
for ifn0=1:n
    tmp=fn0(ifn0);
    nameb0 =['boutonsal.',tmp{1}];
    nameb2=['para.sortie.',tmp{1}];
    strcllbck=[nameb2,'=get(',nameb0,',''value'');'];
    set(eval(nameb0),'value',eval(nameb2));
end
if ~isfield(para,'spct')
    para.spct=0;
end
set(bouton.spct     ,'value',para.spct);
cmd_spec;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update parametros simulacion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(bouton.npplo    ,'string'   ,para.npplo);
set(bouton.fmax     ,'string'   ,para.fmax);
set(bouton.nf       ,'string'   ,para.nf);

%parametros DWN
% set(bouton.DWNkmax  ,'string',para.DWNkmax);
set(bouton.DWNnbptkx,'string',para.DWNnbptkx);
% set(bouton.DWNomei  ,'string',para.DWNomei);
set(bouton.tmax,'string',para.tmaxinteres);
% cmd_DWN;
clearfield;

dibujo_conf_geo(para,bouton.axe_conf_geo);

para.rafraichi = 0;

clear tmp onoff onoff1 onoff2 onoff3 onoff4 onoff5 
clear alpha bet lambda med mu rho strmed tatt icont i
clear strrec fn0 ifn0 n nameb0 nameb2 strcllbck DK DX 
clear strinfoDWN xmax